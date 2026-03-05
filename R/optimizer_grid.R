# optimizer_grid.R — LHS grid search followed by L-BFGS-B refinement.
#
# Strategy:
#   1. Sample n_coarse points from [lower, upper] via Latin Hypercube Sampling.
#   2. Evaluate WLS at each point using variogramSurface (cheap; no optim).
#   3. Shrink bounds to a neighbourhood of the best point.
#   4. Sample n_refine points from the refined bounds.
#   5. Take the overall best point as the starting value for fit.StVariogram.

# Evaluate the WLS criterion for a given parameter vector directly, without
# calling optim (fast, for grid evaluation).
.wls_at_par <- function(par, model_template, stva_emp) {
  mod <- tryCatch(
    gstat:::updateVgmST(model_template, par),
    error = function(e) NULL
  )
  if (is.null(mod)) return(Inf)

  pred_gamma <- tryCatch(
    variogramSurface(mod, stva_emp[, c("timelag", "spacelag")])$gamma,
    error = function(e) rep(NA_real_, nrow(stva_emp))
  )

  if (any(!is.finite(pred_gamma))) return(Inf)

  # WLS: weight = np / gamma_model^2  (same as gstat internals)
  w <- as.numeric(stva_emp$np) / pmax(pred_gamma^2, 1e-12)
  sum(w * (stva_emp$gamma - pred_gamma)^2, na.rm = TRUE)
}

# LHS grid search + refinement.
#
# Arguments:
#   stva_emp       — empirical ST variogram
#   model_template — initial vgmST model (defines model type and structure)
#   bounds         — list(lower, upper) from st_param_bounds()
#   n_coarse       — number of LHS points in the global phase
#   n_refine       — number of LHS points in the local refinement phase
#   refine_frac    — fraction of [lower,upper] range used as refinement radius
#   maxit          — max iterations for the final L-BFGS-B call
#   seed           — random seed for reproducibility
#
# Returns the best-fitting vgmST object.
.optimize_stv_grid <- function(stva_emp, model_template, bounds,
                                n_coarse    = 50L,
                                n_refine    = 30L,
                                refine_frac = 0.25,
                                maxit       = 2500L,
                                seed        = 42L) {
  set.seed(seed)

  k <- length(bounds$lower)

  # ---- Phase 1: coarse LHS ------------------------------------------------
  lhs_coarse  <- lhs::randomLHS(n_coarse, k)
  par_coarse  <- sweep(
    sweep(lhs_coarse, 2, bounds$upper - bounds$lower, "*"),
    2, bounds$lower, "+"
  )

  scores_coarse <- apply(par_coarse, 1L, function(p) {
    .wls_at_par(p, model_template, stva_emp)
  })

  best_coarse_idx <- which.min(scores_coarse)
  best_par        <- par_coarse[best_coarse_idx, ]

  # ---- Phase 2: local refinement ------------------------------------------
  radius    <- refine_frac * (bounds$upper - bounds$lower)
  ref_lower <- pmax(bounds$lower, best_par - radius)
  ref_upper <- pmin(bounds$upper, best_par + radius)
  # Ensure ref_lower < ref_upper
  ref_upper <- pmax(ref_upper, ref_lower + 1e-6)

  lhs_refine  <- lhs::randomLHS(n_refine, k)
  par_refine  <- sweep(
    sweep(lhs_refine, 2, ref_upper - ref_lower, "*"),
    2, ref_lower, "+"
  )

  scores_refine <- apply(par_refine, 1L, function(p) {
    .wls_at_par(p, model_template, stva_emp)
  })

  all_scores <- c(scores_coarse, scores_refine)
  all_pars   <- rbind(par_coarse, par_refine)
  best_par   <- all_pars[which.min(all_scores), ]

  # ---- Phase 3: L-BFGS-B from best grid point ----------------------------
  mod_start <- tryCatch(
    gstat:::updateVgmST(model_template, best_par),
    error = function(e) model_template
  )

  mod_fit <- tryCatch(
    fit.StVariogram(
      object  = stva_emp,
      model   = mod_start,
      method  = "L-BFGS-B",
      lower   = bounds$lower,
      upper   = bounds$upper,
      control = list(maxit = maxit)
    ),
    error = function(e) {
      warning("Grid optimizer: L-BFGS-B refinement failed. Returning grid best.")
      mod_start
    }
  )

  mod_fit
}
