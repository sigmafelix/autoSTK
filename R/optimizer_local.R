# optimizer_local.R — multi-start L-BFGS-B wrapper for vgmST fitting.

# Attempt to update a vgmST object's free parameters from a numeric vector.
# Uses gstat's internal updateVgmST if available, otherwise returns NULL.
.set_stvgm_par <- function(model_template, par) {
  tryCatch(
    gstat:::updateVgmST(model_template, par),
    error = function(e) NULL
  )
}

# Multi-start L-BFGS-B optimisation via gstat::fit.StVariogram.
#
# Arguments:
#   stva_emp       — empirical ST variogram
#   model_template — initial vgmST model
#   bounds         — list(lower, upper) from st_param_bounds()
#   n_restart      — number of optimisation runs (1 = current behaviour)
#   maxit          — maximum L-BFGS-B iterations
#
# Returns the best-fitting vgmST object (lowest MSErr).
.optimize_stv_lbfgsb <- function(stva_emp, model_template, bounds,
                                  n_restart = 1L, maxit = 2500L) {
  best_model <- NULL
  best_mserr <- Inf

  for (i in seq_len(n_restart)) {
    if (i == 1) {
      mod_start <- model_template
    } else {
      # Sample a new starting point from within bounds via LHS
      if (requireNamespace("lhs", quietly = TRUE)) {
        u <- lhs::randomLHS(1L, length(bounds$lower))[1L, ]
      } else {
        u <- runif(length(bounds$lower))
      }
      new_par   <- bounds$lower + u * (bounds$upper - bounds$lower)
      mod_start <- .set_stvgm_par(model_template, new_par)
      if (is.null(mod_start)) {
        # gstat internal not available; perturb by calling with the same start
        # (still useful because fit.StVariogram uses optim which has stochastic
        # elements in some methods, though L-BFGS-B is deterministic — so for
        # n_restart > 1 without updateVgmST we just run once)
        if (i > 1) next
        mod_start <- model_template
      }
    }

    # Clamp starting params inside bounds
    init_par <- extractPar(mod_start)
    init_par <- pmax(pmin(init_par, bounds$upper), bounds$lower)
    mod_start_clamped <- tryCatch(
      .set_stvgm_par(mod_start, init_par),
      error = function(e) mod_start
    )
    if (is.null(mod_start_clamped)) mod_start_clamped <- mod_start

    mod_fit <- tryCatch(
      fit.StVariogram(
        object  = stva_emp,
        model   = mod_start_clamped,
        method  = "L-BFGS-B",
        lower   = bounds$lower,
        upper   = bounds$upper,
        control = list(maxit = maxit)
      ),
      error   = function(e) NULL,
      warning = function(w) {
        suppressWarnings(
          tryCatch(
            fit.StVariogram(
              object  = stva_emp,
              model   = mod_start_clamped,
              method  = "L-BFGS-B",
              lower   = bounds$lower,
              upper   = bounds$upper,
              control = list(maxit = maxit)
            ),
            error = function(e2) NULL
          )
        )
      }
    )

    if (!is.null(mod_fit)) {
      mserr <- attr(mod_fit, "MSErr")
      if (!is.null(mserr) && is.finite(mserr) && mserr < best_mserr) {
        best_mserr <- mserr
        best_model <- mod_fit
      }
    }
  }

  if (is.null(best_model)) {
    warning("All optimisation attempts failed; returning initial model template.")
    best_model <- model_template
  }

  best_model
}
