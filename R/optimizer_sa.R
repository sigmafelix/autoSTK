# optimizer_sa.R — Simulated Annealing optimiser for vgmST fitting.
#
# Strategy:
#   1. Use R's built-in optim(method = "SANN") to minimise the WLS criterion
#      over the parameter space [lower, upper].  SANN uses a Metropolis-type
#      acceptance criterion with a geometric cooling schedule.
#   2. After SANN converges, refine the solution with a short L-BFGS-B run so
#      the final model is at a proper local minimum.
#
# The objective is the same WLS used by .wls_at_par() in optimizer_grid.R,
# so both optimisers are directly comparable.

# Internal: map an unconstrained vector u to [lower, upper] via a scaled
# logistic function so that optim() can explore without hard boundaries.
# The inverse transform maps [lower, upper] back to ℝ.
.sigmoid_to_box <- function(u, lower, upper) {
  lower + (upper - lower) / (1 + exp(-u))
}
.box_to_sigmoid <- function(x, lower, upper) {
  p <- (x - lower) / (upper - lower)
  p <- pmax(pmin(p, 1 - 1e-8), 1e-8)   # clamp away from 0/1
  log(p / (1 - p))
}

# WLS objective in the unconstrained (sigmoid) parameterisation.
.sa_objective <- function(u, lower, upper, model_template, stva_emp) {
  par <- .sigmoid_to_box(u, lower, upper)
  .wls_at_par(par, model_template, stva_emp)
}

#' Simulated Annealing Optimiser for ST Variogram Fitting
#'
#' Uses \code{stats::optim(method = "SANN")} to globally search the WLS
#' objective surface for a \code{vgmST} model, then refines the best
#' solution with a short L-BFGS-B run.  Parameters are transformed to an
#' unconstrained space via a scaled logistic map so that SANN can explore
#' freely without explicit boundary enforcement.
#'
#' @param stva_emp Empirical ST variogram (\code{StVariogram}).
#' @param model_template Initial \code{vgmST} model (defines type/structure).
#' @param bounds List with \code{lower} and \code{upper} named numeric vectors
#'   (from \code{\link{st_param_bounds}}).
#' @param maxit Integer. Maximum number of SANN function evaluations.
#'   Default 5000.
#' @param temp Numeric. Initial temperature for the SANN cooling schedule.
#'   Default 10.
#' @param tmax Integer. Number of function evaluations at each temperature
#'   level before cooling.  Default 10.
#' @param lbfgsb_refine Logical. If \code{TRUE} (default), run a final short
#'   L-BFGS-B step from the SANN solution.
#' @param lbfgsb_maxit Integer. Maximum L-BFGS-B iterations for the
#'   refinement step.  Default 500.
#' @return Best-fitting \code{vgmST} object found.
#' @seealso \code{\link{autofitVariogramST}}, \code{\link{st_param_bounds}}
#' @keywords internal
.optimize_stv_sa <- function(stva_emp, model_template, bounds,
                             maxit         = 5000L,
                             temp          = 10,
                             tmax          = 10L,
                             lbfgsb_refine = TRUE,
                             lbfgsb_maxit  = 500L) {

  lower <- bounds$lower
  upper <- bounds$upper

  # Starting point: current model parameters, clamped to bounds
  init_par <- gstat::extractPar(model_template)
  init_par <- pmax(pmin(init_par, upper), lower)

  # Transform to unconstrained space for SANN
  u_init <- .box_to_sigmoid(init_par, lower, upper)

  # Build a closure so we don't need to pass lower/upper through optim's ...
  # (avoids name collision with optim()'s own lower/upper arguments)
  obj_fn <- function(u) {
    .sa_objective(u, lower, upper, model_template, stva_emp)
  }

  # Run SANN
  sann_result <- tryCatch(
    stats::optim(
      par     = u_init,
      fn      = obj_fn,
      method  = "SANN",
      control = list(
        maxit   = maxit,
        temp    = temp,
        tmax    = tmax,
        fnscale = 1       # minimise
      )
    ),
    error = function(e) {
      warning("SA optimiser (SANN) failed: ", conditionMessage(e),
              "\nReturning initial model.")
      NULL
    }
  )

  # If SANN failed, fall back to the model template
  if (is.null(sann_result)) {
    return(model_template)
  }

  # Recover best parameters
  best_par_u <- sann_result$par
  best_par   <- .sigmoid_to_box(best_par_u, lower, upper)

  mod_sa <- tryCatch(
    rlang::inject(gstat::vgmST(!!!best_par)),
    error = function(e) model_template
  )

  # ---- Optional L-BFGS-B refinement from SANN solution --------------------
  if (lbfgsb_refine) {
    mod_refined <- tryCatch(
      gstat::fit.StVariogram(
        object  = stva_emp,
        model   = mod_sa,
        method  = "L-BFGS-B",
        lower   = lower,
        upper   = upper,
        control = list(maxit = lbfgsb_maxit)
      ),
      error   = function(e) mod_sa,
      warning = function(w) {
        suppressWarnings(
          tryCatch(
            gstat::fit.StVariogram(
              object  = stva_emp,
              model   = mod_sa,
              method  = "L-BFGS-B",
              lower   = lower,
              upper   = upper,
              control = list(maxit = lbfgsb_maxit)
            ),
            error = function(e2) mod_sa
          )
        )
      }
    )

    # Accept refinement only if it improves the objective
    mse_sa      <- attr(mod_sa,      "MSErr")
    mse_refined <- attr(mod_refined, "MSErr")

    if (!is.null(mse_refined) && is.finite(mse_refined) &&
        (is.null(mse_sa) || !is.finite(mse_sa) || mse_refined < mse_sa)) {
      return(mod_refined)
    }
  }

  mod_sa
}
