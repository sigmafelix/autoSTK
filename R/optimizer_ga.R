# optimizer_ga.R — Genetic Algorithm optimiser for vgmST fitting.
#
# Primary path: uses the GA package (Luca Scrucca, 2013) which provides a
# real-valued GA with tournament selection, BLX-alpha crossover, and uniform
# mutation.  The WLS criterion is negated because GA::ga() *maximises*.
#
# Fallback: a lightweight Differential Evolution (DE/rand/1/bin) that runs
# entirely in base R when the GA package is unavailable.  DE is functionally
# similar in exploration breadth while using zero external dependencies.
#
# In both cases the best solution found by the stochastic search is handed
# to a short L-BFGS-B run for local refinement.


# ---- Shared fitness wrapper -----------------------------------------------

# GA::ga() maximises, so we negate the WLS criterion.
.ga_fitness <- function(par, model_template, stva_emp) {
  -.wls_at_par(par, model_template, stva_emp)
}


# ---- Fallback: Differential Evolution (DE/rand/1/bin) ---------------------

# A minimal DE implementation that only depends on base R.
# Population size: pop_size, generations: maxiter, early stop: run
# de_f = 0.8 (scale factor), de_cr = 0.9 (crossover rate)
.de_optimizer <- function(model_template, stva_emp, lower, upper,
                          pop_size = 50L, maxiter = 200L, run = 30L,
                          seed = 42L) {
  set.seed(seed)
  k     <- length(lower)
  de_f  <- 0.8   # DE scale factor
  de_cr <- 0.9   # crossover rate

  # Initialise population with LHS if available, otherwise uniform random
  if (requireNamespace("lhs", quietly = TRUE)) {
    pop <- lhs::randomLHS(pop_size, k)
  } else {
    pop <- matrix(stats::runif(pop_size * k), nrow = pop_size, ncol = k)
  }
  pop <- sweep(sweep(pop, 2, upper - lower, "*"), 2, lower, "+")

  # Evaluate initial fitness (negative WLS = fitness to maximise)
  fitness <- apply(pop, 1L,
    function(p) .ga_fitness(p, model_template, stva_emp))

  best_idx   <- which.max(fitness)
  best_par   <- pop[best_idx, ]
  best_fit   <- fitness[best_idx]
  no_improve <- 0L

  for (gen in seq_len(maxiter)) {
    for (i in seq_len(pop_size)) {
      # Select three distinct individuals != i
      candidates <- setdiff(seq_len(pop_size), i)
      r          <- sample(candidates, 3L)
      vec_a <- pop[r[1L], ]
      vec_b <- pop[r[2L], ]
      vec_c <- pop[r[3L], ]

      # Mutation: donor vector
      donor <- vec_a + de_f * (vec_b - vec_c)
      # Clamp to bounds
      donor <- pmax(pmin(donor, upper), lower)

      # Crossover: binomial
      cross_mask <- stats::runif(k) < de_cr
      if (!any(cross_mask)) {
        cross_mask[sample.int(k, 1L)] <- TRUE
      }
      trial <- ifelse(cross_mask, donor, pop[i, ])

      # Selection
      trial_fit <- .ga_fitness(trial, model_template, stva_emp)
      if (trial_fit > fitness[i]) {
        pop[i, ]   <- trial
        fitness[i] <- trial_fit
      }
    }

    # Track global best
    gen_best_idx <- which.max(fitness)
    if (fitness[gen_best_idx] > best_fit) {
      best_fit   <- fitness[gen_best_idx]
      best_par   <- pop[gen_best_idx, ]
      no_improve <- 0L
    } else {
      no_improve <- no_improve + 1L
    }

    if (no_improve >= run) break
  }

  best_par
}


# ---- Main GA optimiser -------------------------------------------------------

#' Genetic Algorithm Optimiser for ST Variogram Fitting
#'
#' Uses the \pkg{GA} package's real-valued genetic algorithm to globally
#' search the WLS objective surface, then refines the solution with a short
#' L-BFGS-B run.  Falls back to a built-in Differential Evolution
#' (\code{DE/rand/1/bin}) implementation when \pkg{GA} is not installed.
#'
#' @param stva_emp Empirical ST variogram (\code{StVariogram}).
#' @param model_template Initial \code{vgmST} model (defines type/structure).
#' @param bounds List with \code{lower} and \code{upper} named numeric vectors
#'   (from \code{\link{st_param_bounds}}).
#' @param popSize Integer. Population size.  Default 50.
#' @param maxiter Integer. Maximum number of generations.  Default 200.
#' @param run Integer. Stop after this many generations without improvement.
#'   Default 30.
#' @param seed Integer. Random seed for reproducibility.  Default 42.
#' @param lbfgsb_refine Logical. Refine the GA solution with L-BFGS-B.
#'   Default \code{TRUE}.
#' @param lbfgsb_maxit Integer. Maximum L-BFGS-B iterations for refinement.
#'   Default 500.
#' @return Best-fitting \code{vgmST} object found.
#' @seealso \code{\link{autofitVariogramST}}, \code{\link{st_param_bounds}}
#' @keywords internal
.optimize_stv_ga <- function(stva_emp, model_template, bounds,
                             popSize       = 50L,
                             maxiter       = 200L,
                             run           = 30L,
                             seed          = 42L,
                             lbfgsb_refine = TRUE,
                             lbfgsb_maxit  = 500L) {

  lower <- bounds$lower
  upper <- bounds$upper

  # ---- Global search -------------------------------------------------------
  if (requireNamespace("GA", quietly = TRUE)) {
    # Use GA::ga() — real-valued mode
    ga_result <- tryCatch(
      GA::ga(
        type    = "real-valued",
        fitness = function(par) .ga_fitness(par, model_template, stva_emp),
        lower   = lower,
        upper   = upper,
        popSize = popSize,
        maxiter = maxiter,
        run     = run,
        seed    = seed,
        monitor = FALSE
      ),
      error = function(e) {
        warning("GA optimiser (GA package) failed: ", conditionMessage(e),
                "\nFalling back to Differential Evolution.")
        NULL
      }
    )

    if (!is.null(ga_result)) {
      # GA::ga stores the best solution in @solution (one row per tie)
      best_par <- as.numeric(ga_result@solution[1L, ])
    } else {
      best_par <- .de_optimizer(model_template, stva_emp,
                                lower, upper, popSize, maxiter, run, seed)
    }
  } else {
    message("Package 'GA' not available; using built-in Differential Evolution.")
    best_par <- .de_optimizer(model_template, stva_emp,
                              lower, upper, popSize, maxiter, run, seed)
  }

  # ---- Reconstruct vgmST from best parameters -----------------------------
  mod_ga <- tryCatch(
    gstat:::updateVgmST(model_template, best_par),
    error = function(e) {
      warning("Could not update vgmST with GA solution; returning template.")
      model_template
    }
  )

  # ---- Optional L-BFGS-B refinement ---------------------------------------
  if (lbfgsb_refine) {
    mod_refined <- tryCatch(
      fit.StVariogram(
        object  = stva_emp,
        model   = mod_ga,
        method  = "L-BFGS-B",
        lower   = lower,
        upper   = upper,
        control = list(maxit = lbfgsb_maxit)
      ),
      error   = function(e) mod_ga,
      warning = function(w) {
        suppressWarnings(
          tryCatch(
            fit.StVariogram(
              object  = stva_emp,
              model   = mod_ga,
              method  = "L-BFGS-B",
              lower   = lower,
              upper   = upper,
              control = list(maxit = lbfgsb_maxit)
            ),
            error = function(e2) mod_ga
          )
        )
      }
    )

    mse_ga      <- attr(mod_ga,      "MSErr")
    mse_refined <- attr(mod_refined, "MSErr")

    if (!is.null(mse_refined) && is.finite(mse_refined) &&
        (is.null(mse_ga) || !is.finite(mse_ga) || mse_refined < mse_ga)) {
      return(mod_refined)
    }
  }

  mod_ga
}
