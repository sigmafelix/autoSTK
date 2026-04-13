# test-estimate_initial_params.R
# Tests for the data-driven initial parameter estimator added in v2.1.0.

library(testthat)

# ---- Helpers -----------------------------------------------------------------

# Build a synthetic StVariogram-like data.frame with a known spherical shape.
# gamma(h, t) = sph_spatial(h) + sph_temporal(t)  (separable sum, no nugget)
# where sph(x; nugget, psill, range) is the standard spherical model.
.sph <- function(h, nugget, psill, range) {
  ifelse(h <= 0, 0,
    ifelse(h >= range,
      nugget + psill,
      nugget + psill * (1.5 * (h / range) - 0.5 * (h / range)^3)
    )
  )
}

make_synth_stv <- function(nugget    = 0.10,
                            psill     = 0.90,
                            sp_range  = 3000,
                            ts_range  = 3,
                            n_sp_bins = 12L,
                            n_t_bins  = 6L,
                            add_noise = FALSE,
                            seed      = 1L) {
  spl  <- seq(200, 8000, length.out = n_sp_bins)
  tl   <- seq(0, n_t_bins - 1L)

  grid <- expand.grid(spacelag = spl, timelag = tl)
  grid$gamma <- .sph(grid$spacelag, nugget, psill, sp_range) +
                .sph(as.numeric(grid$timelag) * 600, 0, 0.2, ts_range * 600)
  grid$np <- 30L

  if (add_noise) {
    set.seed(seed)
    grid$gamma <- abs(grid$gamma + rnorm(nrow(grid), 0, nugget * 0.2))
  }

  grid
}

# ---- Structure tests ---------------------------------------------------------

test_that("estimate_initial_params returns a named list with 8 required elements", {
  stv <- make_synth_stv()
  res <- estimate_initial_params(stv)

  expect_type(res, "list")
  expected_names <- c("nugget", "psill", "sill",
                       "sp_range", "ts_range",
                       "max_gamma", "max_splag", "max_tlag")
  expect_named(res, expected_names, ignore.order = TRUE)
})

test_that("all estimate_initial_params outputs are finite numerics", {
  stv <- make_synth_stv()
  res <- estimate_initial_params(stv)
  expect_true(all(vapply(res, is.numeric, logical(1L))))
  expect_true(all(vapply(res, is.finite,  logical(1L))))
})

test_that("nugget >= 0 and sill >= nugget", {
  stv <- make_synth_stv(nugget = 0.1, psill = 0.9)
  res <- estimate_initial_params(stv)
  expect_gte(res$nugget, 0)
  expect_gte(res$sill,   res$nugget)
})

test_that("psill = sill - nugget (internal consistency)", {
  stv <- make_synth_stv()
  res <- estimate_initial_params(stv)
  expect_equal(res$psill, res$sill - res$nugget, tolerance = 1e-9)
})

test_that("max_splag and max_tlag match data extent", {
  stv <- make_synth_stv()
  res <- estimate_initial_params(stv)
  expect_equal(res$max_splag, max(stv$spacelag))
  expect_equal(res$max_tlag,  max(as.numeric(stv$timelag)))
})

# ---- Accuracy on a known variogram -------------------------------------------

test_that("sp_range is within 2x of true spatial range for well-resolved variogram", {
  true_sp_range <- 3000
  stv <- make_synth_stv(sp_range = true_sp_range, n_sp_bins = 20L)
  res <- estimate_initial_params(stv)
  # The practical range at 95% of the sill should be within [0.5, 2] x true range
  expect_gte(res$sp_range, true_sp_range * 0.3)
  expect_lte(res$sp_range, true_sp_range * 2.0)
})

test_that("nugget estimate is non-negative and below the sill for nugget=0 data", {
  stv <- make_synth_stv(nugget = 0, psill = 1.0)
  res <- estimate_initial_params(stv)
  expect_gte(res$nugget, 0)
  expect_lte(res$nugget, res$sill)
})

test_that("sill estimate increases when true sill increases", {
  stv_lo <- make_synth_stv(psill = 0.5)
  stv_hi <- make_synth_stv(psill = 2.0)
  expect_gt(estimate_initial_params(stv_hi)$sill,
            estimate_initial_params(stv_lo)$sill)
})

# ---- Edge cases --------------------------------------------------------------

test_that("fallback when variogram never reaches target sill within window", {
  # Force a flat variogram (gamma = nugget everywhere) so the target is
  # never reached and the 0.5 * max_splag fallback activates.
  stv        <- make_synth_stv()
  stv$gamma  <- 0.1   # constant — psill estimated as ~0, target never crossed
  res <- estimate_initial_params(stv)
  # Fallback branch: sp_range = 0.5 * max_splag
  expect_equal(res$sp_range, max(stv$spacelag) * 0.5, tolerance = 1e-6)
})

test_that("single time lag is handled without error", {
  stv <- make_synth_stv(n_t_bins = 1L)
  expect_no_error(res <- estimate_initial_params(stv))
  expect_gte(res$nugget, 0)
})

test_that("single spatial bin is handled without error", {
  stv <- make_synth_stv(n_sp_bins = 1L)
  expect_no_error(res <- estimate_initial_params(stv))
})

test_that("noisy variogram still returns valid estimates", {
  stv <- make_synth_stv(nugget = 0.15, psill = 0.85, add_noise = TRUE)
  res <- estimate_initial_params(stv)
  expect_gte(res$nugget, 0)
  expect_gte(res$sill,   res$nugget)
  expect_true(all(vapply(res, is.finite, logical(1L))))
})

# ---- Missing / bad input -----------------------------------------------------

test_that("error on missing required column 'gamma'", {
  stv <- make_synth_stv()
  stv$gamma <- NULL
  expect_error(estimate_initial_params(stv), "gamma")
})

test_that("error on missing required column 'spacelag'", {
  stv <- make_synth_stv()
  stv$spacelag <- NULL
  expect_error(estimate_initial_params(stv), "spacelag")
})

test_that("NA values in gamma are silently ignored", {
  stv <- make_synth_stv()
  stv$gamma[c(1L, 5L, 10L)] <- NA_real_
  expect_no_error(res <- estimate_initial_params(stv))
  expect_true(all(vapply(res, is.finite, logical(1L))))
})

test_that("all-NA gamma triggers an error", {
  stv        <- make_synth_stv()
  stv$gamma  <- NA_real_
  expect_error(estimate_initial_params(stv), "No valid")
})

# ---- Custom quantile arguments -----------------------------------------------

test_that("higher sill_quantile yields a >= estimate vs lower quantile", {
  stv  <- make_synth_stv(add_noise = TRUE)
  lo   <- estimate_initial_params(stv, sill_quantile = 0.75)
  hi   <- estimate_initial_params(stv, sill_quantile = 0.99)
  expect_gte(hi$sill, lo$sill)
})

test_that("practical_range_target = 0.5 gives smaller sp_range than 0.95", {
  stv    <- make_synth_stv(sp_range = 3000, n_sp_bins = 20L)
  small  <- estimate_initial_params(stv, practical_range_target = 0.50)
  large  <- estimate_initial_params(stv, practical_range_target = 0.95)
  expect_lte(small$sp_range, large$sp_range)
})
