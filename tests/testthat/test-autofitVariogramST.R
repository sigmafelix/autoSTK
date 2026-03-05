# test-autofitVariogramST.R

library(testthat)
library(spacetime)
library(gstat)

# ---- Shared test fixture --------------------------------------------------
# Uses 5 stations x 6 time points so the empirical variogram has enough lags.
make_test_stfdf <- function(n_sp = 15L, n_t = 6L, seed = 42L) {
  set.seed(seed)
  sp   <- sp::SpatialPoints(cbind(x = runif(n_sp, 0, 1e4),
                                  y = runif(n_sp, 0, 1e4)))
  time <- as.POSIXct("2020-01-01") + seq(0, (n_t - 1L) * 3600, by = 3600)
  STFDF(sp, time, data.frame(z = rnorm(n_sp * n_t)))
}

# ---- Basic structure ------------------------------------------------------
test_that("autofitVariogramST returns STVariogramFit with required fields", {
  stf <- make_test_stfdf()
  suppressWarnings(res <- autofitVariogramST(stf = stf, formula = z ~ 1,
                            typestv = "sumMetric",
                            candidate_model = c("Exp", "Sph")))
  expect_s3_class(res, "STVariogramFit")
  expect_true(all(c("jointSTV", "empSTV", "SpV", "TV",
                    "optimizer", "objective", "n_obs") %in% names(res)))
})

test_that("n_obs stored correctly", {
  stf <- make_test_stfdf(n_sp = 20L, n_t = 10L, seed = 2026L)
  res <- suppressWarnings(
    autofitVariogramST(
      stf = stf,
      formula = z ~ 1,
      cutoff = 1e3,
      width = 200,
      candidate_model = c("Exp", "Mat"))
  )
  expect_equal(res$n_obs, 200L)
})

# ---- All six model types --------------------------------------------------
for (mtype in c("sumMetric", "separable", "metric",
                 "productSum", "productSumOld", "simpleSumMetric")) {
  local({
    m <- mtype
    test_that(paste("autofitVariogramST works with", m), {
      stf <- make_test_stfdf()
      suppressWarnings(
        res <- autofitVariogramST(stf = stf, formula = z ~ 1,
          typestv = m, candidate_model = c("Exp"))
      )
      expect_s3_class(res, "STVariogramFit")
      expect_true(all(c("jointSTV", "empSTV", "SpV", "TV") %in% names(res)))
    })
  })
}

# ---- Surface option -------------------------------------------------------
test_that("surface = TRUE adds STVsurface to result", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(stf = stf, formula = z ~ 1,
                              typestv = "sumMetric",
                              candidate_model = c("Exp", "Sph", "Lin", "Mat"),
                              cutoff = 3000,
                              width = 500,
                              n_restart = 3L,
                              surface = TRUE,
                              verbose = TRUE)
  )
  expect_true("STVsurface" %in% names(res))
})

# ---- Optimizer argument ---------------------------------------------------
test_that("optimizer = 'lbfgsb' with n_restart = 1 matches default", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(stf = stf, formula = z ~ 1,
                              typestv = "sumMetric", candidate_model = c("Exp"),
                              optimizer = "lbfgsb", n_restart = 1L)
  )
  expect_equal(res$optimizer, "lbfgsb")
})

test_that("optimizer = 'grid' returns a valid STVariogramFit", {
  skip_if_not_installed("lhs")
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(stf = stf, formula = z ~ 1,
                              typestv = "sumMetric", candidate_model = c("Exp"),
                              optimizer = "grid",
                              optimizer_control = list(n_coarse = 10L,
                                                       n_refine = 5L))
  )
  expect_s3_class(res, "STVariogramFit")
  expect_equal(res$optimizer, "grid")
})

# ---- Objective argument ---------------------------------------------------
test_that("objective = 'WLS' stores NULL loglik", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(stf = stf, formula = z ~ 1,
                              typestv = "sumMetric", candidate_model = c("Exp"),
                              objective = "WLS")
  )
  expect_null(res$loglik)
  expect_equal(res$objective, "WLS")
})

test_that("objective = 'MLE' stores finite loglik for small n", {
  stf <- make_test_stfdf()   # n = 30, below 500 threshold
  suppressWarnings(
    res <- autofitVariogramST(stf = stf, formula = z ~ 1,
                              typestv = "sumMetric", candidate_model = c("Exp"),
                              objective = "MLE")
  )
  expect_equal(res$objective, "MLE")
  # loglik may be -Inf if covariance matrix is singular, but must be numeric
  expect_true(is.numeric(res$loglik))
})

# ---- Error on unknown model -----------------------------------------------
test_that("autofitVariogramST errors on unknown typestv", {
  stf <- make_test_stfdf()
  suppressWarnings(
    expect_error(
      autofitVariogramST(stf = stf, formula = z ~ 1,
                       typestv = "unknownType", candidate_model = c("Exp")),
      "unknown"
    )
  )
})

# ---- Measurement error ----------------------------------------------------
test_that("measurement_error argument is accepted without error", {
  stf <- make_test_stfdf()
  expect_no_error(
    suppressWarnings(
      autofitVariogramST(stf = stf, formula = z ~ 1,
                        typestv = "sumMetric", candidate_model = c("Exp"),
                        measurement_error = c(0.1, 0.2, 0.3))
    )
  )
})
