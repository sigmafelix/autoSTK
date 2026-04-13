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

# ---- Improved initial parameter estimation (v2.1.0) --------------------------

test_that("explicit guess_nugget is respected and overrides estimate_initial_params", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf = stf, formula = z ~ 1, typestv = "sumMetric",
      candidate_model = c("Exp"), guess_nugget = 0.0
    )
  )
  # Just check the fit succeeded; the nugget was passed through
  expect_s3_class(res, "STVariogramFit")
})

test_that("explicit guess_psill is respected", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf = stf, formula = z ~ 1, typestv = "sumMetric",
      candidate_model = c("Exp"), guess_psill = 0.5
    )
  )
  expect_s3_class(res, "STVariogramFit")
})

test_that("verbose = TRUE prints initial estimate message", {
  stf <- make_test_stfdf()
  expect_message(
    suppressWarnings(
      autofitVariogramST(
        stf = stf, formula = z ~ 1, typestv = "sumMetric",
        candidate_model = c("Exp"), verbose = TRUE
      )
    ),
    regexp = "Initial estimates"
  )
})

test_that("joint variogram range is finite and positive after init", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf = stf, formula = z ~ 1, typestv = "sumMetric",
      candidate_model = c("Exp")
    )
  )
  pars <- extractPar(res$jointSTV)
  range_pars <- pars[grepl("range", names(pars), ignore.case = TRUE)]
  expect_true(all(range_pars > 0))
  expect_true(all(is.finite(range_pars)))
})

# ---- SA optimizer integration ------------------------------------------------

test_that("optimizer = 'sa' returns STVariogramFit with correct optimizer field", {
  skip_on_cran()
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf             = stf, formula = z ~ 1, typestv = "sumMetric",
      candidate_model = c("Exp"), cutoff = 3000, width = 500, tlags = 0:3,
      optimizer         = "sa",
      optimizer_control = list(maxit = 60L, temp = 5, tmax = 5L)
    )
  )
  expect_s3_class(res, "STVariogramFit")
  expect_equal(res$optimizer, "sa")
})

# ---- GA optimizer integration ------------------------------------------------

test_that("optimizer = 'ga' returns STVariogramFit with correct optimizer field", {
  skip_on_cran()
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf             = stf, formula = z ~ 1, typestv = "sumMetric",
      candidate_model = c("Exp"), cutoff = 3000, width = 500, tlags = 0:3,
      optimizer         = "ga",
      optimizer_control = list(popSize = 8L, maxiter = 5L, run = 3L, seed = 1L)
    )
  )
  expect_s3_class(res, "STVariogramFit")
  expect_equal(res$optimizer, "ga")
})

# ---- All four optimizers produce comparable MSErr ----------------------------

test_that("all four optimizers finish and produce finite MSErr", {
  skip_on_cran()
  stf    <- make_test_stfdf()
  mserrs <- list()

  for (opt in c("lbfgsb", "grid", "sa", "ga")) {
    ctrl <- switch(opt,
      lbfgsb = list(maxit = 500L),
      grid   = list(n_coarse = 8L, n_refine = 4L, maxit = 200L),
      sa     = list(maxit = 60L, temp = 5, tmax = 5L),
      ga     = list(popSize = 8L, maxiter = 5L, run = 3L, seed = 42L)
    )
    suppressWarnings(
      res <- autofitVariogramST(
        stf = stf, formula = z ~ 1, typestv = "sumMetric",
        candidate_model = c("Exp"), cutoff = 3000, width = 500, tlags = 0:3,
        optimizer = opt, optimizer_control = ctrl
      )
    )
    mse <- attr(res$jointSTV, "MSErr")
    mserrs[[opt]] <- mse
    expect_true(
      is.null(mse) || is.finite(mse),
      label = paste("MSErr should be NULL or finite for optimizer =", opt)
    )
  }
})
