# test-autoKrigeST_cv.R

library(testthat)
library(spacetime)

# Minimal STFDF with enough stations and time steps for CV
make_cv_stfdf <- function(n_sp = 15L, n_t = 6L, seed = 99L) {
  set.seed(seed)
  sp   <- sp::SpatialPoints(cbind(x = runif(n_sp, 0, 1e5),
                                  y = runif(n_sp, 0, 1e5)))
  sp::proj4string(sp) <- sp::CRS("+proj=utm +zone=32 +datum=WGS84")
  time <- as.POSIXct("2020-01-01") + 0:(n_t - 1L) * 86400
  STFDF(sp, time, data.frame(z = rnorm(n_sp * n_t, mean = 10, sd = 2)))
}

# ---- Spatial fold ---------------------------------------------------------
test_that("autoKrigeST.cv spatial fold returns data.frame with RMSE/MAE/BIAS", {
  skip_on_cran()
  stf <- make_cv_stfdf()
  suppressWarnings(
    res <- autoKrigeST.cv(
      data     = stf, formula = z ~ 1, fold_dim = "spatial",
      nfold    = 3L, cutoff = 6e4, width = 1e4, tlags = 0:3, cores = 1L
    )
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("CVFold", "RMSE", "MAE", "BIAS") %in% names(res)))
  expect_equal(nrow(res), 3L)
  expect_true(all(unlist(res[, "RMSE"]) >= 0))
})

# ---- Temporal fold --------------------------------------------------------
test_that("autoKrigeST.cv temporal fold works", {
  skip_on_cran()
  stf <- make_cv_stfdf()
  suppressWarnings(
    res <- autoKrigeST.cv(
      data     = stf, formula = z ~ 1, fold_dim = "temporal",
      nfold    = 3L, cutoff = 6e4, width = 1e4, tlags = 0:3, cores = 1L
    )
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3L)
})

# ---- variogram_from_full --------------------------------------------------
test_that("variogram_from_full = TRUE produces same structure as FALSE", {
  skip_on_cran()
  stf <- make_cv_stfdf()
  suppressWarnings(
    res_full <- autoKrigeST.cv(
      data = stf, formula = z ~ 1, fold_dim = "spatial",
      nfold = 3L, cutoff = 6e4, width = 1e4, tlags = 0:3, cores = 1L,
      variogram_from_full = TRUE
    )
  )
  expect_s3_class(res_full, "data.frame")
  expect_true(all(c("CVFold", "RMSE", "MAE", "BIAS") %in% names(res_full)))
})

# ---- spacetime fold index bug regression ---------------------------------
# Previously `v[v_tindex]` referenced undefined `v`; this test exercises
# the spacetime fold path with nfold = 4 (sqrt = 2, integer)
test_that("autoKrigeST.cv spacetime fold no longer crashes with nfold = 4", {
  skip_on_cran()
  stf <- make_cv_stfdf(n_sp = 12L, n_t = 8L)
  expect_no_error(
    suppressWarnings(
      autoKrigeST.cv(
        data     = stf, formula = z ~ 1, fold_dim = "spacetime",
        nfold    = 4L, cutoff = 6e4, width = 1e4, tlags = 0:3, cores = 1L
      )
    )
  )
})
