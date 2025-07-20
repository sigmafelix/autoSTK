# test_autofitVariogramST.r

library(testthat)
library(spacetime)
library(gstat)

# Helper: create minimal spatio-temporal data
make_test_stfdf <- function() {
  # 3 spatial points
  sp <- sp::SpatialPoints(cbind(
    x = c(0, 1, 2),
    y = c(0, 1, 2)
  ))
  # 3 time points
  time <- as.POSIXct("2020-01-01") + 0:2 * 3600
  # expand.grid for all combinations
  STFDF(sp, time, data.frame(z = rnorm(9)))
}

test_that("autofitVariogramST returns expected structure for sumMetric", {
  stf <- make_test_stfdf()
  res <- autofitVariogramST(
    stf = stf,
    formula = z ~ 1,
    typestv = "sumMetric",
    candidate_model = c("Exp", "Sph"),
    surface = FALSE
  )
  expect_type(res, "list")
  expect_true(all(c("jointSTV", "empSTV", "SpV", "TV") %in% names(res)))
})

test_that("autofitVariogramST works with separable model", {
  stf <- make_test_stfdf()
  res <- autofitVariogramST(
    stf = stf,
    formula = z ~ 1,
    typestv = "separable",
    candidate_model = c("Exp"),
    surface = FALSE
  )
  expect_type(res, "list")
  expect_true(all(c("jointSTV", "empSTV", "SpV", "TV") %in% names(res)))
})

test_that("autofitVariogramST returns surface when requested", {
  stf <- make_test_stfdf()
  res <- autofitVariogramST(
    stf = stf,
    formula = z ~ 1,
    typestv = "sumMetric",
    candidate_model = c("Exp"),
    surface = TRUE
  )
  expect_true("STVsurface" %in% names(res))
})

test_that("autofitVariogramST handles metric type", {
  stf <- make_test_stfdf()
  res <- autofitVariogramST(
    stf = stf,
    formula = z ~ 1,
    typestv = "metric",
    candidate_model = c("Exp"),
    surface = FALSE
  )
  expect_type(res, "list")
  expect_true(all(c("jointSTV", "empSTV", "SpV", "TV") %in% names(res)))
})

test_that("autofitVariogramST works with custom measurement_error", {
  stf <- make_test_stfdf()
  res <- autofitVariogramST(
    stf = stf,
    formula = z ~ 1,
    typestv = "sumMetric",
    candidate_model = c("Exp"),
    measurement_error = c(0.1, 0.2, 0.3),
    surface = FALSE
  )
  expect_type(res, "list")
})

test_that("autofitVariogramST works with productSum", {
  stf <- make_test_stfdf()
  res <- autofitVariogramST(
    stf = stf,
    formula = z ~ 1,
    typestv = "productSum",
    candidate_model = c("Exp"),
    surface = FALSE
  )
  expect_type(res, "list")
})

test_that("autofitVariogramST works with productSumOld", {
  stf <- make_test_stfdf()
  res <- autofitVariogramST(
    stf = stf,
    formula = z ~ 1,
    typestv = "productSumOld",
    candidate_model = c("Exp"),
    surface = FALSE
  )
  expect_type(res, "list")
})

test_that("autofitVariogramST works with simpleSumMetric", {
  stf <- make_test_stfdf()
  res <- autofitVariogramST(
    stf = stf,
    formula = z ~ 1,
    typestv = "simpleSumMetric",
    candidate_model = c("Exp"),
    surface = FALSE
  )
  expect_type(res, "list")
})

# Optionally, test error handling for unknown model
test_that("autofitVariogramST errors on unknown typestv", {
  stf <- make_test_stfdf()
  expect_error(
    autofitVariogramST(
      stf = stf,
      formula = z ~ 1,
      typestv = "unknownType",
      candidate_model = c("Exp"),
      surface = FALSE
    ),
    "unknown"
  )
})
