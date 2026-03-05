# test-generics.R

library(testthat)
library(spacetime)

make_test_stfdf <- function(seed = 7L) {
  set.seed(seed)
  sp   <- sp::SpatialPoints(cbind(x = runif(15L, 0, 1e4),
                                  y = runif(15L, 0, 1e4)))
  time <- as.POSIXct("2020-01-01") + 0:5 * 3600
  STFDF(sp, time, data.frame(z = rnorm(90L)))
}

test_that("print.STVariogramFit runs without error", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf = stf, formula = z ~ 1,
      cutoff = 3e3, width = 500,
      typestv = "sumMetric",
      candidate_model = c("Exp", "Mat", "Wav")
    )
  )
  expect_no_error(print(res))
  expect_invisible(print(res))
})

test_that("summary.STVariogramFit runs without error", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf = stf, formula = z ~ 1,
      cutoff = 3e3, width = 500,
      typestv = "sumMetric", candidate_model = c("Exp", "Mat", "Lin")
    )
  )
  expect_no_error(summary(res))
})

test_that("AIC.STVariogramFit returns NA with warning for WLS fit", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf = stf, formula = z ~ 1,
      cutoff = 3e3, width = 500,
      typestv = "sumMetric", candidate_model = c("Mat", "Exc", "Lin"),
      objective = "WLS"
    )
  )
  expect_warning(a <- AIC(res), "requires objective")
  expect_true(is.na(a))
})

test_that("AIC.STVariogramFit returns numeric for MLE fit", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(stf = stf, formula = z ~ 1,
      typestv = "sumMetric", candidate_model = c("Exc"),
      objective = "MLE")
  )
  a <- AIC(res)
  expect_true(is.numeric(a))
})

test_that("BIC.STVariogramFit returns numeric for MLE fit", {
  stf <- make_test_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf = stf, formula = z ~ 1,
      cutoff = 3e3, width = 500,
      typestv = "sumMetric",
      candidate_model = c("Exc"),
      objective = "MLE")
  )
  b <- BIC(res)
  expect_true(is.numeric(b))
  # BIC >= AIC (for n > e^2 ≈ 7.4; our n=30)
  expect_gte(b, AIC(res))
})
