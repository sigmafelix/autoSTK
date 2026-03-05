# test-coerce_st.R

library(testthat)
library(spacetime)
library(sftime)

make_stfdf <- function() {
  sp   <- sp::SpatialPoints(cbind(x = c(0, 1, 2), y = c(0, 1, 2)))
  time <- as.POSIXct("2020-01-01") + 0:2 * 3600
  STFDF(sp, time, data.frame(z = rnorm(9L)))
}

test_that(".to_stsdf works on STFDF input", {
  stfdf  <- make_stfdf()
  result <- autoSTK:::.to_stsdf(stfdf)
  expect_s4_class(result, "STSDF")
})

test_that(".to_stsdf is identity on STSDF input", {
  stsdf  <- as(make_stfdf(), "STSDF")
  result <- autoSTK:::.to_stsdf(stsdf)
  expect_s4_class(result, "STSDF")
  expect_equal(nrow(result@data), nrow(stsdf@data))
})

test_that(".to_stsdf works on STIDF input", {
  stidf  <- as(make_stfdf(), "STIDF")
  result <- autoSTK:::.to_stsdf(stidf)
  expect_s4_class(result, "STSDF")
})

test_that(".to_stfdf works on STSDF input", {
  stsdf  <- as(make_stfdf(), "STSDF")
  result <- autoSTK:::.to_stfdf(stsdf)
  expect_s4_class(result, "STFDF")
})

test_that(".to_stfdf is identity on STFDF input", {
  stfdf  <- make_stfdf()
  result <- autoSTK:::.to_stfdf(stfdf)
  expect_s4_class(result, "STFDF")
})

test_that(".to_stsdf errors on unsupported class", {
  expect_error(autoSTK:::.to_stsdf(list()), "Cannot coerce")
})

test_that("autofitVariogram accepts sf and sftime objects", {
  skip_on_cran()
  set.seed(42)

  d <- data.frame(
    x = runif(30L, 0, 5000),
    y = runif(30L, 0, 5000),
    z = rnorm(30L),
    time = as.POSIXct("2020-01-01", tz = "UTC") + seq_len(30L)
  )

  sf_obj <- sf::st_as_sf(d, coords = c("x", "y"), crs = 3857)
  fit_sf <- autoSTK:::autofitVariogram(z ~ 1, input_data = sf_obj, model = "Sph")
  expect_s3_class(fit_sf, "autofitVariogram")

  sftime_obj <- sftime::st_as_sftime(sf_obj, time_column_name = "time")
  fit_sftime <- autoSTK:::autofitVariogram(z ~ 1, input_data = sftime_obj, model = "Sph")
  expect_s3_class(fit_sftime, "autofitVariogram")
})
