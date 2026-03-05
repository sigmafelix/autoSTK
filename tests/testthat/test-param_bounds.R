# test-param_bounds.R

library(testthat)
library(gstat)
library(spacetime)

make_emp_stv <- function() {
  sp   <- sp::SpatialPoints(cbind(x = runif(15L, 0, 1e4),
                                  y = runif(15L, 0, 1e4)))
  sp::proj4string(sp) <- sp::CRS("+proj=utm +zone=32 +datum=WGS84")
  time <- as.POSIXct("2020-01-01") + 0:5 * 3600
  stf  <- STFDF(sp, time, data.frame(z = rnorm(90L)))
  autoSTK::setSTI(
    stf = stf,
    formula = z ~ 1,
    tlags = 0:3,
    cutoff = 6000,
    width = 1000,
    wireframe = FALSE,
    cores = 1
  )
}

make_summetric_model <- function() {
  sp_vgm <- vgm(psill = 1.0, model = "Exp", range = 2500, nugget = 0.1)
  tm_vgm <- vgm(psill = 1.0, model = "Exp", range = 2.0, nugget = 0.1)
  vgmST("sumMetric",
        space = sp_vgm,
        time  = tm_vgm,
        joint = vgm(psill = 1.0, model = "Exp", range = 2000, nugget = 0.1),
        stAni = 500)
}

test_that("st_param_bounds returns lower < upper for all parameters", {
  skip_on_cran()
  set.seed(1L)
  stva <- make_emp_stv()
  mod  <- make_summetric_model()
  b    <- st_param_bounds(mod, stva)
  expect_true(all(b$lower >= 0))
  expect_true(all(b$upper > b$lower))
  expect_equal(length(b$lower), length(b$upper))
  expect_equal(length(b$lower), length(extractPar(mod)))
})

test_that("st_param_bounds names match extractPar names", {
  skip_on_cran()
  set.seed(1L)
  stva <- make_emp_stv()
  mod  <- make_summetric_model()
  b    <- st_param_bounds(mod, stva)
  ep   <- extractPar(mod)
  expect_equal(names(b$lower), names(ep))
  expect_equal(names(b$upper), names(ep))
})
