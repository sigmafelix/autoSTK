# test-optimizer_sa.R
# Tests for the Simulated Annealing (SA) optimizer added in v2.1.0.

library(testthat)
library(spacetime)
library(gstat)

# ---- Shared fixtures ---------------------------------------------------------

# 15 stations over a 10 km x 10 km grid — same size used in existing tests.
make_stfdf <- function(n_sp = 15L, n_t = 6L, seed = 11L) {
  set.seed(seed)
  sp <- sp::SpatialPoints(cbind(
    x = runif(n_sp, 0, 1e4),
    y = runif(n_sp, 0, 1e4)
  ))
  sp::proj4string(sp) <- sp::CRS("+proj=utm +zone=32 +datum=WGS84")
  time <- as.POSIXct("2020-01-01", tz = "UTC") +
            seq(0, (n_t - 1L) * 3600, by = 3600)
  STFDF(sp, time, data.frame(z = rnorm(n_sp * n_t)))
}

make_emp_stv <- function(stf = NULL) {
  if (is.null(stf)) stf <- make_stfdf()
  autoSTK::setSTI(stf, z ~ 1, tlags = 0:3,
                  cutoff = 6000, width = 1000, wireframe = FALSE, cores = 1L)
}

make_summetric_vgm <- function() {
  gstat::vgmST(
    "sumMetric",
    space = gstat::vgm(psill = 0.8, model = "Exp", range = 2500, nugget = 0.1),
    time  = gstat::vgm(psill = 0.8, model = "Exp", range = 2.5,  nugget = 0.1),
    joint = gstat::vgm(psill = 0.5, model = "Exp", range = 1500, nugget = 0.05),
    stAni = 500
  )
}

# ---- Return type and structure -----------------------------------------------

test_that(".optimize_stv_sa returns a StVariogramModel", {
  skip_on_cran()
  stva  <- make_emp_stv()
  mod   <- make_summetric_vgm()
  bnds  <- st_param_bounds(mod, stva)

  res <- suppressWarnings(
    autoSTK:::.optimize_stv_sa(
      stva_emp       = stva,
      model_template = mod,
      bounds         = bnds,
      maxit          = 50L,
      temp           = 5,
      tmax           = 5L,
      lbfgsb_refine  = FALSE
    )
  )

  expect_true(inherits(res, "StVariogramModel"))
})

test_that(".optimize_stv_sa result has extractable parameters", {
  skip_on_cran()
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  res  <- suppressWarnings(
    autoSTK:::.optimize_stv_sa(stva_emp = stva, model_template = mod,
                                bounds = bnds, maxit = 50L,
                                lbfgsb_refine = FALSE)
  )
  pars <- extractPar(res)
  expect_true(is.numeric(pars))
  expect_gt(length(pars), 0L)
  expect_true(all(is.finite(pars)))
})

# ---- Parameter bounds enforcement --------------------------------------------

test_that("SA solution parameters stay within supplied bounds", {
  skip_on_cran()
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  res  <- suppressWarnings(
    autoSTK:::.optimize_stv_sa(stva_emp = stva, model_template = mod,
                                bounds = bnds, maxit = 80L,
                                lbfgsb_refine = FALSE)
  )
  pars <- extractPar(res)
  expect_true(all(pars >= bnds$lower - 1e-4))
  expect_true(all(pars <= bnds$upper + 1e-4))
})

# ---- lbfgsb_refine flag ------------------------------------------------------

test_that("SA with lbfgsb_refine = TRUE also returns a valid model", {
  skip_on_cran()
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  res <- suppressWarnings(
    autoSTK:::.optimize_stv_sa(stva_emp = stva, model_template = mod,
                                bounds = bnds, maxit = 50L,
                                lbfgsb_refine = TRUE, lbfgsb_maxit = 100L)
  )
  expect_true(inherits(res, "StVariogramModel"))
  pars <- extractPar(res)
  expect_true(is.numeric(pars) && length(pars) > 0L)
})

# ---- Integration through autofitVariogramST ----------------------------------

test_that("autofitVariogramST with optimizer='sa' returns STVariogramFit", {
  skip_on_cran()
  stf <- make_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf               = stf,
      formula           = z ~ 1,
      typestv           = "sumMetric",
      candidate_model   = c("Exp", "Sph"),
      optimizer         = "sa",
      optimizer_control = list(maxit = 80L, temp = 5, tmax = 5L)
    )
  )
  expect_s3_class(res, "STVariogramFit")
  expect_equal(res$optimizer, "sa")
  expect_true(all(c("jointSTV", "empSTV", "SpV", "TV") %in% names(res)))
})

test_that("autofitVariogramST SA: extractPar gives named numeric vector", {
  skip_on_cran()
  stf <- make_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf             = stf, formula = z ~ 1, typestv = "separable",
      candidate_model = c("Exp", "Sph"),
      optimizer = "sa",
      optimizer_control = list(maxit = 60L)
    )
  )
  pars <- extractPar(res$jointSTV)
  expect_true(is.numeric(pars))
  expect_false(is.null(names(pars)))
})

# ---- Sigmoid transform round-trip --------------------------------------------

test_that("sigmoid <-> box transforms are inverse of each other", {
  lower <- c(0, 0, 100,   0,   0,   0, 10)
  upper <- c(2, 2, 5000, 10, 10, 10, 200)
  set.seed(7L)
  x <- lower + runif(length(lower)) * (upper - lower)

  u         <- autoSTK:::.box_to_sigmoid(x, lower, upper)
  x_recover <- autoSTK:::.sigmoid_to_box(u, lower, upper)

  expect_equal(x_recover, x, tolerance = 1e-9)
})

test_that("sigmoid transform outputs are finite for in-bounds input", {
  lower <- rep(0,    5L)
  upper <- rep(1000, 5L)
  x     <- c(100, 200, 500, 800, 900)
  u     <- autoSTK:::.box_to_sigmoid(x, lower, upper)
  expect_true(all(is.finite(u)))
})

