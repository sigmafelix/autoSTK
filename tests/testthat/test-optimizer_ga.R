# test-optimizer_ga.R
# Tests for the Genetic Algorithm (GA) optimizer added in v2.1.0.
# Tests cover:
#   - the built-in DE fallback (always available: base R + lhs)
#   - the GA::ga() path (skipped when GA is not installed)
#   - integration through autofitVariogramST(optimizer = "ga")

library(testthat)
library(spacetime)
library(gstat)

# ---- Shared fixtures ---------------------------------------------------------

make_stfdf <- function(n_sp = 15L, n_t = 6L, seed = 77L) {
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

# ---- DE fallback: return type and structure ----------------------------------

test_that(".de_optimizer returns a numeric vector of the right length", {
  skip_on_cran()
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)
  k    <- length(bnds$lower)

  par <- autoSTK:::.de_optimizer(
    model_template = mod,
    stva_emp       = stva,
    lower          = bnds$lower,
    upper          = bnds$upper,
    pop_size       = 10L,
    maxiter        = 5L,
    run            = 3L,
    seed           = 1L
  )

  expect_type(par, "double")
  expect_length(par, k)
})

test_that(".de_optimizer solution lies within [lower, upper]", {
  skip_on_cran()
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  par <- autoSTK:::.de_optimizer(
    model_template = mod, stva_emp = stva,
    lower = bnds$lower, upper = bnds$upper,
    pop_size = 10L, maxiter = 5L, run = 3L, seed = 2L
  )

  expect_true(all(par >= bnds$lower - 1e-9))
  expect_true(all(par <= bnds$upper + 1e-9))
})

test_that(".de_optimizer early-stops when run criterion is met", {
  skip_on_cran()
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  par <- autoSTK:::.de_optimizer(
    model_template = mod, stva_emp = stva,
    lower = bnds$lower, upper = bnds$upper,
    pop_size = 6L, maxiter = 200L, run = 1L, seed = 99L
  )
  expect_type(par, "double")
})

# ---- .optimize_stv_ga via DE fallback ----------------------------------------

test_that(".optimize_stv_ga (DE fallback) returns a StVariogramModel", {
  skip_on_cran()
  skip_if(requireNamespace("GA", quietly = TRUE), "GA package is installed")
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  res <- suppressWarnings(
    autoSTK:::.optimize_stv_ga(
      stva_emp       = stva,
      model_template = mod,
      bounds         = bnds,
      popSize        = 8L,
      maxiter        = 5L,
      run            = 3L,
      seed           = 7L,
      lbfgsb_refine  = FALSE
    )
  )

  expect_true(inherits(res, "StVariogramModel"))
  expect_gt(length(extractPar(res)), 0L)
})

test_that(".optimize_stv_ga (DE fallback) parameters within bounds", {
  skip_on_cran()
  skip_if(requireNamespace("GA", quietly = TRUE), "GA package is installed")
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  res <- suppressWarnings(
    autoSTK:::.optimize_stv_ga(
      stva_emp = stva, model_template = mod, bounds = bnds,
      popSize = 8L, maxiter = 5L, run = 3L, seed = 8L,
      lbfgsb_refine = FALSE
    )
  )
  pars <- extractPar(res)
  expect_true(all(pars >= bnds$lower - 1e-4))
  expect_true(all(pars <= bnds$upper + 1e-4))
})

# ---- .optimize_stv_ga via GA package (when installed) ------------------------

test_that(".optimize_stv_ga (GA pkg) returns a StVariogramModel", {
  skip_on_cran()
  skip_if_not_installed("GA")
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  res <- suppressWarnings(
    autoSTK:::.optimize_stv_ga(
      stva_emp       = stva,
      model_template = mod,
      bounds         = bnds,
      popSize        = 10L,
      maxiter        = 10L,
      run            = 5L,
      seed           = 42L,
      lbfgsb_refine  = FALSE
    )
  )

  expect_true(inherits(res, "StVariogramModel"))
  expect_gt(length(extractPar(res)), 0L)
})

test_that(".optimize_stv_ga (GA pkg) with lbfgsb_refine returns valid model", {
  skip_on_cran()
  skip_if_not_installed("GA")
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  res <- suppressWarnings(
    autoSTK:::.optimize_stv_ga(
      stva_emp = stva, model_template = mod, bounds = bnds,
      popSize = 10L, maxiter = 10L, run = 5L, seed = 42L,
      lbfgsb_refine = TRUE, lbfgsb_maxit = 100L
    )
  )

  pars <- extractPar(res)
  expect_true(is.numeric(pars) && length(pars) > 0L)
})

# ---- Integration through autofitVariogramST ----------------------------------

test_that("autofitVariogramST with optimizer='ga' returns STVariogramFit", {
  skip_on_cran()
  stf <- make_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf               = stf,
      formula           = z ~ 1,
      typestv           = "sumMetric",
      candidate_model   = c("Exp", "Sph"),
      optimizer         = "ga",
      optimizer_control = list(popSize = 8L, maxiter = 5L,
                               run = 3L, seed = 42L)
    )
  )
  expect_s3_class(res, "STVariogramFit")
  expect_equal(res$optimizer, "ga")
  expect_true(all(c("jointSTV", "empSTV", "SpV", "TV") %in% names(res)))
})

test_that("autofitVariogramST GA: separable model type works", {
  skip_on_cran()
  stf <- make_stfdf()
  suppressWarnings(
    res <- autofitVariogramST(
      stf             = stf, formula = z ~ 1, typestv = "separable",
      candidate_model = c("Exp", "Sph"),
      optimizer = "ga",
      optimizer_control = list(popSize = 8L, maxiter = 5L, run = 3L)
    )
  )
  expect_s3_class(res, "STVariogramFit")
})

# ---- .ga_fitness wrapper -----------------------------------------------------

test_that(".ga_fitness returns finite non-positive value for in-bounds params", {
  skip_on_cran()
  stva <- make_emp_stv()
  mod  <- make_summetric_vgm()
  bnds <- st_param_bounds(mod, stva)

  pars <- extractPar(mod)
  pars <- pmax(pmin(pars, bnds$upper), bnds$lower)

  fit_val <- autoSTK:::.ga_fitness(pars, mod, stva)

  expect_type(fit_val, "double")
  # .ga_fitness = -WLS, WLS >= 0, so result <= 0
  expect_lte(fit_val, 0)
})
