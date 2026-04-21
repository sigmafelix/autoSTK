# test-cubble_to_sftime.R
# Tests for cubble_to_sftime() and the cubble dispatch in
# .to_stsdf() / .to_stfdf().
# All tests are guarded by skip_if_not_installed("cubble") because cubble is
# in Suggests, not Imports.

library(testthat)
library(sf)
library(sftime)

# ---- Fixture -----------------------------------------------------------------

# Build a minimal nested cubble: 4 stations x 3 time steps.
# Returns a nested-face cubble_df.
make_cubble_fixture <- function(n_sp = 4L, n_t = 3L, seed = 42L) {
  skip_if_not_installed("cubble")

  set.seed(seed)

  sp_df <- sf::st_as_sf(
    data.frame(
      id = paste0("s", seq_len(n_sp)),
      x  = runif(n_sp, 0, 1e4),
      y  = runif(n_sp, 0, 1e4)
    ),
    coords = c("x", "y"),
    crs    = 3857
  )

  times <- as.POSIXct("2020-01-01", tz = "UTC") + seq_len(n_t) * 3600

  t_df <- data.frame(
    id    = rep(paste0("s", seq_len(n_sp)), each = n_t),
    time  = rep(times, times = n_sp),
    value = rnorm(n_sp * n_t)
  )

  cubble::make_cubble( # nolint: object_usage_linter
    spatial  = sp_df,
    temporal = t_df,
    key      = id,   # NSE — id is a column name, not a variable
    index    = time  # NSE — time is a column name, not a variable
  )
}

# ---- Output class / dimensions -----------------------------------------------

test_that("cubble_to_sftime returns an sftime object", {
  skip_if_not_installed("cubble")
  cb  <- make_cubble_fixture()
  out <- cubble_to_sftime(cb)
  expect_s3_class(out, "sftime")
})

test_that("output has n_sp * n_t rows", {
  skip_if_not_installed("cubble")
  n_sp <- 4L
  n_t  <- 3L
  cb  <- make_cubble_fixture(n_sp = n_sp, n_t = n_t)
  out <- cubble_to_sftime(cb)
  expect_equal(nrow(out), n_sp * n_t)
})

test_that("output retains the data column 'value'", {
  skip_if_not_installed("cubble")
  cb  <- make_cubble_fixture()
  out <- cubble_to_sftime(cb)
  expect_true("value" %in% names(out))
})

test_that("output has a point geometry column", {
  skip_if_not_installed("cubble")
  cb  <- make_cubble_fixture()
  out <- cubble_to_sftime(cb)
  geom <- sf::st_geometry(out)
  expect_s3_class(geom, "sfc_POINT")
})

test_that("CRS is preserved from the spatial face", {
  skip_if_not_installed("cubble")
  cb  <- make_cubble_fixture()
  out <- cubble_to_sftime(cb)
  expect_equal(sf::st_crs(out)$epsg, 3857L)
})

# ---- Temporal face as input --------------------------------------------------

test_that("cubble_to_sftime works when given the temporal face directly", {
  skip_if_not_installed("cubble")
  cb      <- make_cubble_fixture()
  cb_temp <- cubble::face_temporal(cb)
  # temporal face: is_cubble_spatial = FALSE, is_cubble_temporal = TRUE
  expect_true(cubble::is_cubble_temporal(cb_temp))
  out <- cubble_to_sftime(cb_temp, key_col = "id", time_col = "time")
  expect_s3_class(out, "sftime")
  expect_equal(nrow(out), 4L * 3L)
})

# ---- Auto-detection ----------------------------------------------------------

test_that("time_col auto-detection finds the POSIXct column", {
  skip_if_not_installed("cubble")
  cb <- make_cubble_fixture()
  # Pass NULL explicitly to force auto-detection
  expect_message(
    out <- cubble_to_sftime(cb, time_col = NULL),
    regexp = "time"   # message mentions detected column name
  )
  expect_s3_class(out, "sftime")
})

test_that("error when no POSIXct/Date column exists and time_col is NULL", {
  skip_if_not_installed("cubble")

  # Build cubble where temporal face has only numeric columns (no datetime)
  set.seed(1L)
  sp_df <- sf::st_as_sf(
    data.frame(id = c("a", "b"), x = c(0, 1), y = c(0, 1)),
    coords = c("x", "y"), crs = 3857
  )
  t_df <- data.frame(
    id    = rep(c("a", "b"), each = 2L),
    step  = rep(1:2, times = 2L),   # integer, not POSIXct
    value = rnorm(4L)
  )
  cb <- cubble::make_cubble( # nolint: object_usage_linter
    spatial  = sp_df,
    temporal = t_df,
    key      = id,   # NSE
    index    = step  # NSE
  )

  expect_error(
    cubble_to_sftime(cb, time_col = NULL),
    regexp = "time column"
  )
})

test_that("explicit time_col overrides auto-detection", {
  skip_if_not_installed("cubble")
  cb  <- make_cubble_fixture()
  out <- cubble_to_sftime(cb, time_col = "time")
  expect_s3_class(out, "sftime")
})

# ---- Coercion helpers (.to_stsdf / .to_stfdf) --------------------------------

test_that(".to_stsdf dispatches cubble_df to STSDF", {
  skip_if_not_installed("cubble")
  cb  <- make_cubble_fixture()
  out <- autoSTK:::.to_stsdf(cb)
  expect_s4_class(out, "STSDF")
})

test_that(".to_stfdf dispatches cubble_df to STFDF", {
  skip_if_not_installed("cubble")
  cb  <- make_cubble_fixture()
  out <- autoSTK:::.to_stfdf(cb)
  expect_s4_class(out, "STFDF")
})

test_that(".to_stsdf from cubble has the correct observation count", {
  skip_if_not_installed("cubble")
  n_sp <- 4L
  n_t  <- 3L
  cb  <- make_cubble_fixture(n_sp = n_sp, n_t = n_t)
  out <- autoSTK:::.to_stsdf(cb)
  expect_equal(nrow(out@data), n_sp * n_t)
})

# ---- Error message mentions cubble when class unsupported --------------------

test_that(".to_stsdf error on unknown class mentions cubble_df", {
  expect_error(
    autoSTK:::.to_stsdf(list()),
    regexp = "cubble_df"
  )
})

test_that(".to_stfdf error on unknown class mentions cubble_df", {
  expect_error(
    autoSTK:::.to_stfdf(42L),
    regexp = "cubble_df"
  )
})
