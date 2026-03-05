# parallel_utils.R — thin parallelism abstraction.
#
# Uses furrr (future-based) when available, falls back to lapply.
# Callers pass .parallel = TRUE to opt in; the user controls the plan via
# future::plan() before calling package functions.

.par_map <- function(xs, fn, ..., .parallel = FALSE) {
  if (.parallel && requireNamespace("furrr", quietly = TRUE)) {
    furrr::future_map(xs, fn, ...,
                      .options = furrr::furrr_options(seed = TRUE))
  } else {
    lapply(xs, fn, ...)
  }
}
