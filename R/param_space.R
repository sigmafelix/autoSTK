# param_space.R — semantically-aware parameter bounds and data-driven
# initial parameter estimation for vgmST models.
#
# The key improvement over the existing 0.25x / 1.5x heuristic:
# parameter names returned by gstat::extractPar() are used to classify each
# parameter as a sill/psill, nugget, range, anisotropy ratio, or k-factor,
# and each class gets bounds derived from the empirical variogram statistics.
#
# estimate_initial_params() replaces the ad-hoc arithmetic guesses in
# autofitVariogramST with a robust, data-driven approach that locates the
# practical range (lag where variogram first reaches 95 % of its sill) for
# both the spatial and temporal marginals.


#' Estimate Initial ST Variogram Parameters from Empirical Data
#'
#' Derives nugget, partial sill, spatial range, and temporal range directly
#' from the empirical spatio-temporal variogram instead of relying on fixed
#' arithmetic heuristics.  The approach is:
#'
#' \enumerate{
#'   \item \strong{Nugget} — median \eqn{\gamma} over the first decile of
#'     spatial lags (extrapolation toward zero-lag).
#'   \item \strong{Sill} — 95th percentile of \eqn{\gamma} values over all
#'     spatial lags at the minimum time lag (spatial marginal plateau).
#'   \item \strong{Spatial practical range} — smallest spatial lag at which
#'     the spatial marginal exceeds 95\% of the estimated sill; falls back to
#'     50\% of the maximum spatial lag when the variogram never reaches the
#'     sill within the sampling window.
#'   \item \strong{Temporal practical range} — same rule applied to the
#'     temporal marginal at the minimum spatial lag.
#' }
#'
#' @param stva_emp An empirical \code{StVariogram} (data.frame with columns
#'   \code{spacelag}, \code{timelag}, \code{gamma}, \code{np}).
#' @param nugget_quantile Numeric in (0, 1). Quantile of near-zero-lag gamma
#'   used as the nugget estimate.  Default 0.25 (lower quartile of the first
#'   decile of spatial bins).
#' @param sill_quantile Numeric in (0, 1). Quantile of the spatial-marginal
#'   gamma values used as the plateau estimate.  Default 0.95.
#' @param practical_range_target Numeric in (0, 1). Fraction of the sill at
#'   which the practical range is read off.  Default 0.95.
#'
#' @return A named list:
#'   \describe{
#'     \item{nugget}{Estimated nugget (>= 0).}
#'     \item{psill}{Estimated partial sill (sill - nugget, >= 0).}
#'     \item{sill}{Estimated total sill.}
#'     \item{sp_range}{Estimated spatial practical range.}
#'     \item{ts_range}{Estimated temporal practical range.}
#'     \item{max_gamma}{Maximum observed gamma.}
#'     \item{max_splag}{Maximum spatial lag in the empirical variogram.}
#'     \item{max_tlag}{Maximum temporal lag (numeric).}
#'   }
#' @export
estimate_initial_params <- function(stva_emp,
                                    nugget_quantile         = 0.25,
                                    sill_quantile           = 0.95,
                                    practical_range_target  = 0.95) {

  # ---- Guard against missing columns ----------------------------------------
  required_cols <- c("spacelag", "timelag", "gamma")
  missing_cols  <- setdiff(required_cols, names(stva_emp))
  if (length(missing_cols) > 0L)
    stop("stva_emp is missing columns: ", paste(missing_cols, collapse = ", "))

  gamma    <- stva_emp$gamma
  spl      <- stva_emp$spacelag
  tl       <- as.numeric(stva_emp$timelag)

  # Remove NA/infinite rows
  valid    <- is.finite(gamma) & is.finite(spl) & is.finite(tl)
  gamma    <- gamma[valid]
  spl      <- spl[valid]
  tl       <- tl[valid]

  if (length(gamma) == 0L)
    stop("No valid (finite) gamma values found in stva_emp.")

  # ---- Spatial marginal (minimum time lag) ----------------------------------
  min_tl      <- min(tl)
  sp_marg_idx <- tl == min_tl
  sp_gamma    <- gamma[sp_marg_idx]
  sp_lags     <- spl[sp_marg_idx]

  # Sort by spatial lag for range look-up
  sp_ord   <- order(sp_lags)
  sp_gamma <- sp_gamma[sp_ord]
  sp_lags  <- sp_lags[sp_ord]

  # ---- Temporal marginal (minimum spatial lag) ------------------------------
  min_spl     <- min(spl)
  ts_marg_idx <- spl == min_spl
  ts_gamma    <- gamma[ts_marg_idx]
  ts_lags     <- tl[ts_marg_idx]

  ts_ord   <- order(ts_lags)
  ts_gamma <- ts_gamma[ts_ord]
  ts_lags  <- ts_lags[ts_ord]

  # ---- Nugget ---------------------------------------------------------------
  # Use the lower quantile of gamma within the first decile of spatial bins
  # at the minimum time lag.  This is more robust than 0.5 * min(gamma).
  n_sp         <- length(sp_lags)
  decile_idx   <- seq_len(max(1L, ceiling(n_sp * 0.10)))
  near_zero    <- sp_gamma[decile_idx]
  nugget_est   <- max(0,
                      stats::quantile(near_zero, nugget_quantile, na.rm = TRUE))

  # ---- Sill (plateau of spatial marginal) -----------------------------------
  sill_est     <- as.numeric(
    stats::quantile(sp_gamma, sill_quantile, na.rm = TRUE)
  )
  # Make sure sill >= nugget
  sill_est     <- max(sill_est, nugget_est + 1e-6)
  psill_est    <- sill_est - nugget_est

  # ---- Practical ranges -----------------------------------------------------
  gamma_target <- nugget_est + practical_range_target * psill_est

  # Spatial practical range
  above_sp  <- which(sp_gamma >= gamma_target)
  if (length(above_sp) > 0L) {
    sp_range_est <- sp_lags[above_sp[1L]]
  } else {
    sp_range_est <- max(sp_lags) * 0.5
  }
  sp_range_est <- max(sp_range_est, 1e-6)

  # Temporal practical range — same logic on temporal marginal
  if (length(ts_gamma) > 0L) {
    ts_sill        <- as.numeric(
      stats::quantile(ts_gamma, sill_quantile, na.rm = TRUE)
    )
    ts_sill        <- max(ts_sill, nugget_est + 1e-6)
    ts_psill       <- ts_sill - nugget_est
    ts_target      <- nugget_est + practical_range_target * ts_psill
    above_ts       <- which(ts_gamma >= ts_target)
    if (length(above_ts) > 0L) {
      ts_range_est <- ts_lags[above_ts[1L]]
    } else {
      ts_range_est <- max(ts_lags) * 0.5
    }
    ts_range_est   <- max(ts_range_est, 1e-6)
  } else {
    ts_range_est   <- max(tl) * 0.5
  }

  list(
    nugget    = nugget_est,
    psill     = psill_est,
    sill      = sill_est,
    sp_range  = as.numeric(sp_range_est),
    ts_range  = as.numeric(ts_range_est),
    max_gamma = max(gamma, na.rm = TRUE),
    max_splag = max(spl,   na.rm = TRUE),
    max_tlag  = max(tl,    na.rm = TRUE)
  )
}

#' Compute parameter bounds for vgmST model optimisation
#'
#' @param model_template A \code{vgmST} object (initial guess).
#' @param stva_emp Empirical ST variogram (\code{StVariogram}/data.frame) from
#'   \code{setSTI} or \code{gstat::variogramST}.
#' @param sill_scale Numeric (2). Upper bound multiplier applied to
#'   \code{max(stva_emp$gamma)} for sill/psill/nugget parameters.
#' @param range_scale Numeric (3). Upper bound multiplier applied to the
#'   maximum spatial lag for range parameters.
#' @param ani_scale Numeric (20). Controls the breadth of the stAni search
#'   interval relative to the data extent.
#' @return A named list with elements \code{lower} and \code{upper} (both
#'   named numeric vectors matching \code{gstat::extractPar(model_template)}).
#' @importFrom gstat extractPar
#' @export
st_param_bounds <- function(model_template, stva_emp,
                             sill_scale  = 2.0,
                             range_scale = 3.0,
                             ani_scale   = 20.0) {
  init_par  <- gstat::extractPar(model_template)
  par_names <- names(init_par)

  maxgamma <- max(stva_emp$gamma, na.rm = TRUE)
  maxspl   <- max(stva_emp$spacelag, na.rm = TRUE)
  maxtl    <- as.numeric(max(stva_emp$timelag, na.rm = TRUE))

  lower <- setNames(numeric(length(init_par)), par_names)
  upper <- setNames(numeric(length(init_par)), par_names)

  for (i in seq_along(init_par)) {
    nm <- if (!is.null(par_names)) par_names[i] else ""
    v  <- init_par[i]

    if (grepl("sill|psill|nugget", nm, ignore.case = TRUE)) {
      lower[i] <- 0
      upper[i] <- maxgamma * sill_scale

    } else if (grepl("^range\\.", nm, ignore.case = TRUE) ||
               identical(nm, "range")) {
      # Spatial range
      lower[i] <- 1e-6
      upper[i] <- maxspl * range_scale

    } else if (grepl("^range\\.t", nm, ignore.case = TRUE)) {
      # Temporal range
      lower[i] <- 1e-6
      upper[i] <- maxtl * range_scale

    } else if (grepl("^(stAni|anis)", nm, ignore.case = TRUE) ||
               identical(nm, "anis")) {
      # Anisotropy ratio: space units per time unit
      # Lower: at least a small positive ratio
      # Upper: allow a wide range up to ani_scale * max_space / max_time
      lower[i] <- max(maxspl / (maxtl * ani_scale), 1e-6)
      upper[i] <- maxspl * ani_scale / max(maxtl, 1e-6)

    } else {
      # k (productSum) or unknown: relative bounds from initial value
      if (v > 0) {
        lower[i] <- v * 0.01
        upper[i] <- v * 5.0
      } else {
        lower[i] <- 1e-6
        upper[i] <- maxgamma * sill_scale
      }
    }
  }

  # Safety: ensure lower >= 0 and upper > lower
  lower <- pmax(lower, 0)
  upper <- pmax(upper, lower + 1e-6)

  list(lower = lower, upper = upper)
}
