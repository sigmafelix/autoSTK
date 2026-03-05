# diagnostics.R — separability test and anisotropy diagnostics.

#' Likelihood Ratio Test for Spatiotemporal Separability
#'
#' Tests whether a separable covariance structure is adequate versus the more
#' general sumMetric model using a likelihood ratio test (LRT).  Both models
#' are fitted with \code{objective = "MLE"}.
#'
#' @param stf Spatio-temporal data (STFDF, STSDF, STIDF, or sftime).
#' @param formula Formula, e.g. \code{PM10 ~ 1}.
#' @param ... Additional arguments forwarded to \code{autofitVariogramST}
#'   (e.g. \code{tlags}, \code{cutoff}, \code{width}, \code{cores}).
#' @return A list with elements:
#'   \describe{
#'     \item{statistic}{LRT chi-squared statistic.}
#'     \item{df}{Degrees of freedom (difference in number of parameters).}
#'     \item{p_value}{p-value from chi-squared distribution.}
#'     \item{separable_fit}{Fitted \code{STVariogramFit} for the separable model.}
#'     \item{summetric_fit}{Fitted \code{STVariogramFit} for the sumMetric model.}
#'   }
#' @export
test_separability <- function(stf, formula, ...) {
  message("Fitting separable model...")
  fit_sep <- tryCatch(
    autofitVariogramST(stf, formula, typestv = "separable",
                       objective = "MLE", ...),
    error = function(e) stop("Separable fit failed: ", conditionMessage(e))
  )

  message("Fitting sumMetric model...")
  fit_sum <- tryCatch(
    autofitVariogramST(stf, formula, typestv = "sumMetric",
                       objective = "MLE", ...),
    error = function(e) stop("sumMetric fit failed: ", conditionMessage(e))
  )

  if (is.null(fit_sep$loglik) || is.null(fit_sum$loglik))
    stop("Both models must be fitted with objective = 'MLE' for the LRT.")

  lrt_stat <- 2 * (fit_sum$loglik - fit_sep$loglik)
  df_diff  <- length(extractPar(fit_sum$jointSTV)) -
              length(extractPar(fit_sep$jointSTV))

  if (df_diff <= 0) {
    warning("sumMetric has fewer or equal parameters than separable; LRT ",
            "may not be valid.")
    df_diff <- abs(df_diff) + 1L
  }

  p_value <- stats::pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)

  list(
    statistic      = lrt_stat,
    df             = df_diff,
    p_value        = p_value,
    separable_fit  = fit_sep,
    summetric_fit  = fit_sum
  )
}

#' Sensitivity Plot of the Spatio-Temporal Anisotropy Ratio
#'
#' Evaluates \code{gstat::estiStAni} at a sequence of spatial intervals and
#' plots how the estimated anisotropy ratio varies.  Useful for diagnosing
#' whether the single-interval estimate used in \code{autofitVariogramST} is
#' stable.
#'
#' @param stva_emp Empirical ST variogram from \code{setSTI}.
#' @param spatial_vgm Fitted spatial marginal variogram model.
#' @param temporal_vgm Fitted temporal marginal variogram model.
#' @param n_intervals Integer. Number of evaluation intervals.
#' @param plot Logical. Whether to produce a base-R plot.
#' @return Invisibly, a data.frame with columns \code{interval} and
#'   \code{stAni}.
#' @export
plot_aniso_sensitivity <- function(stva_emp, spatial_vgm, temporal_vgm,
                                   n_intervals = 20L, plot = TRUE) {
  maxspl    <- max(stva_emp$spacelag, na.rm = TRUE)
  intervals <- seq(0.05, 0.95, length.out = n_intervals) * median(stva_emp$spacelag)

  ratios <- vapply(intervals, function(iv) {
    tryCatch(
      estiStAni(stva_emp,
                interval     = c(iv, 2 * iv),
                spatialVgm   = spatial_vgm,
                temporalVgm  = temporal_vgm),
      error   = function(e) NA_real_,
      warning = function(w) suppressWarnings(
        tryCatch(
          estiStAni(stva_emp, interval = c(iv, 2 * iv),
                    spatialVgm = spatial_vgm, temporalVgm = temporal_vgm),
          error = function(e2) NA_real_
        )
      )
    )
  }, numeric(1L))

  out <- data.frame(interval = intervals, stAni = ratios)

  if (plot) {
    plot(out$interval, out$stAni, type = "b",
         xlab = "Estimation interval (space units)",
         ylab = "Estimated stAni",
         main = "Sensitivity of space-time anisotropy ratio")
  }

  invisible(out)
}
