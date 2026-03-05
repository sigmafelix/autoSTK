# generics.R — S3 print, summary, AIC, BIC for STVariogramFit
#              and STModelSelection.

#' @export
print.STVariogramFit <- function(x, ...) {
  cat("Spatio-Temporal Variogram Fit\n")
  cat("  ST model type :", attr(x$jointSTV, "stModel"), "\n")
  cat("  Optimizer     :", x$optimizer, "\n")
  cat("  Objective     :", x$objective, "\n")
  mse <- attr(x$jointSTV, "MSErr")
  if (!is.null(mse)) cat("  MSErr         :", signif(mse, 5L), "\n")
  if (!is.null(x$loglik))
    cat("  log-likelihood:", signif(x$loglik, 6L), "\n")
  invisible(x)
}

#' @export
summary.STVariogramFit <- function(object, ...) {
  cat("=== STVariogramFit Summary ===\n\n")
  cat("ST model type :", attr(object$jointSTV, "stModel"), "\n")
  cat("Optimizer     :", object$optimizer, "\n")
  cat("Objective     :", object$objective, "\n\n")

  mse <- attr(object$jointSTV, "MSErr")
  if (!is.null(mse)) cat("MSErr         :", signif(mse, 5L), "\n")
  if (!is.null(object$loglik)) {
    cat("log-likelihood:", signif(object$loglik, 6L), "\n")
    k <- length(extractPar(object$jointSTV))
    n <- object$n_obs
    if (!is.null(n)) {
      aic <- -2 * object$loglik + 2 * k
      bic <- -2 * object$loglik + log(n) * k
      cat("AIC           :", signif(aic, 6L), "\n")
      cat("BIC           :", signif(bic, 6L), "\n")
    }
  }

  cat("\nJoint ST variogram parameters:\n")
  print(extractPar(object$jointSTV))
  cat("\nSpatial marginal:\n")
  print(object$SpV$var_model)
  cat("\nTemporal marginal:\n")
  print(object$TV$var_model)
  invisible(object)
}

#' @export
AIC.STVariogramFit <- function(object, ..., k = 2) {
  if (is.null(object$loglik)) {
    warning("AIC requires objective = 'MLE'. Returning NA.")
    return(NA_real_)
  }
  n_par <- length(extractPar(object$jointSTV))
  -2 * object$loglik + k * n_par
}

#' @export
BIC.STVariogramFit <- function(object, ...) {
  if (is.null(object$loglik)) {
    warning("BIC requires objective = 'MLE'. Returning NA.")
    return(NA_real_)
  }
  if (is.null(object$n_obs)) {
    warning("BIC requires n_obs to be stored in the fit. Returning NA.")
    return(NA_real_)
  }
  n_par <- length(extractPar(object$jointSTV))
  -2 * object$loglik + log(object$n_obs) * n_par
}

#' @export
print.STModelSelection <- function(x, ...) {
  cat("Spatio-Temporal Covariance Model Selection\n")
  cat("Criterion:", x$criterion, "\n\n")
  print(x$comparison)
  cat("\nBest model:", as.character(x$comparison$model[1L]), "\n")
  invisible(x)
}
