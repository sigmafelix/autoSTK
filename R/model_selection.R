# model_selection.R — information criteria and automated model selection.

#' Compute AIC and BIC for a fitted STVariogramFit (MLE-based)
#'
#' @param fit_obj A \code{STVariogramFit} object with an \code{objective =
#'   "MLE"} fit (i.e., \code{fit_obj$loglik} must not be \code{NULL}).
#' @param n_obs Integer. Number of observations used for fitting.
#' @return A named list: \code{AIC}, \code{BIC}, \code{k} (number of free
#'   parameters), \code{loglik}.
#' @export
ic_stv <- function(fit_obj, n_obs) {
  if (is.null(fit_obj$loglik))
    stop("ic_stv() requires a fit with objective = 'MLE'.")
  k   <- length(extractPar(fit_obj$jointSTV))
  ll  <- fit_obj$loglik
  list(
    AIC    = -2 * ll + 2 * k,
    BIC    = -2 * ll + log(n_obs) * k,
    k      = k,
    loglik = ll
  )
}

#' Automated Spatio-Temporal Covariance Model Selection
#'
#' Fits multiple ST covariance model types and ranks them by a user-chosen
#' criterion.
#'
#' @param stf A spatio-temporal data object (STFDF, STSDF, STIDF, or sftime).
#' @param formula A formula (e.g., \code{PM10 ~ 1}).
#' @param candidates Character vector of ST model types to try. Defaults to
#'   all six supported types.
#' @param criterion One of \code{"AIC"}, \code{"BIC"} (require
#'   \code{objective = "MLE"}), \code{"CV_RMSE"} (not yet implemented in v2.0;
#'   use \code{"WLS"} as a fast proxy), or \code{"WLS"}.
#' @param optimizer Passed to \code{autofitVariogramST}.
#' @param cores Passed to \code{autofitVariogramST}.
#' @param ... Additional arguments forwarded to \code{autofitVariogramST}.
#' @return An \code{STModelSelection} object (list) with elements
#'   \code{best_model}, \code{comparison} (data.frame), and \code{all_fits}.
#' @export
selectModelST <- function(stf, formula,
                           candidates = c("separable", "productSum",
                                          "sumMetric", "simpleSumMetric",
                                          "metric"),
                           criterion  = c("WLS", "AIC", "BIC"),
                           optimizer  = "lbfgsb",
                           cores      = 1L,
                           ...) {
  criterion <- match.arg(criterion)

  objective <- switch(criterion,
    AIC = "MLE",
    BIC = "MLE",
    WLS = "WLS"
  )

  # Parallel-aware map over candidate model types
  fits <- .par_map(candidates, function(m) {
    tryCatch(
      autofitVariogramST(
        stf       = stf,
        formula   = formula,
        typestv   = m,
        optimizer = optimizer,
        objective = objective,
        cores     = cores,
        ...
      ),
      error = function(e) {
        message("Model '", m, "' failed: ", conditionMessage(e))
        NULL
      }
    )
  }, .parallel = FALSE)
  names(fits) <- candidates
  fits <- Filter(Negate(is.null), fits)

  if (length(fits) == 0L)
    stop("All candidate models failed to fit.")

  n_obs <- nrow(.to_stsdf(stf)@data)

  scores <- vapply(fits, function(f) {
    switch(criterion,
      AIC = tryCatch(ic_stv(f, n_obs)$AIC, error = function(e) Inf),
      BIC = tryCatch(ic_stv(f, n_obs)$BIC, error = function(e) Inf),
      WLS = {
        mse <- attr(f$jointSTV, "MSErr")
        if (is.null(mse) || !is.finite(mse)) Inf else mse
      }
    )
  }, numeric(1L))

  comparison <- data.frame(
    model     = names(scores),
    score     = unname(scores),
    criterion = criterion,
    stringsAsFactors = FALSE
  )
  comparison <- comparison[order(comparison$score), ]

  structure(
    list(
      best_model = fits[[comparison$model[1L]]],
      comparison = comparison,
      all_fits   = fits,
      criterion  = criterion
    ),
    class = "STModelSelection"
  )
}
