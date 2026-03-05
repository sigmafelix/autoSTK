# objectives.R — WLS and Gaussian MLE objective functions for vgmST fitting.

# ---- WLS ------------------------------------------------------------------

# Evaluate the WLS (weighted least squares) score for a vgmST model against
# the empirical ST variogram.  Weights are n_h / gamma_model(h,u)^2 (Cressie
# 1985).  Returns a scalar (MSErr from fit.StVariogram) or Inf on failure.
.wls_stv <- function(model, stva_emp) {
  tryCatch({
    fitted <- fit.StVariogram(stva_emp, model, method = "L-BFGS-B")
    mse <- attr(fitted, "MSErr")
    if (is.null(mse) || !is.finite(mse)) Inf else mse
  }, error = function(e) Inf)
}

# ---- MLE ------------------------------------------------------------------

# Compute the Gaussian log-likelihood for a fitted vgmST model given
# observation data.  Only suitable for small n (n_max default 500) because it
# builds an n x n covariance matrix.
#
# Arguments:
#   model        — fitted vgmST object
#   data         — STSDF (or coercible); must already be coerced before call
#   response_col — character name of the response column in data@data
#   n_max        — hard cap on n; returns -Inf if exceeded
#
# Returns a scalar log-likelihood.
.loglik_stv <- function(model, stsdf, response_col, n_max = 500L) {
  y <- stsdf@data[[response_col]]
  n <- length(y)

  if (n > n_max) {
    warning("MLE: n (", n, ") exceeds n_max (", n_max,
            "). Returning -Inf. Consider objective = 'WLS' for large datasets.")
    return(-Inf)
  }

  # Build pairwise lag matrix for all observations
  lag_df <- .obs_lag_pairs(stsdf)

  # Evaluate variogram at all pairs
  gamma_vals <- tryCatch(
    variogramSurface(model, lag_df)$gamma,
    error = function(e) return(rep(Inf, nrow(lag_df)))
  )

  if (any(!is.finite(gamma_vals))) return(-Inf)

  Gamma <- matrix(gamma_vals, n, n)

  # Total sill (process variance = C(0,0)) by evaluating at very large lag
  C00 <- tryCatch(
    variogramSurface(model,
                     data.frame(spacelag = 1e10, timelag = 1e10))$gamma,
    error = function(e) NA_real_
  )
  if (!is.finite(C00) || C00 <= 0) return(-Inf)

  # Covariance matrix: C(h,u) = C(0,0) - gamma(h,u)
  Sigma <- C00 - Gamma
  # Diagonal is C(0,0) since gamma(0,0)=0
  diag(Sigma) <- C00

  tryCatch({
    L        <- t(chol(Sigma))                     # lower triangular, L L' = Sigma
    log_det  <- 2 * sum(log(diag(L)))

    ones     <- rep(1, n)
    z_ones   <- forwardsolve(L, ones)
    z_y      <- forwardsolve(L, y)

    mu_hat   <- sum(z_ones * z_y) / sum(z_ones^2)  # GLS mean (intercept model)
    r        <- y - mu_hat
    z_r      <- forwardsolve(L, r)
    quad     <- sum(z_r^2)

    -0.5 * (n * log(2 * pi) + log_det + quad)
  }, error = function(e) -Inf)
}
