# Compute AIC and BIC for a fitted STVariogramFit (MLE-based)

Compute AIC and BIC for a fitted STVariogramFit (MLE-based)

## Usage

``` r
ic_stv(fit_obj, n_obs)
```

## Arguments

- fit_obj:

  A `STVariogramFit` object with an `objective = "MLE"` fit (i.e.,
  `fit_obj$loglik` must not be `NULL`).

- n_obs:

  Integer. Number of observations used for fitting.

## Value

A named list: `AIC`, `BIC`, `k` (number of free parameters), `loglik`.
