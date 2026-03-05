# Automated Spatio-Temporal Covariance Model Selection

Fits multiple ST covariance model types and ranks them by a user-chosen
criterion.

## Usage

``` r
selectModelST(
  stf,
  formula,
  candidates = c("separable", "productSum", "sumMetric", "simpleSumMetric", "metric"),
  criterion = c("WLS", "AIC", "BIC"),
  optimizer = "lbfgsb",
  cores = 1L,
  ...
)
```

## Arguments

- stf:

  A spatio-temporal data object (STFDF, STSDF, STIDF, or sftime).

- formula:

  A formula (e.g., `PM10 ~ 1`).

- candidates:

  Character vector of ST model types to try. Defaults to all six
  supported types.

- criterion:

  One of `"AIC"`, `"BIC"` (require `objective = "MLE"`), `"CV_RMSE"`
  (not yet implemented in v2.0; use `"WLS"` as a fast proxy), or
  `"WLS"`.

- optimizer:

  Passed to `autofitVariogramST`.

- cores:

  Passed to `autofitVariogramST`.

- ...:

  Additional arguments forwarded to `autofitVariogramST`.

## Value

An `STModelSelection` object (list) with elements `best_model`,
`comparison` (data.frame), and `all_fits`.
