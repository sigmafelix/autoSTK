# Likelihood Ratio Test for Spatiotemporal Separability

Tests whether a separable covariance structure is adequate versus the
more general sumMetric model using a likelihood ratio test (LRT). Both
models are fitted with `objective = "MLE"`.

## Usage

``` r
test_separability(stf, formula, ...)
```

## Arguments

- stf:

  Spatio-temporal data (STFDF, STSDF, STIDF, or sftime).

- formula:

  Formula, e.g. `PM10 ~ 1`.

- ...:

  Additional arguments forwarded to `autofitVariogramST` (e.g. `tlags`,
  `cutoff`, `width`, `cores`).

## Value

A list with elements:

- statistic:

  LRT chi-squared statistic.

- df:

  Degrees of freedom (difference in number of parameters).

- p_value:

  p-value from chi-squared distribution.

- separable_fit:

  Fitted `STVariogramFit` for the separable model.

- summetric_fit:

  Fitted `STVariogramFit` for the sumMetric model.
