# Fit a Spatio-Temporal Variogram

A named wrapper around
[`autofitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/autofitVariogramST.md)
that makes the fit/predict separation explicit. Returns a
`STVariogramFit` object which can be passed directly to
[`predictKrigeST`](https://sigmafelix.github.io/autoSTK/reference/predictKrigeST.md).

## Usage

``` r
fitVariogramST(formula, data, ...)
```

## Arguments

- formula:

  A formula (e.g. `PM10 ~ 1`).

- data:

  A spatio-temporal data object (STFDF, STSDF, STIDF, or sftime).

- ...:

  Arguments forwarded to
  [`autofitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/autofitVariogramST.md).

## Value

A `STVariogramFit` object.

## See also

[`predictKrigeST`](https://sigmafelix.github.io/autoSTK/reference/predictKrigeST.md),
[`autoKrigeST`](https://sigmafelix.github.io/autoSTK/reference/autoKrigeST.md)
