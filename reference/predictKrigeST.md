# Predict Using a Fitted Spatio-Temporal Kriging Model

Takes a `STVariogramFit` object (from
[`fitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/fitVariogramST.md)
or
[`autofitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/autofitVariogramST.md))
and performs spatiotemporal Kriging at the locations in `newdata`.

## Usage

``` r
predictKrigeST(
  fit,
  data,
  newdata,
  formula,
  nmax = Inf,
  predict_chunk = NULL,
  ...
)
```

## Arguments

- fit:

  A `STVariogramFit` object.

- data:

  The original data used to fit the variogram (STFDF, STSDF, STIDF, or
  sftime).

- newdata:

  Prediction grid (STFDF or sftime).

- formula:

  Formula used during fitting.

- nmax:

  Passed to
  [`gstat::krigeST`](https://r-spatial.github.io/gstat/reference/krigeST.html).

- predict_chunk:

  Optional integer; splits `newdata` into chunks for large grids to
  avoid memory issues.

- ...:

  Additional arguments forwarded to
  [`gstat::krigeST`](https://r-spatial.github.io/gstat/reference/krigeST.html).

## Value

An `autoKrigeST` object (list) with elements `krige_output` (STFDF) and
`var_model`.

## See also

[`fitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/fitVariogramST.md),
[`autoKrigeST`](https://sigmafelix.github.io/autoSTK/reference/autoKrigeST.md)
