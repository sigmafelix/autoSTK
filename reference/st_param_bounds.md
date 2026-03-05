# Compute parameter bounds for vgmST model optimisation

Compute parameter bounds for vgmST model optimisation

## Usage

``` r
st_param_bounds(
  model_template,
  stva_emp,
  sill_scale = 2,
  range_scale = 3,
  ani_scale = 20
)
```

## Arguments

- model_template:

  A `vgmST` object (initial guess).

- stva_emp:

  Empirical ST variogram (`StVariogram`/data.frame) from `setSTI` or
  [`gstat::variogramST`](https://r-spatial.github.io/gstat/reference/variogramST.html).

- sill_scale:

  Numeric (2). Upper bound multiplier applied to `max(stva_emp$gamma)`
  for sill/psill/nugget parameters.

- range_scale:

  Numeric (3). Upper bound multiplier applied to the maximum spatial lag
  for range parameters.

- ani_scale:

  Numeric (20). Controls the breadth of the stAni search interval
  relative to the data extent.

## Value

A named list with elements `lower` and `upper` (both named numeric
vectors matching `gstat::extractPar(model_template)`).
