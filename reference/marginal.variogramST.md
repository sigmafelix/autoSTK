# Compute the marginal spatial or temporal sample variogram

Compute the marginal spatial or temporal sample variogram

## Usage

``` r
marginal.variogramST(stv, bound, spatial = TRUE)
```

## Arguments

- stv:

  A StVariogram and data.frame object.

- bound:

  numeric. The maximum distance that will be used to compute a spatial
  variogram.

- spatial:

  boolean. if TRUE, the spatial marginal variogram will be obtained,
  temporal otherwise.

## Value

A gstatVariogram object.
