# A convenience function for the sample spatiotemporal variogram

A convenience function for the sample spatiotemporal variogram

## Usage

``` r
setSTI(
  stf,
  formula,
  tlags = 0:6,
  cutoff = 30000,
  width = 1000,
  assumeRegular = FALSE,
  pseudo = TRUE,
  na.omit = TRUE,
  wireframe = FALSE,
  plot3d = FALSE,
  cores = 1
)
```

## Arguments

- stf:

  A ST\*DF object.

- formula:

  formula (inherits the same parameter in variogramST)

- tlags:

  temporal lags to compute semivariance (inherits the same parameter in
  variogramST)

- cutoff:

  the maximum bound of the set of spatial lags (inherits the same
  parameter in variogramST)

- width:

  integer (1). spatial lag (inherits the same parameter in variogramST)

- assumeRegular:

  Boolean. Assuming regular grid?

- pseudo:

  Boolean. See ?gstat::variogramST

- na.omit:

  Boolean. Omit NA values.

- wireframe:

  Boolean. Whether you plot a StVariogram in wireframe or not. If not,
  the return will be in class of data.frame, not a list

- plot3d:

  Boolean. Wheter you make a three-dimensional graph with rgl package

- logarithm:

  Boolean. log-transformation

## Value

Depends on the arguments wireframe (if TRUE, list of length 2) and
plot3d (if TRUE, list of length 3), a StVariogram object otherwise.
