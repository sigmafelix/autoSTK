# Sensitivity Plot of the Spatio-Temporal Anisotropy Ratio

Evaluates
[`gstat::estiStAni`](https://r-spatial.github.io/gstat/reference/estiStAni.html)
at a sequence of spatial intervals and plots how the estimated
anisotropy ratio varies. Useful for diagnosing whether the
single-interval estimate used in `autofitVariogramST` is stable.

## Usage

``` r
plot_aniso_sensitivity(
  stva_emp,
  spatial_vgm,
  temporal_vgm,
  n_intervals = 20L,
  plot = TRUE
)
```

## Arguments

- stva_emp:

  Empirical ST variogram from `setSTI`.

- spatial_vgm:

  Fitted spatial marginal variogram model.

- temporal_vgm:

  Fitted temporal marginal variogram model.

- n_intervals:

  Integer. Number of evaluation intervals.

- plot:

  Logical. Whether to produce a base-R plot.

## Value

Invisibly, a data.frame with columns `interval` and `stAni`.
