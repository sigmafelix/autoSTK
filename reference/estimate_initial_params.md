# Estimate Initial ST Variogram Parameters from Empirical Data

Derives nugget, partial sill, spatial range, and temporal range directly
from the empirical spatio-temporal variogram instead of relying on fixed
arithmetic heuristics. The approach is:

## Usage

``` r
estimate_initial_params(
  stva_emp,
  nugget_quantile = 0.25,
  sill_quantile = 0.95,
  practical_range_target = 0.95
)
```

## Arguments

- stva_emp:

  An empirical `StVariogram` (data.frame with columns `spacelag`,
  `timelag`, `gamma`, `np`).

- nugget_quantile:

  Numeric in (0, 1). Quantile of near-zero-lag gamma used as the nugget
  estimate. Default 0.25 (lower quartile of the first decile of spatial
  bins).

- sill_quantile:

  Numeric in (0, 1). Quantile of the spatial-marginal gamma values used
  as the plateau estimate. Default 0.95.

- practical_range_target:

  Numeric in (0, 1). Fraction of the sill at which the practical range
  is read off. Default 0.95.

## Value

A named list:

- nugget:

  Estimated nugget (\>= 0).

- psill:

  Estimated partial sill (sill - nugget, \>= 0).

- sill:

  Estimated total sill.

- sp_range:

  Estimated spatial practical range.

- ts_range:

  Estimated temporal practical range.

- max_gamma:

  Maximum observed gamma.

- max_splag:

  Maximum spatial lag in the empirical variogram.

- max_tlag:

  Maximum temporal lag (numeric).

## Details

1.  **Nugget** — median \\\gamma\\ over the first decile of spatial lags
    (extrapolation toward zero-lag).

2.  **Sill** — 95th percentile of \\\gamma\\ values over all spatial
    lags at the minimum time lag (spatial marginal plateau).

3.  **Spatial practical range** — smallest spatial lag at which the
    spatial marginal exceeds 95% of the estimated sill; falls back to
    50% of the maximum spatial lag when the variogram never reaches the
    sill within the sampling window.

4.  **Temporal practical range** — same rule applied to the temporal
    marginal at the minimum spatial lag.
