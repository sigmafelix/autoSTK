# Automatic Spatio-Temporal Variogram Fitting

Fits a spatio-temporal variogram model to a spatio-temporal data frame
using automatic parameter estimation. Supports multiple optimisers and
objective functions (v2.0+).

## Usage

``` r
autofitVariogramST(
  stf,
  formula,
  typestv = "sumMetric",
  candidate_model = c("Ste", "Exc", "Exp", "Wav"),
  guess_nugget = NULL,
  guess_psill = NULL,
  tlags = 0:6,
  cutoff = 20000,
  width = 500,
  aniso_method = "vgm",
  type_joint = "Exp",
  prodsum_k = NULL,
  surface = FALSE,
  measurement_error = c(0, 0, 0),
  cores = 1L,
  verbose = FALSE,
  optimizer = c("lbfgsb", "grid"),
  objective = c("WLS", "MLE"),
  n_restart = 1L,
  optimizer_control = list()
)
```

## Arguments

- stf:

  A spatio-temporal data frame (STFDF, STSDF, STIDF, or sftime).

- formula:

  A formula specifying the response and explanatory variables.

- typestv:

  Character. ST variogram model type. One of `"sumMetric"`,
  `"separable"`, `"productSum"`, `"productSumOld"`, `"simpleSumMetric"`,
  `"metric"`.

- candidate_model:

  Character vector of candidate spatial/temporal variogram model
  families (e.g. `c("Ste", "Exc", "Exp", "Wav")`).

- guess_nugget:

  Optional initial guess for the nugget. Estimated if NULL.

- guess_psill:

  Optional initial guess for the partial sill. Estimated if NULL.

- tlags:

  Integer vector of temporal lags.

- cutoff:

  Numeric. Maximum spatial lag distance.

- width:

  Numeric. Width of spatial lag bins.

- aniso_method:

  Character. Method for estimating ST anisotropy ratio (`"vgm"` or
  `"linear"`).

- type_joint:

  Character. Model family for the joint variogram component.

- prodsum_k:

  Numeric. k parameter for productSum models.

- surface:

  Logical. If TRUE, also compute and return the variogram surface.

- measurement_error:

  Numeric vector of length 3: spatial, temporal, joint.

- cores:

  Integer. Number of cores for `variogramST`.

- verbose:

  Logical. If TRUE, print diagnostic messages.

- optimizer:

  Character. Optimisation strategy: `"lbfgsb"` (default, current
  behaviour) or `"grid"` (LHS grid search + L-BFGS-B).

- objective:

  Character. Fitting criterion: `"WLS"` (weighted least squares,
  default) or `"MLE"` (Gaussian log-likelihood; only for n \< 500).

- n_restart:

  Integer. Number of optimisation restarts for `optimizer = "lbfgsb"`.
  Default 1 matches v1.x behaviour.

- optimizer_control:

  Named list of extra arguments forwarded to the optimiser (e.g.
  `list(n_coarse = 100L)` for `"grid"`).

## Value

A `STVariogramFit` object (list) with elements:

- jointSTV:

  Fitted ST variogram model.

- empSTV:

  Empirical ST variogram.

- SpV:

  Fitted spatial marginal variogram.

- TV:

  Fitted temporal marginal variogram.

- STVsurface:

  (Optional) Variogram surface if `surface = TRUE`.

- optimizer:

  Character. Optimiser used.

- objective:

  Character. Objective function used.

- loglik:

  Numeric. Log-likelihood (only when `objective = "MLE"`).

- n_obs:

  Integer. Number of observations (for AIC/BIC).

## See also

[`fitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/fitVariogramST.md),
[`selectModelST`](https://sigmafelix.github.io/autoSTK/reference/selectModelST.md),
[`test_separability`](https://sigmafelix.github.io/autoSTK/reference/test_separability.md)

## Examples

``` r
library(spacetime)
library(sp)
library(sftime)
data(air)
stations <- st_as_sf(stations)
stations <- st_transform(stations, 'EPSG:3857')
airdf <- data.frame(PM10 = as.vector(air))
stations_full <- do.call(c, rep(stations, length(dates)))
dates_full <- rep(dates, each = nrow(stations))
rural <- cbind(airdf, time = dates_full, stations_full)
rural <- st_as_sftime(rural, sf_column_name = 'geometry')
rr <- rural[match(rural$time, dates[3001:3060], nomatch = FALSE) > 0, ]
rr <- as(as(as(rr, "STIDF"), "STFDF"), "STSDF")
rrstv <- autofitVariogramST(stf = rr, formula = PM10 ~ 1, surface = TRUE)
#> Error in apply(do.call(cbind, lapply(ret, function(x) x$np)), 1, sum,     na.rm = TRUE): dim(X) must have a positive length
rrstv
#> Error: object 'rrstv' not found
```
