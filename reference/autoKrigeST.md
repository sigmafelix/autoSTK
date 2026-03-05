# Automatic Spatio-Temporal Kriging

Convenience wrapper that fits the variogram and predicts in one call.
Internally calls
[`fitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/fitVariogramST.md)
then
[`predictKrigeST`](https://sigmafelix.github.io/autoSTK/reference/predictKrigeST.md).

## Usage

``` r
autoKrigeST(
  formula,
  input_data,
  new_data,
  type_stv = "sumMetric",
  block = 0,
  model = c("Sph", "Exp", "Gau", "Ste"),
  kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
  fix.values = c(NA, NA, NA),
  newdata_mode = "rect",
  newdata_npoints = 3000,
  GLS.model = NA,
  tlags = 0:6,
  cutoff = 20000,
  width = 500,
  forward = 6,
  predict_chunk = NULL,
  nmax = Inf,
  aniso_method = "vgm",
  type_joint = "Exp",
  prodsum_k = 0.25,
  surface = FALSE,
  start_vals = c(NA, NA, NA),
  miscFitOptions = list(),
  measurement_error = c(0, 0, 0),
  cores = 1L,
  verbose = TRUE,
  optimizer = "lbfgsb",
  objective = "WLS",
  n_restart = 1L,
  optimizer_control = list()
)
```

## Arguments

- formula:

  A formula (e.g. `PM10 ~ 1`).

- input_data:

  Training data (STFDF, STSDF, STIDF, or sftime).

- new_data:

  Prediction grid. If missing, generated automatically via
  [`create_new_data.ST`](https://sigmafelix.github.io/autoSTK/reference/create_new_data.ST.md).

- type_stv:

  Character. ST variogram model type. Passed as `typestv` to
  [`autofitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/autofitVariogramST.md).

- block:

  Numeric. Block size for block Kriging.

- model:

  Character vector. Candidate marginal variogram model families.

- kappa:

  Numeric vector. Kappa values for Matern-family models.

- fix.values:

  Numeric vector (3). Fixed initial values for nugget, range, sill.

- newdata_mode:

  Character. `"rect"` or `"chull"`.

- newdata_npoints:

  Integer. Number of prediction grid points.

- GLS.model:

  Variogram model for GLS sample variogram (rarely needed).

- tlags:

  Integer vector. Temporal lags.

- cutoff:

  Numeric. Maximum spatial lag.

- width:

  Numeric. Spatial lag bin width.

- forward:

  Integer. Number of future time steps to predict.

- predict_chunk:

  Integer or NULL. Chunk size for large prediction grids.

- nmax:

  Numeric. Maximum number of ST neighbours.

- aniso_method:

  Character. Anisotropy estimation method.

- type_joint:

  Character. Joint variogram component model family.

- prodsum_k:

  Numeric. k parameter for productSum models.

- surface:

  Logical. Return variogram surface?

- start_vals:

  Numeric vector (3). Not used in v2.0 (kept for API compat).

- miscFitOptions:

  List. Not used in v2.0 (kept for API compat).

- measurement_error:

  Numeric vector (3). Measurement error components.

- cores:

  Integer. Cores for variogramST computation.

- verbose:

  Logical. Print progress messages.

- optimizer:

  Character. `"lbfgsb"` or `"grid"`.

- objective:

  Character. `"WLS"` or `"MLE"`.

- n_restart:

  Integer. Number of optimisation restarts.

- optimizer_control:

  List. Extra control arguments for the optimiser.

## Value

An `autoKrigeST` object.

## Examples

``` r
data(air)
#> Warning: data set ‘air’ not found
deair <- STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
#> Error in STFDF(stations, dates, data.frame(PM10 = as.vector(air))): could not find function "STFDF"
deair_rs <- deair[, 3751:3800]
#> Error: object 'deair' not found
## Not run:
# akst <- autoKrigeST(formula = PM10 ~ 1, input_data = deair_rs,
#                     cutoff = 300000, width = 30000, tlags = 0:7, cores = 4)
```
