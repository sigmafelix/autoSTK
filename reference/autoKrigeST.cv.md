# Cross-validation of spatiotemporal Kriging

Cross-validation of spatiotemporal Kriging

## Usage

``` r
autoKrigeST.cv(
  data,
  fold_dim = c("spatial", "temporal", "random", "spacetime"),
  nfold = 10L,
  formula,
  type_stv = "sumMetric",
  block = 0,
  model = c("Sph", "Exp", "Gau", "Ste"),
  kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
  fix.values = c(NA, NA, NA),
  tlags = 0:6,
  cutoff = 20000,
  width = 500,
  nmax = Inf,
  aniso_method = "vgm",
  type_joint = "Exp",
  prodsum_k = 0.25,
  surface = FALSE,
  start_vals = c(NA, NA, NA),
  miscFitOptions = list(),
  measurement_error = c(0, 0, 0),
  cores = 1,
  seed = 130425L,
  variogram_from_full = FALSE,
  optimizer = "lbfgsb",
  objective = "WLS",
  n_restart = 1L,
  optimizer_control = list()
)
```

## Arguments

- data:

  a \`STFDF\`-class object

- fold_dim:

  character. the dimension at which you want to cross-validate (spatial,
  temporal, and random)

- nfold:

  integer. the number of folds. 10 as the default

- formula:

  formula. e.g., y~1

- type_stv:

  character. One of 'sumMetric', 'metric', 'productSum', and 'separable'

- block:

  numeric. passed to conduct block spatiotemporal Kriging.

- model:

  character vector. Default is c("Sph", "Exp", "Gau", "Ste"), but users
  can specify the list of theoretical variograms by referring
  gstat::vgm.

- kappa:

  numeric vector. Kappa values tested for Matern-family variogram
  models.

- fix.values:

  numeric vector. Initial values in order of nugget, range, and sill,
  respectively.

- tlags:

  integer vector (increasing, preferably to be consecutive). temporal
  lags.

- cutoff:

  numeric. The maximum distance at which the sample variogram will be
  computed.

- width:

  numeric. The interval at which the variogram cloud will be summarized.

- nmax:

  integer or positive infinite. The maximum number of spatiotemporal
  neighbors to conduct the local spatiotemporal Kriging.

- aniso_method:

  character. One of 'vgm', 'linear', 'range', and 'metric'. Please refer
  to ?gstat::estiStAni.

- type_joint:

  character. The model form of joint spatiotemporal variogram.

- prodsum_k:

  numeric. The parameter for the case when 'productSum' is chosen for
  type_stv.

- surface:

  logical. If TRUE, also return variogram surface in fitting path.

- start_vals:

  numeric vector (3). The initial values to optimize the spatiotemporal
  variogram model.

- miscFitOptions:

  list. See ?automap::autofitVariogram.

- measurement_error:

  numeric vector of length 3: spatial, temporal, and joint measurement
  error terms.

- cores:

  integer. The number of threads that will be used to compute the sample
  spatiotemporal variogram.

- seed:

  integer. Random seed for fold generation.

- variogram_from_full:

  logical. If TRUE, fit one ST variogram on full data and reuse across
  folds.

- optimizer:

  character. Optimization method passed to ST variogram fitting.

- objective:

  character. Objective for ST variogram fitting ('WLS' or 'MLE').

- n_restart:

  integer. Number of restarts for optimizer = 'lbfgsb'.

- optimizer_control:

  list. Extra control options forwarded to the selected optimizer.

## Value

The cross-validated spatiotemporal Kriging results.

## Examples

``` r
# \donttest{
library(sp)
library(gstat)
library(spacetime)
library(stars)
#> Loading required package: abind
#> Loading required package: sf
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
data(air)
deair <- STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
deair_sf <- st_as_stars(deair, crs = "+proj=longlat +ellps=sphere")
deair_sf <- st_transform(deair_sf, 3857)
deair_r <- as(deair_sf, "STFDF")
deair_r@sp@proj4string <- CRS("EPSG:3857")
deair_rs <- deair_r[, 3751:3800]
## autoKrigeST.cv test
if (FALSE) { # \dontrun{
akst_cv_t <- autoKrigeST.cv(
  formula = PM10 ~ 1, data = deair_rs, nfold = 3, fold_dim = "temporal",
  cutoff = 300000, width = 30000, tlags = 0:7, cores = 8
)
akst_cv_s <- autoKrigeST.cv(
  formula = PM10 ~ 1, data = deair_rs, nfold = 3, fold_dim = "spatial",
  cutoff = 300000, width = 30000, tlags = 0:7, cores = 8
)
# akst_cv_r = autoKrigeST.cv(formula = PM10~1, data = deair_rs,  nfold = 3, fold_dim = 'random',
#                          cutoff = 300000, width = 30000, tlags = 0:7, cores = 8)
akst_cv_spt <- autoKrigeST.cv(
  formula = PM10 ~ 1, data = deair_rs, nfold = 4, fold_dim = "spacetime",
  cutoff = 300000, width = 30000, tlags = 0:7, cores = 8
)
} # }
# }
```
