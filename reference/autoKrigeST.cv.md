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

- start_vals:

  numeric vector (3). The initial values to optimize the spatiotemporal
  variogram model.

- miscFitOptions:

  list. See ?automap::autofitVariogram.

- cores:

  integer. The number of threads that will be used to compute the sample
  spatiotemporal variogram.

- newdata_mode:

  character. One of 'rect' (rectangular grid) and 'chull' (convex hull)

- newdata_npoints:

  integer. The number of points that will be generated in the range of
  geometry the user specified (one of rectangle or convex hull)

- GLS.model:

  a variogram model. The default value is NA. If a variogram model is
  passed, a Generalized Lease Squares sample variogram will be
  calculated.

- predict_chunk:

  integer. The number of data points per chunk in the new data for the
  large data. It should be meticulously chosen according to the user's
  machine specification.

## Value

The cross-validated spatiotemporal Kriging results.

## Examples

``` r
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
akst_cv_t <- autoKrigeST.cv(
  formula = PM10 ~ 1, data = deair_rs, nfold = 3, fold_dim = "temporal",
  cutoff = 300000, width = 30000, tlags = 0:7, cores = 8
)
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> [[1]]
#>   model     psill   range
#> 1   Nug -2.750512     0.0
#> 2   Exp 17.024991 34616.6
#> 
#> [[2]]
#>   model     psill    range kappa
#> 1   Nug -66.90675     0.00  0.00
#> 2   Ste  81.39284 26171.65  0.05
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -31.04160     0.00   0.0
#> 2   Ste  45.48776 34018.16   0.1
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -13.20526     0.00   0.0
#> 2   Ste  27.58955 41755.95   0.2
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -7.334976     0.00   0.0
#> 2   Ste 21.673364 45536.86   0.3
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -4.447594     0.00   0.0
#> 2   Ste 18.750677 47671.01   0.4
#> 
#> [[7]]
#>   model     psill    range kappa
#> 1   Nug -2.748901     0.00   0.0
#> 2   Ste 17.023821 48964.78   0.5
#> 
#> [[8]]
#>   model     psill    range kappa
#> 1   Nug -1.641193     0.00   0.0
#> 2   Ste 15.893040 49781.52   0.6
#> 
#> [[9]]
#>   model      psill    range kappa
#> 1   Nug -0.8713713     0.00   0.0
#> 2   Ste 15.1030390 50289.48   0.7
#> 
#> [[10]]
#>   model      psill    range kappa
#> 1   Nug -0.3068278     0.00   0.0
#> 2   Ste 14.5218746 50627.03   0.8
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model       psill    range
#> 1   Nug  0.03479162     0.00
#> 2   Sph 13.81354645 91346.26
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 1            Sph     0 1.200893e-06
#> 16           Ste    10 1.268739e-06
#> 15           Ste     5 1.291266e-06
#> 14           Ste     2 1.336985e-06
#> 13           Ste   1.9 1.339898e-06
#> 12           Ste   1.8 1.342981e-06
#> 11           Ste   1.7 1.346251e-06
#> 10           Ste   1.6 1.349723e-06
#> 9            Ste   1.5 1.353415e-06
#> 8            Ste   1.4 1.357347e-06
#> 7            Ste   1.3 1.361541e-06
#> 6            Ste   1.2 1.366021e-06
#> 5            Ste   1.1 1.370816e-06
#> 4            Ste     1 1.375955e-06
#> 3            Ste   0.9 1.381474e-06
#> 2            Gau     0 1.405466e-06
#> [[1]]
#>   model      psill    range
#> 1   Nug -0.2573231 0.000000
#> 2   Exp 26.9824013 1.815439
#> 
#> [[2]]
#>   model     psill    range kappa
#> 1   Nug -96.84206 0.000000  0.00
#> 2   Ste 124.96292 1.631926  0.05
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -42.86098 0.000000   0.0
#> 2   Ste  70.70258 2.035699   0.1
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -16.01864 0.000000   0.0
#> 2   Ste  43.44711 2.365654   0.2
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -7.178896 0.000000   0.0
#> 2   Ste 34.310855 2.491334   0.3
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -2.825176 0.000000   0.0
#> 2   Ste 29.730585 2.545152   0.4
#> 
#> [[7]]
#>   model     psill    range kappa
#> 1   Nug -0.257322 0.000000   0.0
#> 2   Ste 26.982401 2.567418   0.5
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model     psill   range kappa
#> 1   Nug  4.638563 0.00000     0
#> 2   Ste 21.537599 2.54947     1
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 7            Ste     1 4.156670e+00
#> 6            Ste   0.9 4.197673e+00
#> 8            Ste   1.1 4.467457e+00
#> 5            Ste   0.8 4.685713e+00
#> 9            Ste   1.2 5.054896e+00
#> 4            Ste   0.7 5.742546e+00
#> 10           Ste   1.3 5.859299e+00
#> 11           Ste   1.4 6.832977e+00
#> 3            Ste   0.6 7.525431e+00
#> 12           Ste   1.5 7.937644e+00
#> 13           Ste   1.6 9.142452e+00
#> 14           Ste   1.7 1.042246e+01
#> 15           Ste   1.8 1.175747e+01
#> 16           Ste   1.9 1.313110e+01
#> 17           Ste     2 1.453004e+01
#> 18           Ste     5 4.910915e+01
#> 2            Gau     0 2.085264e+05
#> 19           Ste    10 2.102906e+05
#> 1            Sph     0 3.041361e+05
#> Warning: singular model in variogram fit
#> Initial estimates — nugget: 4.078  psill: 13.11  sp_range: 2.85e+05  ts_range: 5
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> Warning: strictly irregular time steps were assumed to be regular
#> [[1]]
#>   model     psill   range
#> 1   Nug -16.25538     0.0
#> 2   Sph  30.56212 38792.4
#> 
#> [[2]]
#>   model     psill    range
#> 1   Nug -2.393876     0.00
#> 2   Exp 22.535408 61550.85
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -64.60273     0.00  0.00
#> 2   Ste  87.52757 89310.82  0.05
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -29.54177     0.00   0.0
#> 2   Ste  51.65907 95755.73   0.1
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -12.24803     0.00   0.0
#> 2   Ste  33.47536 95232.43   0.2
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -6.650159     0.00   0.0
#> 2   Ste 27.371457 92336.67   0.3
#> 
#> [[7]]
#>   model     psill    range kappa
#> 1   Nug -3.949852     0.00   0.0
#> 2   Ste 24.334894 89514.55   0.4
#> 
#> [[8]]
#>   model     psill   range kappa
#> 1   Nug -2.393734     0.0   0.0
#> 2   Ste 22.535390 87048.1   0.5
#> 
#> [[9]]
#>   model     psill    range kappa
#> 1   Nug -1.400336     0.00   0.0
#> 2   Ste 21.355961 84930.65   0.6
#> 
#> [[10]]
#>   model      psill    range kappa
#> 1   Nug -0.7222315     0.00   0.0
#> 2   Ste 20.5301341 83110.12   0.7
#> 
#> [[11]]
#>   model      psill    range kappa
#> 1   Nug -0.2368332     0.00   0.0
#> 2   Ste 19.9240495 81534.37   0.8
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model     psill    range kappa
#> 1   Nug  1.417093     0.00     0
#> 2   Ste 17.623935 71944.83     2
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 13           Ste     2 2.298537e-06
#> 12           Ste   1.9 2.298799e-06
#> 11           Ste   1.8 2.299154e-06
#> 10           Ste   1.7 2.299619e-06
#> 9            Ste   1.6 2.300214e-06
#> 8            Ste   1.5 2.300964e-06
#> 7            Ste   1.4 2.301902e-06
#> 6            Ste   1.3 2.303065e-06
#> 14           Ste     5 2.303427e-06
#> 5            Ste   1.2 2.304502e-06
#> 4            Ste   1.1 2.306274e-06
#> 3            Ste     1 2.308460e-06
#> 2            Ste   0.9 2.311158e-06
#> 15           Ste    10 2.312606e-06
#> 1            Gau     0 2.475958e-06
#> [[1]]
#>   model     psill    range
#> 1   Nug -12.83443 0.000000
#> 2   Exp  41.31969 1.153024
#> 
#> [[2]]
#>   model     psill     range kappa
#> 1   Nug -231.4688 0.0000000  0.00
#> 2   Ste  260.1569 0.6962493  0.05
#> 
#> [[3]]
#>   model     psill     range kappa
#> 1   Nug -109.9179 0.0000000   0.0
#> 2   Ste  138.5765 0.9442824   0.1
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -49.19882 0.000000   0.0
#> 2   Ste  77.80413 1.237524   0.2
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -28.98493 0.000000   0.0
#> 2   Ste  57.54428 1.416767   0.3
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -18.88809 0.000000   0.0
#> 2   Ste  47.40774 1.540123   0.4
#> 
#> [[7]]
#>   model     psill    range kappa
#> 1   Nug -12.80455 0.000000   0.0
#> 2   Ste  41.29167 1.631681   0.5
#> 
#> [[8]]
#>   model     psill    range kappa
#> 1   Nug -8.795074 0.000000   0.0
#> 2   Ste 37.250798 1.700137   0.6
#> 
#> [[9]]
#>   model     psill    range kappa
#> 1   Nug -5.942383 0.000000   0.0
#> 2   Ste 34.369849 1.753627   0.7
#> 
#> [[10]]
#>   model     psill    range kappa
#> 1   Nug -3.788924 0.000000   0.0
#> 2   Ste 32.192739 1.797518   0.8
#> 
#> [[11]]
#>   model     psill    range kappa
#> 1   Nug -2.152007 0.000000   0.0
#> 2   Ste 30.531573 1.831682   0.9
#> 
#> [[12]]
#>   model      psill    range kappa
#> 1   Nug -0.8582684 0.000000     0
#> 2   Ste 29.2147616 1.859259     1
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model      psill    range kappa
#> 1   Nug  0.2325107 0.000000   0.0
#> 2   Ste 28.1066306 1.884528   1.1
#> 
#> Tested models, best first:
#>    Tested.models kappa     SSerror
#> 3            Ste   1.1    458.3484
#> 4            Ste   1.2    460.8384
#> 5            Ste   1.3    463.1737
#> 6            Ste   1.4    465.3684
#> 7            Ste   1.5    467.4353
#> 8            Ste   1.6    469.3852
#> 9            Ste   1.7    471.2279
#> 10           Ste   1.8    472.9722
#> 11           Ste   1.9    474.6257
#> 12           Ste     2    476.1955
#> 13           Ste     5    502.3187
#> 14           Ste    10    516.7612
#> 2            Gau     0  68466.9666
#> 1            Sph     0 139928.8370
#> Warning: No convergence after 200 iterations: try different initial values?
#> Initial estimates — nugget: 3.213  psill: 19.61  sp_range: 2.85e+05  ts_range: 10.39
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> [[1]]
#>   model     psill    range
#> 1   Nug -25.69519     0.00
#> 2   Sph  39.08181 33264.26
#> 
#> [[2]]
#>   model     psill    range
#> 1   Nug -3.093181     0.00
#> 2   Exp 23.762645 72576.72
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -63.03713      0.0  0.00
#> 2   Ste  91.17363 205666.8  0.05
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -29.15582      0.0   0.0
#> 2   Ste  54.05569 155600.2   0.1
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -12.50099      0.0   0.0
#> 2   Ste  35.18308 127871.7   0.2
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -7.138641      0.0   0.0
#> 2   Ste 28.821557 115628.6   0.3
#> 
#> [[7]]
#>   model     psill    range kappa
#> 1   Nug -4.567709      0.0   0.0
#> 2   Ste 25.647678 108002.1   0.4
#> 
#> [[8]]
#>   model     psill    range kappa
#> 1   Nug -3.093243      0.0   0.0
#> 2   Ste 23.762634 102637.8   0.5
#> 
#> [[9]]
#>   model     psill    range kappa
#> 1   Nug -2.159182     0.00   0.0
#> 2   Ste 22.523996 98533.57   0.6
#> 
#> [[10]]
#>   model     psill    range kappa
#> 1   Nug -1.523342     0.00   0.0
#> 2   Ste 21.655084 95319.22   0.7
#> 
#> [[11]]
#>   model     psill    range kappa
#> 1   Nug -1.071365     0.00   0.0
#> 2   Ste 21.015716 92680.59   0.8
#> 
#> [[12]]
#>   model      psill    range kappa
#> 1   Nug -0.7378885     0.00   0.0
#> 2   Ste 20.5283368 90477.51   0.9
#> 
#> [[13]]
#>   model      psill    range kappa
#> 1   Nug -0.4847955     0.00     0
#> 2   Ste 20.1463069 88606.91     1
#> 
#> [[14]]
#>   model      psill    range kappa
#> 1   Nug -0.2892328     0.00   0.0
#> 2   Ste 19.8398311 86980.45   1.1
#> 
#> [[15]]
#>   model      psill    range kappa
#> 1   Nug -0.1332368     0.00   0.0
#> 2   Ste 19.5896943 85587.19   1.2
#> 
#> [[16]]
#>   model        psill    range kappa
#> 1   Nug -0.008212176     0.00   0.0
#> 2   Ste 19.381992341 84354.94   1.3
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model      psill    range kappa
#> 1   Nug  0.0936723     0.00   0.0
#> 2   Ste 19.2072174 83262.23   1.4
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 2            Ste   1.4 1.651075e-06
#> 3            Ste   1.5 1.657897e-06
#> 4            Ste   1.6 1.664345e-06
#> 5            Ste   1.7 1.670439e-06
#> 6            Ste   1.8 1.676201e-06
#> 7            Ste   1.9 1.681651e-06
#> 8            Ste     2 1.686811e-06
#> 9            Ste     5 1.768422e-06
#> 10           Ste    10 1.809209e-06
#> 1            Gau     0 1.861869e-06
#> [[1]]
#>   model    psill    range kappa
#> 1   Nug -59.1403   0.0000  0.00
#> 2   Ste 151.9497 172.7039  0.05
#> 
#> [[2]]
#>   model     psill   range kappa
#> 1   Nug -20.96228   0.000   0.0
#> 2   Ste 155.46550 428.539   0.1
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -1.113918 0.000000   0.0
#> 2   Ste 46.552538 6.939322   0.2
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model     psill    range kappa
#> 1   Nug  5.356541 0.000000   0.0
#> 2   Ste 36.424222 5.285908   0.3
#> 
#> Tested models, best first:
#>    Tested.models kappa    SSerror
#> 4            Ste   0.3   2371.533
#> 5            Ste   0.4   2475.144
#> 6            Ste   0.5   2568.451
#> 2            Exp     0   2568.451
#> 7            Ste   0.6   2653.277
#> 8            Ste   0.7   2730.873
#> 9            Ste   0.8   2802.182
#> 10           Ste   0.9   2867.955
#> 11           Ste     1   2928.811
#> 12           Ste   1.1   2985.270
#> 13           Ste   1.2   3037.776
#> 14           Ste   1.3   3086.707
#> 15           Ste   1.4   3132.406
#> 16           Ste   1.5   3175.162
#> 17           Ste   1.6   3215.234
#> 18           Ste   1.7   3252.852
#> 19           Ste   1.8   3288.226
#> 20           Ste   1.9   3321.531
#> 21           Ste     2   3352.938
#> 22           Ste     5   3835.264
#> 23           Ste    10   4057.689
#> 3            Gau     0 400179.935
#> 1            Sph     0 602544.390
#> Initial estimates — nugget: 1.984  psill: 21.84  sp_range: 2.85e+05  ts_range: 6
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
akst_cv_s <- autoKrigeST.cv(
  formula = PM10 ~ 1, data = deair_rs, nfold = 3, fold_dim = "spatial",
  cutoff = 300000, width = 30000, tlags = 0:7, cores = 8
)
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> [[1]]
#>   model      psill    range
#> 1   Nug -0.8276907      0.0
#> 2   Sph 18.6071037 136638.3
#> 
#> [[2]]
#>   model    psill    range
#> 1   Nug -2.11424     0.00
#> 2   Exp 21.15690 61623.69
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -55.58551     0.00  0.00
#> 2   Ste  76.71769 88365.37  0.05
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -25.25458     0.00   0.0
#> 2   Ste  45.77216 95245.12   0.1
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -10.37835     0.00   0.0
#> 2   Ste  30.22669 95097.59   0.2
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -5.634511     0.00   0.0
#> 2   Ste 25.099996 92229.57   0.3
#> 
#> [[7]]
#>   model     psill    range kappa
#> 1   Nug -3.385326     0.00   0.0
#> 2   Ste 22.602586 89486.91   0.4
#> 
#> [[8]]
#>   model     psill    range kappa
#> 1   Nug -2.115449     0.00   0.0
#> 2   Ste 21.155860 87116.49   0.5
#> 
#> [[9]]
#>   model     psill    range kappa
#> 1   Nug -1.323256     0.00   0.0
#> 2   Ste 20.229071 85089.87   0.6
#> 
#> [[10]]
#>   model      psill    range kappa
#> 1   Nug -0.7950941     0.00   0.0
#> 2   Ste 19.5949448 83367.91   0.7
#> 
#> [[11]]
#>   model      psill    range kappa
#> 1   Nug -0.4260886     0.00   0.0
#> 2   Ste 19.1396680 81891.11   0.8
#> 
#> [[12]]
#>   model     psill    range kappa
#> 1   Nug -0.158891     0.00   0.0
#> 2   Ste 18.800553 80614.61   0.9
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model      psill    range kappa
#> 1   Nug  0.9752227     0.00     0
#> 2   Ste 16.8849618 66392.92    10
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 14           Ste    10 8.971262e-07
#> 13           Ste     5 9.165445e-07
#> 12           Ste     2 9.820790e-07
#> 11           Ste   1.9 9.876806e-07
#> 10           Ste   1.8 9.938275e-07
#> 9            Ste   1.7 1.000594e-06
#> 8            Ste   1.6 1.008069e-06
#> 7            Ste   1.5 1.016355e-06
#> 6            Ste   1.4 1.025576e-06
#> 5            Ste   1.3 1.035880e-06
#> 4            Ste   1.2 1.047442e-06
#> 3            Ste   1.1 1.060479e-06
#> 2            Ste     1 1.075248e-06
#> 1            Gau     0 1.257860e-06
#> [[1]]
#>   model      psill    range kappa
#> 1   Nug -0.5192606  0.00000  0.00
#> 2   Ste 39.7243332 45.91883  0.05
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model    psill    range kappa
#> 1   Nug  5.64060  0.00000   0.0
#> 2   Ste 39.68648 12.39721   0.3
#> 
#> Tested models, best first:
#>    Tested.models kappa     SSerror
#> 6            Ste   0.3    370.4376
#> 7            Ste   0.4    403.7349
#> 8            Ste   0.5    434.3052
#> 2            Exp     0    434.3052
#> 9            Ste   0.6    462.6071
#> 10           Ste   0.7    488.8888
#> 11           Ste   0.8    513.3336
#> 12           Ste   0.9    536.0952
#> 13           Ste     1    557.3098
#> 14           Ste   1.1    577.1004
#> 15           Ste   1.2    595.5794
#> 16           Ste   1.3    612.8496
#> 17           Ste   1.4    629.0051
#> 18           Ste   1.5    644.1327
#> 19           Ste   1.6    658.3114
#> 20           Ste   1.7    671.6141
#> 21           Ste   1.8    684.1072
#> 22           Ste   1.9    695.8514
#> 23           Ste     2    706.9031
#> 24           Ste     5    869.9700
#> 5            Ste   0.2   1414.5593
#> 25           Ste    10  10385.0381
#> 4            Ste   0.1  52644.2061
#> 3            Gau     0 229397.8256
#> 1            Sph     0 331960.0448
#> Initial estimates — nugget: 1.571  psill: 20.17  sp_range: 2.85e+05  ts_range: 6
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> [[1]]
#>   model     psill    range
#> 1   Nug -2.615519     0.00
#> 2   Exp 20.730707 54632.28
#> 
#> [[2]]
#>   model     psill    range kappa
#> 1   Nug -63.12143     0.00  0.00
#> 2   Ste  83.00943 66380.39  0.05
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -29.07991     0.00   0.0
#> 2   Ste  48.51694 76001.24   0.1
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -12.25866     0.00   0.0
#> 2   Ste  31.13703 80124.36   0.2
#> 
#> [[5]]
#>   model     psill   range kappa
#> 1   Nug -6.794077     0.0   0.0
#> 2   Ste 25.327641 79830.4   0.3
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -4.147447     0.00   0.0
#> 2   Ste 22.441278 78626.77   0.4
#> 
#> [[7]]
#>   model     psill    range kappa
#> 1   Nug -2.615513     0.00   0.0
#> 2   Ste 20.730705 77261.79   0.5
#> 
#> [[8]]
#>   model     psill    range kappa
#> 1   Nug -1.633027     0.00   0.0
#> 2   Ste 19.608820 75942.04   0.6
#> 
#> [[9]]
#>   model      psill    range kappa
#> 1   Nug -0.9591942     0.00   0.0
#> 2   Ste 18.8225508 74725.77   0.7
#> 
#> [[10]]
#>   model      psill    range kappa
#> 1   Nug -0.4745252     0.00   0.0
#> 2   Ste 18.2449115 73623.87   0.8
#> 
#> [[11]]
#>   model      psill    range kappa
#> 1   Nug -0.1132527     0.00   0.0
#> 2   Ste 17.8052583 72631.09   0.9
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model     psill    range kappa
#> 1   Nug  1.652072     0.00     0
#> 2   Ste 15.182624 59275.66    10
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 15           Ste    10 2.663324e-06
#> 14           Ste     5 2.668975e-06
#> 13           Ste     2 2.688794e-06
#> 12           Ste   1.9 2.690378e-06
#> 11           Ste   1.8 2.692096e-06
#> 10           Ste   1.7 2.693966e-06
#> 9            Ste   1.6 2.696005e-06
#> 8            Ste   1.5 2.698237e-06
#> 7            Ste   1.4 2.700689e-06
#> 6            Ste   1.3 2.703391e-06
#> 5            Ste   1.2 2.706384e-06
#> 4            Ste   1.1 2.709715e-06
#> 3            Ste     1 2.713442e-06
#> 2            Gau     0 2.721971e-06
#> 1            Sph     0 1.229503e-04
#> [[1]]
#>   model      psill    range
#> 1   Nug -0.1765302 0.000000
#> 2   Exp 30.0561499 1.341714
#> 
#> [[2]]
#>   model     psill    range kappa
#> 1   Nug -120.6635 0.000000  0.00
#> 2   Ste  151.4580 1.039794  0.05
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -53.50513 0.000000   0.0
#> 2   Ste  84.13547 1.343326   0.1
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -19.98833 0.000000   0.0
#> 2   Ste  50.36466 1.637831   0.2
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -8.915837 0.000000   0.0
#> 2   Ste 39.091119 1.776766   0.3
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -3.442554 0.000000   0.0
#> 2   Ste 33.453016 1.851744   0.4
#> 
#> [[7]]
#>   model      psill    range kappa
#> 1   Nug -0.1750564 0.000000   0.0
#> 2   Ste 30.0551240 1.897632   0.5
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model     psill    range kappa
#> 1   Nug  1.969491 0.000000   0.0
#> 2   Ste 27.797669 1.924807   0.6
#> 
#> Tested models, best first:
#>    Tested.models kappa    SSerror
#> 3            Ste   0.6   1263.208
#> 4            Ste   0.7   1295.431
#> 5            Ste   0.8   1325.210
#> 6            Ste   0.9   1352.810
#> 7            Ste     1   1378.462
#> 8            Ste   1.1   1402.362
#> 9            Ste   1.2   1424.680
#> 10           Ste   1.3   1445.568
#> 11           Ste   1.4   1465.156
#> 12           Ste   1.5   1483.559
#> 13           Ste   1.6   1500.880
#> 14           Ste   1.7   1517.212
#> 15           Ste   1.8   1532.634
#> 16           Ste   1.9   1547.218
#> 17           Ste     2   1561.032
#> 18           Ste     5   1784.856
#> 19           Ste    10   1902.103
#> 2            Gau     0 301561.583
#> 1            Sph     0 466295.065
#> Warning: No convergence after 200 iterations: try different initial values?
#> Initial estimates — nugget: 3.102  psill: 18.28  sp_range: 2.85e+05  ts_range: 6
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> [[1]]
#>   model     psill   range kappa
#> 1   Nug -40.87672       0  0.00
#> 2   Ste  84.84012 3199734  0.05
#> 
#> [[2]]
#>   model     psill   range kappa
#> 1   Nug -16.11544       0   0.0
#> 2   Ste  82.73422 9807864   0.1
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -5.300322      0.0   0.0
#> 2   Ste 27.261806 185703.8   0.2
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -1.706565      0.0   0.0
#> 2   Ste 21.833933 142349.5   0.3
#> 
#> [[5]]
#>   model       psill    range kappa
#> 1   Nug -0.01900666      0.0   0.0
#> 2   Ste 19.23248338 123154.5   0.4
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model      psill    range kappa
#> 1   Nug  0.9046112      0.0   0.0
#> 2   Ste 17.7248541 111272.2   0.5
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 4            Ste   0.5 2.540922e-06
#> 2            Exp     0 2.540927e-06
#> 5            Ste   0.6 2.578568e-06
#> 6            Ste   0.7 2.611059e-06
#> 7            Ste   0.8 2.639130e-06
#> 8            Ste   0.9 2.663403e-06
#> 9            Ste     1 2.684441e-06
#> 10           Ste   1.1 2.702709e-06
#> 11           Ste   1.2 2.718650e-06
#> 12           Ste   1.3 2.732603e-06
#> 13           Ste   1.4 2.744884e-06
#> 14           Ste   1.5 2.755750e-06
#> 15           Ste   1.6 2.765396e-06
#> 16           Ste   1.7 2.774007e-06
#> 17           Ste   1.8 2.781734e-06
#> 18           Ste   1.9 2.788684e-06
#> 19           Ste     2 2.794965e-06
#> 20           Ste     5 2.865410e-06
#> 21           Ste    10 2.886436e-06
#> 3            Gau     0 3.065001e-06
#> 1            Sph     0 6.280397e-05
#> [[1]]
#>   model     psill    range kappa
#> 1   Nug -102.4619 0.000000  0.00
#> 2   Ste  134.6849 1.299998  0.05
#> 
#> [[2]]
#>   model     psill    range kappa
#> 1   Nug -43.74052 0.000000   0.0
#> 2   Ste  75.70389 1.639029   0.1
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -14.50830 0.000000   0.0
#> 2   Ste  46.07718 1.930503   0.2
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -4.851016 0.000000   0.0
#> 2   Ste 36.131337 2.050458   0.3
#> 
#> [[5]]
#>   model       psill    range kappa
#> 1   Nug -0.06727243 0.000000   0.0
#> 2   Ste 31.12755033 2.108462   0.4
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model     psill    range
#> 1   Nug  2.744822 0.000000
#> 2   Exp 28.129338 1.508827
#> 
#> Tested models, best first:
#>    Tested.models kappa    SSerror
#> 2            Exp     0   1472.275
#> 4            Ste   0.5   1472.275
#> 5            Ste   0.6   1512.310
#> 6            Ste   0.7   1549.144
#> 7            Ste   0.8   1583.145
#> 8            Ste   0.9   1614.625
#> 9            Ste     1   1643.849
#> 10           Ste   1.1   1671.050
#> 11           Ste   1.2   1696.424
#> 12           Ste   1.3   1720.146
#> 13           Ste   1.4   1742.369
#> 14           Ste   1.5   1763.230
#> 15           Ste   1.6   1782.841
#> 16           Ste   1.7   1801.314
#> 17           Ste   1.8   1818.742
#> 18           Ste   1.9   1835.209
#> 19           Ste     2   1850.789
#> 20           Ste     5   2100.833
#> 21           Ste    10   2229.502
#> 3            Gau     0 269163.756
#> 1            Sph     0 419211.366
#> Warning: No convergence after 200 iterations: try different initial values?
#> Initial estimates — nugget: 4.767  psill: 16.91  sp_range: 2.85e+05  ts_range: 6
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
# akst_cv_r = autoKrigeST.cv(formula = PM10~1, data = deair_rs,  nfold = 3, fold_dim = 'random',
#                          cutoff = 300000, width = 30000, tlags = 0:7, cores = 8)
akst_cv_spt <- autoKrigeST.cv(
  formula = PM10 ~ 1, data = deair_rs, nfold = 4, fold_dim = "spacetime",
  cutoff = 300000, width = 30000, tlags = 0:7, cores = 8
)
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> [[1]]
#>   model    psill    range
#> 1   Nug -3.65287     0.00
#> 2   Exp 19.43348 42769.85
#> 
#> [[2]]
#>   model     psill    range kappa
#> 1   Nug -65.16195     0.00  0.00
#> 2   Ste  80.94663 34265.97  0.05
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -30.43573     0.00   0.0
#> 2   Ste  46.22595 44236.55   0.1
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -13.34534     0.00   0.0
#> 2   Ste  29.12807 53343.58   0.2
#> 
#> [[5]]
#>   model    psill    range kappa
#> 1   Nug -7.84018     0.00   0.0
#> 2   Ste 23.60870 57207.87   0.3
#> 
#> [[6]]
#>   model     psill   range kappa
#> 1   Nug -5.181789     0.0   0.0
#> 2   Ste 20.952421 59279.1   0.4
#> 
#> [[7]]
#>   model     psill    range kappa
#> 1   Nug -3.659542     0.00   0.0
#> 2   Ste 19.433313 60390.13   0.5
#> 
#> [[8]]
#>   model     psill    range kappa
#> 1   Nug -2.695907     0.00   0.0
#> 2   Ste 18.472747 60985.16   0.6
#> 
#> [[9]]
#>   model     psill    range kappa
#> 1   Nug -2.039653     0.00   0.0
#> 2   Ste 17.825584 61371.03   0.7
#> 
#> [[10]]
#>   model     psill    range kappa
#> 1   Nug -1.579721     0.00   0.0
#> 2   Ste 17.366293 61476.91   0.8
#> 
#> [[11]]
#>   model    psill    range kappa
#> 1   Nug -1.24213     0.00   0.0
#> 2   Ste 17.02915 61482.77   0.9
#> 
#> [[12]]
#>   model      psill    range kappa
#> 1   Nug -0.9873308     0.00     0
#> 2   Ste 16.7745542 61427.96     1
#> 
#> [[13]]
#>   model      psill    range kappa
#> 1   Nug -0.7935451     0.00   0.0
#> 2   Ste 16.5750723 61269.06   1.1
#> 
#> [[14]]
#>   model      psill    range kappa
#> 1   Nug -0.6375449     0.00   0.0
#> 2   Ste 16.4207306 61180.21   1.2
#> 
#> [[15]]
#>   model      psill    range kappa
#> 1   Nug -0.5126778     0.00   0.0
#> 2   Ste 16.2971506 61076.84   1.3
#> 
#> [[16]]
#>   model      psill   range kappa
#> 1   Nug -0.4121481     0.0   0.0
#> 2   Ste 16.1956309 60943.4   1.4
#> 
#> [[17]]
#>   model      psill    range kappa
#> 1   Nug -0.3294466     0.00   0.0
#> 2   Ste 16.1118807 60809.57   1.5
#> 
#> [[18]]
#>   model      psill    range kappa
#> 1   Nug -0.2586308     0.00   0.0
#> 2   Ste 16.0447870 60731.11   1.6
#> 
#> [[19]]
#>   model      psill    range kappa
#> 1   Nug -0.2007314     0.00   0.0
#> 2   Ste 15.9860056 60604.82   1.7
#> 
#> [[20]]
#>   model      psill    range kappa
#> 1   Nug -0.1533175     0.00   0.0
#> 2   Ste 15.9331901 60434.09   1.8
#> 
#> [[21]]
#>   model     psill    range kappa
#> 1   Nug -0.110983     0.00   0.0
#> 2   Ste 15.890440 60322.84   1.9
#> 
#> [[22]]
#>   model       psill    range kappa
#> 1   Nug -0.07445916     0.00     0
#> 2   Ste 15.85338245 60215.33     2
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model      psill    range kappa
#> 1   Nug  0.3042431     0.00     0
#> 2   Ste 15.4390120 57631.89    10
#> 
#> Tested models, best first:
#>   Tested.models kappa      SSerror
#> 4           Ste    10 7.792503e-07
#> 3           Ste     5 7.926989e-07
#> 2           Gau     0 8.609380e-07
#> 1           Sph     0 4.496065e-05
#> Selected:
#>   model     psill    range kappa
#> 1   Nug  5.231168 0.000000   0.0
#> 2   Ste 21.116614 5.869454   1.2
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 16           Ste   1.2     19.21318
#> 15           Ste   1.1     19.26208
#> 17           Ste   1.3     19.28715
#> 18           Ste   1.4     19.45189
#> 14           Ste     1     19.47666
#> 19           Ste   1.5     19.68317
#> 13           Ste   0.9     19.91483
#> 20           Ste   1.6     19.96252
#> 21           Ste   1.7     20.27584
#> 22           Ste   1.8     20.61231
#> 12           Ste   0.8     20.65633
#> 23           Ste   1.9     20.96363
#> 24           Ste     2     21.32343
#> 11           Ste   0.7     21.81374
#> 10           Ste   0.6     23.55240
#> 9            Ste   0.5     26.13190
#> 2            Exp     0     26.13190
#> 8            Ste   0.4     30.02168
#> 7            Ste   0.3   6583.71122
#> 6            Ste   0.2  33527.28298
#> 5            Ste   0.1  55029.32626
#> 4            Ste  0.05  78005.35490
#> 25           Ste     5  89092.69779
#> 3            Gau     0  90282.05015
#> 26           Ste    10  90750.94820
#> 1            Sph     0 121122.75117
#> Initial estimates — nugget: 0.9947  psill: 19.72  sp_range: 1.35e+05  ts_range: 7
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> [[1]]
#>   model     psill   range
#> 1   Nug -1.455498     0.0
#> 2   Exp 24.424008 75394.7
#> 
#> [[2]]
#>   model     psill    range kappa
#> 1   Nug -56.21026      0.0  0.00
#> 2   Ste  86.20723 233738.8  0.05
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -24.98261      0.0   0.0
#> 2   Ste  51.83605 170452.8   0.1
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -9.789865      0.0   0.0
#> 2   Ste 34.519328 135651.2   0.2
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -4.980515      0.0   0.0
#> 2   Ste 28.808045 121181.6   0.3
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -2.720539      0.0   0.0
#> 2   Ste 26.029799 112546.4   0.4
#> 
#> [[7]]
#>   model     psill  range kappa
#> 1   Nug -1.455845      0   0.0
#> 2   Ste 24.423320 106611   0.5
#> 
#> [[8]]
#>   model      psill    range kappa
#> 1   Nug -0.6708335      0.0   0.0
#> 2   Ste 23.4021321 102311.6   0.6
#> 
#> [[9]]
#>   model      psill   range kappa
#> 1   Nug -0.1564081     0.0   0.0
#> 2   Ste 22.6957934 98803.6   0.7
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model     psill    range kappa
#> 1   Nug  1.483245     0.00     0
#> 2   Ste 19.640590 72586.15    10
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 17           Ste    10 1.774686e-06
#> 16           Ste     5 1.798760e-06
#> 15           Ste     2 1.874137e-06
#> 14           Ste   1.9 1.880333e-06
#> 13           Ste   1.8 1.887099e-06
#> 12           Ste   1.7 1.894509e-06
#> 11           Ste   1.6 1.902651e-06
#> 10           Ste   1.5 1.911625e-06
#> 9            Ste   1.4 1.921555e-06
#> 8            Ste   1.3 1.932580e-06
#> 7            Ste   1.2 1.944873e-06
#> 6            Ste   1.1 1.958640e-06
#> 5            Ste     1 1.974128e-06
#> 4            Ste   0.9 1.991639e-06
#> 3            Ste   0.8 2.011544e-06
#> 2            Gau     0 2.092623e-06
#> 1            Sph     0 4.671628e-05
#> [[1]]
#>   model     psill    range kappa
#> 1   Nug -63.34609 0.000000  0.00
#> 2   Ste 102.26997 2.303536  0.05
#> 
#> [[2]]
#>   model     psill    range kappa
#> 1   Nug -20.71496 0.000000   0.0
#> 2   Ste  58.95728 2.641531   0.1
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model      psill    range kappa
#> 1   Nug  0.4627998 0.000000   0.0
#> 2   Ste 36.9288911 2.804719   0.2
#> 
#> Tested models, best first:
#>    Tested.models kappa     SSerror
#> 4            Ste   0.2    382.7498
#> 5            Ste   0.3    397.8473
#> 6            Ste   0.4    411.5777
#> 7            Ste   0.5    424.1450
#> 2            Exp     0    424.1451
#> 8            Ste   0.6    435.7020
#> 9            Ste   0.7    446.3699
#> 10           Ste   0.8    456.2483
#> 11           Ste   0.9    465.4211
#> 12           Ste     1    473.9600
#> 13           Ste   1.1    481.9267
#> 14           Ste   1.2    489.3753
#> 15           Ste   1.3    496.3528
#> 16           Ste   1.4    502.9010
#> 17           Ste   1.5    509.0571
#> 18           Ste   1.6    514.8540
#> 19           Ste   1.7    520.3199
#> 20           Ste   1.8    525.4825
#> 21           Ste   1.9    530.3648
#> 22           Ste     2    534.9884
#> 23           Ste     5    609.2719
#> 24           Ste    10    646.8962
#> 3            Gau     0  69973.4425
#> 1            Sph     0 118841.3350
#> Warning: No convergence after 200 iterations: try different initial values?
#> Initial estimates — nugget: 2.147  psill: 23.44  sp_range: 1.35e+05  ts_range: 6
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> [[1]]
#>   model      psill    range
#> 1   Nug -0.3478508     0.00
#> 2   Exp 14.4241077 30592.15
#> 
#> [[2]]
#>   model     psill    range kappa
#> 1   Nug -60.65205     0.00  0.00
#> 2   Ste  74.85623 21536.34  0.05
#> 
#> [[3]]
#>   model     psill    range kappa
#> 1   Nug -27.10682     0.00   0.0
#> 2   Ste  41.28572 28344.67   0.1
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -10.35474     0.00   0.0
#> 2   Ste  24.49475 35485.16   0.2
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -4.767371     0.00   0.0
#> 2   Ste 18.880830 39345.13   0.3
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -2.020112     0.00   0.0
#> 2   Ste 16.110473 41628.15   0.4
#> 
#> [[7]]
#>   model      psill    range kappa
#> 1   Nug -0.3488305     0.00   0.0
#> 2   Ste 14.4249836 43260.36   0.5
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model    psill    range kappa
#> 1   Nug 4.972598     0.00     0
#> 2   Ste 8.953094 45412.34    10
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 19           Ste    10 2.545107e-06
#> 18           Ste     5 2.561715e-06
#> 2            Gau     0 2.569578e-06
#> 17           Ste     2 2.592650e-06
#> 16           Ste   1.9 2.594435e-06
#> 15           Ste   1.8 2.596290e-06
#> 14           Ste   1.7 2.598232e-06
#> 13           Ste   1.6 2.600272e-06
#> 12           Ste   1.5 2.602399e-06
#> 11           Ste   1.4 2.604633e-06
#> 10           Ste   1.3 2.606971e-06
#> 9            Ste   1.2 2.609429e-06
#> 8            Ste   1.1 2.612011e-06
#> 7            Ste     1 2.614712e-06
#> 6            Ste   0.9 2.617552e-06
#> 5            Ste   0.8 2.620528e-06
#> 4            Ste   0.7 2.623652e-06
#> 3            Ste   0.6 2.626925e-06
#> 1            Sph     0 1.904228e-05
#> [[1]]
#>   model     psill    range
#> 1   Nug -6.660689 0.000000
#> 2   Exp 29.857051 1.060647
#> 
#> [[2]]
#>   model     psill     range kappa
#> 1   Nug -152.7338 0.0000000  0.00
#> 2   Ste  176.0640 0.6612364  0.05
#> 
#> [[3]]
#>   model     psill     range kappa
#> 1   Nug -71.39884 0.0000000   0.0
#> 2   Ste  94.70947 0.8922805   0.1
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -30.79993 0.000000   0.0
#> 2   Ste  54.07638 1.159837   0.2
#> 
#> [[5]]
#>   model    psill    range kappa
#> 1   Nug -17.3353 0.000000   0.0
#> 2   Ste  40.5815 1.318528   0.3
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -10.64572 0.000000   0.0
#> 2   Ste  33.86533 1.424511   0.4
#> 
#> [[7]]
#>   model    psill   range kappa
#> 1   Nug -6.65965 0.00000   0.0
#> 2   Ste 29.85618 1.50005   0.5
#> 
#> [[8]]
#>   model    psill    range kappa
#> 1   Nug -4.02521 0.000000   0.0
#> 2   Ste 27.20115 1.556023   0.6
#> 
#> [[9]]
#>   model     psill    range kappa
#> 1   Nug -2.165168 0.000000   0.0
#> 2   Ste 25.321795 1.598408   0.7
#> 
#> [[10]]
#>   model      psill    range kappa
#> 1   Nug -0.7801943 0.000000   0.0
#> 2   Ste 23.9203097 1.631827   0.8
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model     psill    range kappa
#> 1   Nug  4.342118 0.000000   0.0
#> 2   Ste 18.685775 1.759063   1.8
#> 
#> Tested models, best first:
#>    Tested.models kappa     SSerror
#> 12           Ste   1.8    53.05329
#> 11           Ste   1.7    54.08952
#> 10           Ste   1.6    55.19776
#> 9            Ste   1.5    56.38534
#> 8            Ste   1.4    57.66061
#> 7            Ste   1.3    59.03308
#> 6            Ste   1.2    60.51366
#> 5            Ste   1.1    62.11491
#> 4            Ste     1    63.85136
#> 3            Ste   0.9    65.73988
#> 15           Ste     5 27481.80129
#> 14           Ste     2 44238.40992
#> 13           Ste   1.9 45700.02406
#> 16           Ste    10 48836.12247
#> 2            Gau     0 48883.18752
#> 1            Sph     0 75959.57404
#> Warning: No convergence after 200 iterations: try different initial values?
#> Initial estimates — nugget: 7.901  psill: 12.35  sp_range: 7.5e+04  ts_range: 3
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: 'tzone' attributes are inconsistent
#> Warning: 'tzone' attributes are inconsistent
#> Fitting the optimal spatiotemporal variogram model...
#> [[1]]
#>   model     psill    range
#> 1   Nug -2.365506      0.0
#> 2   Exp 25.642446 100601.9
#> 
#> [[2]]
#>   model      psill    range
#> 1   Nug -0.7446184     0.00
#> 2   Gau 18.1369191 58192.72
#> 
#> [[3]]
#>   model     psill   range kappa
#> 1   Nug -55.99046       0  0.00
#> 2   Ste 111.18933 3258554  0.05
#> 
#> [[4]]
#>   model     psill    range kappa
#> 1   Nug -28.95534        0   0.0
#> 2   Ste 138.80712 20051140   0.1
#> 
#> [[5]]
#>   model     psill    range kappa
#> 1   Nug -12.28544      0.0   0.0
#> 2   Ste  40.24831 219584.3   0.2
#> 
#> [[6]]
#>   model     psill    range kappa
#> 1   Nug -6.691436      0.0   0.0
#> 2   Ste 32.032168 173153.3   0.3
#> 
#> [[7]]
#>   model     psill    range kappa
#> 1   Nug -3.957002      0.0   0.0
#> 2   Ste 28.027382 153717.2   0.4
#> 
#> [[8]]
#>   model     psill    range kappa
#> 1   Nug -2.371743      0.0   0.0
#> 2   Ste 25.640268 142125.6   0.5
#> 
#> [[9]]
#>   model    psill    range kappa
#> 1   Nug -1.33796      0.0   0.0
#> 2   Ste 24.06098 134508.8   0.6
#> 
#> [[10]]
#>   model      psill    range kappa
#> 1   Nug -0.6212813      0.0   0.0
#> 2   Ste 22.9389204 128934.3   0.7
#> 
#> [[11]]
#>   model      psill    range kappa
#> 1   Nug -0.0888399      0.0   0.0
#> 2   Ste 22.1043501 124853.4   0.8
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model      psill    range kappa
#> 1   Nug  0.3099989      0.0   0.0
#> 2   Ste 21.4580865 121507.8   0.9
#> 
#> Tested models, best first:
#>    Tested.models kappa      SSerror
#> 2            Ste   0.9 1.070373e-06
#> 3            Ste     1 1.094446e-06
#> 4            Ste   1.1 1.116504e-06
#> 5            Ste   1.2 1.136749e-06
#> 6            Ste   1.3 1.155366e-06
#> 7            Ste   1.4 1.172515e-06
#> 8            Ste   1.5 1.188338e-06
#> 9            Ste   1.6 1.202957e-06
#> 10           Ste   1.7 1.216493e-06
#> 11           Ste   1.8 1.229043e-06
#> 12           Ste   1.9 1.240691e-06
#> 13           Ste     2 1.251528e-06
#> 14           Ste     5 1.396046e-06
#> 15           Ste    10 1.447883e-06
#> 1            Sph     0 5.222607e-05
#> [[1]]
#>   model     psill     range
#> 1   Nug -35.04191 0.0000000
#> 2   Exp  71.22018 0.6303269
#> 
#> [[2]]
#>   model     psill     range kappa
#> 1   Nug -377.8581 0.0000000  0.00
#> 2   Ste  414.5584 0.4058415  0.05
#> 
#> [[3]]
#>   model     psill     range kappa
#> 1   Nug -188.0038 0.0000000   0.0
#> 2   Ste  224.6148 0.5428775   0.1
#> 
#> [[4]]
#>   model     psill     range kappa
#> 1   Nug -93.40076 0.0000000   0.0
#> 2   Ste 129.85523 0.6950429   0.2
#> 
#> [[5]]
#>   model     psill     range kappa
#> 1   Nug -61.98218 0.0000000   0.0
#> 2   Ste  98.30685 0.7812225   0.3
#> 
#> [[6]]
#>   model     psill     range kappa
#> 1   Nug -45.58019 0.0000000   0.0
#> 2   Ste  81.81503 0.8420977   0.4
#> 
#> [[7]]
#>   model     psill     range kappa
#> 1   Nug -35.73414 0.0000000   0.0
#> 2   Ste  71.89191 0.8854161   0.5
#> 
#> [[8]]
#>   model     psill     range kappa
#> 1   Nug -29.00777 0.0000000   0.0
#> 2   Ste  65.10368 0.9191955   0.6
#> 
#> [[9]]
#>   model     psill     range kappa
#> 1   Nug -24.36686 0.0000000   0.0
#> 2   Ste  60.40251 0.9435495   0.7
#> 
#> [[10]]
#>   model     psill     range kappa
#> 1   Nug -20.72106 0.0000000   0.0
#> 2   Ste  56.70931 0.9645869   0.8
#> 
#> [[11]]
#>   model     psill     range kappa
#> 1   Nug -17.84717 0.0000000   0.0
#> 2   Ste  53.79436 0.9820887   0.9
#> 
#> [[12]]
#>   model     psill     range kappa
#> 1   Nug -15.46975 0.0000000     0
#> 2   Ste  51.38295 0.9975484     1
#> 
#> [[13]]
#>   model     psill    range kappa
#> 1   Nug -13.28114 0.000000   0.0
#> 2   Ste  49.17345 1.013893   1.1
#> 
#> [[14]]
#>   model     psill   range kappa
#> 1   Nug -11.73814 0.00000   0.0
#> 2   Ste  47.59916 1.02396   1.2
#> 
#> [[15]]
#>   model     psill    range kappa
#> 1   Nug -10.36175 0.000000   0.0
#> 2   Ste  46.19721 1.033602   1.3
#> 
#> [[16]]
#>   model     psill    range kappa
#> 1   Nug -9.273051 0.000000   0.0
#> 2   Ste 45.080900 1.040568   1.4
#> 
#> [[17]]
#>   model    psill    range kappa
#> 1   Nug -8.46608 0.000000   0.0
#> 2   Ste 44.24221 1.044439   1.5
#> 
#> [[18]]
#>   model     psill    range kappa
#> 1   Nug -7.598676 0.000000   0.0
#> 2   Ste 43.353099 1.050305   1.6
#> 
#> [[19]]
#>   model     psill    range kappa
#> 1   Nug -6.811088 0.000000   0.0
#> 2   Ste 42.546413 1.055847   1.7
#> 
#> [[20]]
#>   model     psill    range kappa
#> 1   Nug -5.862324 0.000000   0.0
#> 2   Ste 41.591944 1.064986   1.8
#> 
#> [[21]]
#>   model     psill    range kappa
#> 1   Nug -5.281685 0.000000   0.0
#> 2   Ste 40.992649 1.068693   1.9
#> 
#> [[22]]
#>   model     psill    range kappa
#> 1   Nug -4.713347 0.000000     0
#> 2   Ste 40.409159 1.072804     2
#> 
#> ^^^ ABOVE MODELS WERE REMOVED ^^^
#> 
#> Warning: Some models where removed for being either NULL or having a negative sill/range/nugget, 
#>  set verbose == TRUE for more information
#> Selected:
#>   model     psill    range kappa
#> 1   Nug  1.888992 0.000000     0
#> 2   Ste 33.592457 1.121365     5
#> 
#> Tested models, best first:
#>   Tested.models kappa    SSerror
#> 3           Ste     5   2241.981
#> 4           Ste    10   2253.273
#> 2           Gau     0 101422.240
#> 1           Sph     0 251078.727
#> Initial estimates — nugget: 1.349  psill: 24.19  sp_range: 2.25e+05  ts_range: 5
#> Warning: All optimisation attempts failed; returning initial model template.
#> Predicting 6 time step(s)...
#> Warning: The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.
#>  The unit 'days' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
```
