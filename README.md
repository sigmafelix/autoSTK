# `autoSTK`

`autoSTK` provides automatic spatio-temporal Kriging inspired by
`automap` (Hiemstra et al. 2010). It extends the `automap` workflow to
spatio-temporal variogram fitting, Kriging prediction, and cross-validation
for common R spatio-temporal data classes.

## Installation

``` r
remotes::install_github("sigmafelix/autoSTK")
```

Before installing, make sure the core spatial dependencies are available:
`gstat`, `spacetime`, `sp`, `sf`, `sftime`, and `automap`. The `cubble`
package is suggested and is only required when working directly with
`cubble_df` objects.

## Main Features

- Automatic spatio-temporal variogram fitting with `autofitVariogramST()`.
- Automatic spatio-temporal Kriging with `autoKrigeST()`.
- Explicit fit/predict workflow with `fitVariogramST()` and
  `predictKrigeST()`.
- Cross-validation with `autoKrigeST.cv()` using spatial, temporal,
  spatio-temporal, or random folds.
- Support for `STFDF`, `STSDF`, `STIDF`, `sftime`, and `cubble` inputs in the
  main fitting and prediction workflow.
- Conversion helper `cubble_to_sftime()` for nested or temporal `cubble`
  objects.
- Multiple spatio-temporal variogram structures, including `sumMetric`,
  `metric`, `productSum`, and `separable`.
- Multiple optimizers, including `lbfgsb`, `grid`, `sa`, and `ga`.

## Supported Functions

- `autofitVariogramST()`: automatically fits a spatio-temporal variogram.
- `fitVariogramST()`: fits a variogram and returns a reusable
  `STVariogramFit` object.
- `predictKrigeST()`: predicts from a fitted `STVariogramFit` object.
- `autoKrigeST()`: fits the variogram and predicts in one call.
- `autoKrigeST.cv()`: cross-validates spatio-temporal Kriging results.
- `cubble_to_sftime()`: converts `cubble_df` objects to `sftime`.

## Data Classes

The main Kriging workflow accepts both legacy `spacetime` classes and modern
simple-feature classes:

| Input class | Typical use |
| --- | --- |
| `STFDF` | Full space-time grids, such as all stations observed at all dates. |
| `STSDF` | Sparse space-time grids, such as missing station-date observations. |
| `STIDF` | Irregular space-time observations. |
| `sftime` | Simple-feature workflows with one geometry and one time column. |
| `cubble_df` | Tidy spatio-temporal data stored in spatial and temporal faces. |

Internally, `autoSTK` converts modern inputs to the `spacetime` classes
required by `gstat::variogramST()` and `gstat::krigeST()`.

## Use Cases

### Forecast Air Quality From Station Networks

Use `autoKrigeST()` when you have repeated measurements at monitoring
stations and want predictions on a future spatio-temporal grid. This is the
classic `STFDF`/`STSDF` workflow.

``` r
library(autoSTK)
library(gstat)
library(sp)
library(spacetime)
library(stars)

data(air)

deair <- STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
deair_sf <- st_as_stars(deair) |>
  st_transform("+proj=longlat +ellps=sphere") |>
  st_transform(3857)

deair_r <- as(deair_sf, "STFDF")
deair_r@sp@proj4string <- CRS(
  "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs +type=crs"
)

deair_rs <- deair_r[, 3701:3800]
deair_rss <- as(deair_rs, "STSDF")

akst_stk <- autoKrigeST(
  formula = PM10 ~ 1,
  input_data = deair_rss,
  cutoff = 300000,
  width = 30000,
  type_stv = "sumMetric",
  model = c("Exc", "Mat", "Ste", "Exp", "Wav"),
  tlags = 0:7,
  cores = 8,
  optimizer = "sa",
  verbose = TRUE
)

akst_stk_stars <- st_as_stars(akst_stk$krige_output)
plot(akst_stk_stars[1, ])
```

![](https://i.imgur.com/xDadzlJ.png)

### Work Directly With `sftime`

Use `sftime` when your data already lives in an `sf`-style workflow with one
geometry column and one time column. `autoKrigeST()` accepts `sftime` for both
training data and prediction data.

``` r
library(autoSTK)
library(sf)
library(sftime)

# station
data(air)
stations_sf <- sf::st_as_sf(stations) |>
  sf::st_transform(3857)

# expand to station-time records and keep 2009-05-10 to 2009-05-21
time_idx <- dates >= as.Date("2009-05-10") & dates <= as.Date("2009-05-21")
dates_sub <- dates[time_idx]
air_sub <- air[, time_idx, drop = FALSE]

station_long <- do.call(rbind, lapply(seq_along(dates_sub), function(i) {
  out <- stations_sf
  out$id <- seq_len(nrow(stations_sf))
  out$time <- dates_sub[i]
  out$PM10 <- as.numeric(air_sub[, i])
  out
}))
station_long <- station_long[!is.na(station_long$PM10), ]

# station_long has columns:
#   id, time, PM10, geometry
# where geometry is an sfc_POINT column and time is Date/POSIXct.
station_sf <- st_as_sf(
  station_long,
  sf_column_name = "geometry",
  crs = 3857
)

station_st <- st_as_sftime(
  station_sf,
  time_column_name = "time",
  sf_column_name = "geometry"
)

fit <- fitVariogramST(
  formula = PM10 ~ 1,
  data = station_st,
  typestv = "sumMetric",
  candidate_model = c("Sph", "Exp", "Gau", "Ste"),
  cutoff = 300000,
  width = 30000,
  tlags = 0:7
)

prediction_grid_st <-
  {
    # 1) spatial grid over the station_st extent
    station_sf <- sf::st_as_sf(station_st)
    bb <- sf::st_bbox(station_sf)
    bbox_poly <- sf::st_as_sfc(bb)
    grid_xy <- sf::st_make_grid(bbox_poly, cellsize = 20000, what = "centers")
    grid_sf <- sf::st_sf(grid_id = seq_along(grid_xy), geometry = grid_xy)

    # 2) expand grid across all times in station_st and fill with 1
    times <- sort(unique(sftime::st_time(station_st)))
    grid_long <- do.call(rbind, lapply(times, function(tt) {
      out <- grid_sf
      out$time <- tt
      out$const1 <- 1
      out
    }))

    # 3) convert to sftime for predictKrigeST(newdata = ...)
    sftime::st_as_sftime(grid_long, time_column_name = "time")
  }


pred <- predictKrigeST(
  fit = fit,
  data = station_st,
  newdata = prediction_grid_st,
  formula = PM10 ~ 1,
  nmax = 40,
  cores = 8L
)

pred_stars <- st_as_stars(pred$krige_output)
plot(pred_stars[3, ])

```

`fitVariogramST()` and `predictKrigeST()` are useful when you want to reuse the
same fitted variogram across several prediction grids.

### Use `cubble` Data

Use `cubble` when your workflow stores station metadata in the spatial face
and repeated measurements in the temporal face. `autoSTK` can convert either
the nested or temporal face to `sftime`.

``` r
library(autoSTK)
library(cubble)

# Build a cubble object from the long sf table created above.
cb <- cubble::as_cubble(station_long, key = id, index = time)

# Convert nested/temporal cubble faces to sftime.
pm10_st <- cubble_to_sftime(cb, time_col = "time", key_col = "id")

akst_cb <- autoKrigeST(
  formula = PM10 ~ 1,
  input_data = pm10_st,
  cutoff = 300000,
  width = 30000,
  tlags = 0:7,
  type_stv = "sumMetric"
)
```

You can also pass a `cubble_df` directly to the main Kriging workflow when
supplying a prediction grid:

``` r
akst_cb <- autoKrigeST(
  formula = PM10 ~ 1,
  input_data = cb,
  new_data = prediction_grid_st,
  cutoff = 300000,
  width = 30000,
  tlags = 0:7
)
```

For cross-validation, convert `cubble` or `sftime` data to a `spacetime`
object first, because `autoKrigeST.cv()` builds folds from the `@sp` and
`@time` slots.

``` r
pm10_stfdf <- as(as(pm10_st, "STIDF"), "STFDF")

akst_cv_cb <- autoKrigeST.cv(
  formula = PM10 ~ 1,
  data = pm10_stfdf,
  nfold = 3,
  fold_dim = "temporal",
  cutoff = 300000,
  width = 30000,
  tlags = 0:7
)
```

## Cross-Validation

`autoKrigeST.cv()` evaluates prediction performance by holding out different
parts of the spatio-temporal data set.

| `fold_dim` | What is held out | Good for |
| --- | --- | --- |
| `"spatial"` | Groups of monitoring locations | Testing transfer to unseen stations. |
| `"temporal"` | Groups of dates/times | Testing forecasting over unseen periods. |
| `"spacetime"` | Blocks of locations and times | Testing combined spatial and temporal extrapolation. |
| `"random"` | Random observations | Testing interpolation under scattered missingness. |

``` r
akst_cv_t <- autoKrigeST.cv(
  formula = PM10 ~ 1,
  data = deair_rs,
  nfold = 3,
  fold_dim = "temporal",
  cutoff = 300000,
  width = 30000,
  tlags = 0:7,
  cores = 8,
  optimizer = "sa"
)

akst_cv_s <- autoKrigeST.cv(
  formula = PM10 ~ 1,
  data = deair_rs,
  nfold = 3,
  fold_dim = "spatial",
  cutoff = 300000,
  width = 30000,
  tlags = 0:7,
  cores = 8
)

akst_cv_spt <- autoKrigeST.cv(
  formula = PM10 ~ 1,
  data = deair_rs,
  nfold = 4,
  fold_dim = "spacetime",
  cutoff = 300000,
  width = 30000,
  tlags = 0:7,
  cores = 8
)
```

The result is a data frame with one row per fold:

``` r
akst_cv_t
#>   CVFold     RMSE      MAE      BIAS
#> 1      1  ...
#> 2      2  ...
#> 3      3  ...
```

The columns are:

- `CVFold`: fold identifier.
- `RMSE`: root mean squared prediction error.
- `MAE`: mean absolute prediction error.
- `BIAS`: mean prediction minus mean observed value.

Summarise the cross-validation results to compare model settings:

``` r
summary_cv <- data.frame(
  fold_dim = "temporal",
  mean_RMSE = mean(akst_cv_t$RMSE, na.rm = TRUE),
  mean_MAE = mean(akst_cv_t$MAE, na.rm = TRUE),
  mean_BIAS = mean(akst_cv_t$BIAS, na.rm = TRUE)
)

summary_cv
```

Lower `RMSE` and `MAE` indicate smaller prediction errors. `BIAS` values close
to zero indicate little systematic over- or under-prediction.

For `fold_dim = "spacetime"`, `nfold` must be a perfect square, such as `4`,
`9`, or `16`, because the folds are built from spatial groups crossed with
temporal groups.

## Notes

- Projected coordinate reference systems are recommended because variogram
  distances are computed in map units.
- Universal Kriging formulas such as `PM10 ~ elevation + traffic` require
  `new_data` with matching covariate columns.
- `predict_chunk` can be used for large prediction grids to reduce memory
  pressure.
- `variogram_from_full = TRUE` in `autoKrigeST.cv()` fits one variogram on the
  full data set and reuses it across folds, which can be much faster when that
  assumption is appropriate.
