# `autoSTK`

- Description: Automatic spatio-temporal kriging inspired by `automap`
  (Hiemstra et al. 2010)
- Inherits most of `automap` functions, but extensively revised
- Includes spatiotemporal variogram fitting and Kriging analysis

# Installation

- `remotes::install_github('sigmafelix/autoSTK')` in R (The `remotes`
  package is required)
- Before you install the package, please make it sure you have `gstat`,
  `spacetime`, and `automap` packages in your machine

# Supported functions

- `autofitVariogramST`: automatically fit the spatiotemporal variogram
- `autoKrigeST`: automatically estimate the spatiotemporal variables for
  specified period ahead (the default value is 6 temporal units)
  - Currently support `STFDF` and `STSDF` class and ordinary
    spatiotemporal Kriging (e.g., `y~1`)
- `autoKrigeST.cv`: automatically cross-validate the spatiotemporal data
  by spatial, temporal, spatiotemporal, and random slicing

# Update plan

- Add support for `STIDF` class and universal spatiotemporal Kriging:
  necessitates the `gstat` function fix for `variogramST.STIDF`
- Transition to `sf` and `sftime` in preparation of the retirement of
  `sp` in October 2023

# **Main features**

- Split data into spatial and temporal dimensions which are compatible
  to be fitted as components of a spatio-temporal variogram
- Find the optimal theoretical variograms with BFHS algorithm following
  the `autofitVariogram` function of `automap` package; but several
  changes were applied

# Please note

- I strongly recommend users to convert `STIDF` to `STSDF` before
  running functions
  - When the input is `sftime`, convert the input into `STIDF` then into
    `STSDF`. It is done automatically in most cases.

# Working example

- Please consult the help page of the main functions `autoKrigeST` and
  `autoKrigeST.cv`.

``` r
library(autoSTK)
library(gstat)
library(spacetime)
library(stars)
library(sp)

data(air)
deair = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
deair_sf = st_as_stars(deair) %>%
    st_transform('+proj=longlat +ellps=sphere')
deair_sf = st_transform(deair_sf, 3857)
deair_r = as(deair_sf, 'STFDF')
deair_r@sp@proj4string = CRS('+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs +type=crs')

deair_rs = deair_r[,3701:3800]
deair_rss = as(deair_rs, 'STSDF')

## autoKrigeST
akst_stk = autoKrigeST(formula = PM10~1, 
                       input_data = deair_rss, 
                       cutoff = 300000, width = 30000,
                       type_stv = "sumMetric",
                       model = c("Exc", "Mat", "Ste", "Exp", "Wav"),
                       tlags = 0:7, cores = 8,
                       verbose = TRUE)

akst_stk_stars = st_as_stars(akst_stk[[1]])
plot(akst_stk_stars[1,])
```

![](https://i.imgur.com/xDadzlJ.png)

``` r

## autoKrigeST.cv
akst_cv_t = autoKrigeST.cv(formula = PM10~1, data = deair_rs,  nfold = 3, fold_dim = 'temporal', 
                         cutoff = 300000, width = 30000, tlags = 0:7, cores = 8)
akst_cv_s = autoKrigeST.cv(formula = PM10~1, data = deair_rs,  nfold = 3, fold_dim = 'spatial', 
                         cutoff = 300000, width = 30000, tlags = 0:7, cores = 8)
akst_cv_spt = autoKrigeST.cv(formula = PM10~1, data = deair_rs,  nfold = 4, fold_dim = 'spacetime', 
                         cutoff = 300000, width = 30000, tlags = 0:7, cores = 8)
```
