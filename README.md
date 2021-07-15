# `autoSTK`
- Description: Automatic spatio-temporal kriging inspired by `automap` (Hiemstra et al. 2010)
- Inherits most of `automap` functions, but extensively revised
- Includes spatiotemporal variogram fitting and Kriging analysis

# Installation
- `remotes::install_github('sigmafelix/autoSTK')` in R (The `remotes` package is required)
- Before you install the package, please make it sure you have `gstat`, `spacetime`, and `automap` packages in your machine

# Supported functions
- `autofitVariogramST`: automatically fit the spatiotemporal variogram
- `autoKrigeST`: automatically estimate the spatiotemporal variables for specified period ahead (the default value is 6 temporal units)
    - Currently support `STFDF` and `STSDF` class and ordinary spatiotemporal Kriging (e.g., `y~1`)
- `autoKrigeST.cv`: automatically cross-validate the spatiotemporal data by spatial, temporal, spatiotemporal, and random slicing

# Update plan
- Add support for `STIDF` class and universal spatiotemporal Kriging: necessitates the `gstat` function fix for `variogramST.STIDF`
- Add support for the spatiotemporal separability test (referring to [`covatest`](https://cran.r-project.org/web/packages/covatest/index.html) (De Iaco, 2020))

# __Main features__
+ Split data into spatial and temporal dimensions which are compatible to be fitted as components of a spatio-temporal variogram
+ Find the optimal theoretical variograms with BFHS algorithm following the `autofitVariogram` function of `automap` package; but several changes were applied

# Please note
- Some functions related to `npsp` package will retire as the `npsp` package is no longer available on CRAN.

# Working example
- Please consult the help page of the main functions `autoKrigeST` and `autoKrigeST.cv`.

``` r
library(autoSTK)
library(gstat)
#> Warning: package 'gstat' was built under R version 4.0.5
library(spacetime)
#> Warning: package 'spacetime' was built under R version 4.0.5
library(stars)
#> Warning: package 'stars' was built under R version 4.0.5
#> Loading required package: abind
#> Loading required package: sf
#> Warning: package 'sf' was built under R version 4.0.5
#> Linking to GEOS 3.9.0, GDAL 3.2.1, PROJ 7.2.1
library(sp)
#> Warning: package 'sp' was built under R version 4.0.3

data(air)
deair = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
deair_sf = st_as_stars(deair) %>%
    st_transform('+proj=longlat +ellps=sphere')
deair_sf = st_transform(deair_sf, 3857)
deair_r = as(deair_sf, 'STFDF')
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded ellps WGS 84 in Proj4 definition: +proj=merc +a=6378137
#> +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null
#> +wktext +no_defs +type=crs
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum World Geodetic System 1984 in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded ellps WGS 84 in Proj4 definition: +proj=merc +a=6378137
#> +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null
#> +wktext +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum World Geodetic System 1984 in Proj4 definition
deair_r@sp@proj4string = CRS('+init=epsg:3857')
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded ellps WGS 84 in Proj4 definition: +proj=merc +a=6378137
#> +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null
#> +wktext +no_defs
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum WGS_1984 in Proj4 definition
deair_rs = deair_r[,3701:3800]
deair_rss = as(deair_rs, 'STSDF')

## autoKrigeST
akst_stk = autoKrigeST(formula = PM10~1, 
                       input_data = deair_rss, 
                       cutoff = 300000, width = 30000, tlags = 0:7, cores = 8)
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded ellps WGS 84 in Proj4 definition: +proj=merc +a=6378137
#> +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null
#> +wktext +no_defs +type=crs
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum World Geodetic System 1984 in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded ellps WGS 84 in Proj4 definition: +proj=merc +a=6378137
#> +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null
#> +wktext +no_defs +type=crs
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

<sup>Created on 2021-07-15 by the [reprex package](https://reprex.tidyverse.org) (v2.0.0)</sup>
