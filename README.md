# `autoSTK`
- Description: Automatic spatio-temporal kriging inspired by `automap` (Hiemstra et al. 2010)
- Inherits most of `automap` functions, but extensively revised
- Includes spatiotemporal variogram fitting and Kriging analysis

# Installation
`remotes::install_github('sigmafelix/autoSTK')` in R (The `remotes` package is required)

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
