# autoSTK
Description: Automatic spatio-temporal kriging inspired by `automap` (Hiemstra et al. 2010)

# Installation
`devtools::install_github('sigmafelix/autoSTK/autoSTK')` in R (The `devtools` package is required)

# Supported functions
- `autofitVariogramST`: automatically fit the spatiotemporal variogram
- `autoKrigeST`: automatically estimate the spatiotemporal variables for specified period (the default value is 6 temporal units)
- Currently support `STFDF` class and ordinary spatiotemporal Kriging (e.g., `y~1`)

# Update plan
- Add support for `STIDF` class and universal spatiotemporal Kriging
- Add support for the spatiotemporal separability test (referring to [`covatest`](https://cran.r-project.org/web/packages/covatest/index.html) (De Iaco, 2020))

# __Main features__
+ Split data into spatial and temporal dimensions which are compatible to be fitted as components of a spatio-temporal variogram
+ Find the optimal theoretical variograms with BFHS algorithm following the `autofitVariogram` function of `automap` package; but several changes were applied

