# autoSTK
Description: Automatic spatio-temporal kriging inspired by `automap` (Hiemstra et al. 2010)
     
__Main features__
+ Split data into spatial and temporal dimensions which are compatible to be fitted as components of a spatio-temporal variogram
+ Find the optimal theoretical variograms with BFHS algorithm following the `autofitVariogram` function of `automap` package; but several changes were applied

