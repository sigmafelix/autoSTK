### plot.autofitVariogram.R
### Original Author: Paul Hiemstra (paul@numbertheory.nl)
### Partial fix: Insang Song (sigmafelix@hotmail.com)
### Fix
#### Easy implementation of the main title
#### Change of the parameter direction to the numeric vector of length 3

plot.autofitVariogram <- 
function (x, plotit = TRUE, title="Experimental variogram and fitted variogram model", ...) 
{
  shift = 0.03
  labels = as.character(x$exp_var$np)
  vario = xyplot(gamma ~ dist, data = x$exp_var, panel = automap:::autokrige.vgm.panel, 
                 labels = labels, shift = shift, model = x$var_model, 
                 direction = c(1, 0, 0), 
                 ylim = c(min(0, 1.04 * min(x$exp_var$gamma)), 1.04 * 
                            max(x$exp_var$gamma)), xlim = c(0, 1.04 * max(x$exp_var$dist)), 
                 xlab = "Distance", ylab = "Semi-variance", main = title, 
                 mode = "direct", ...)
  if (plotit) 
    print(vario)
  else vario
}
