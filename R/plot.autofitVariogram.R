### plot.autofitVariogram.R
### Original Author: Paul Hiemstra (paul@numbertheory.nl)
### Partial fix: Insang Song (sigmafelix@hotmail.com)
### Fix
#### Easy implementation of the main title
#### Change of the parameter direction to the numeric vector of length 3


autokrige.vgm.panel <-
  function(x, y, model, subscripts, ...) {
    lattice::panel.xyplot(x, y, ...)
    model_line <- tryCatch(
      gstat::variogramLine(model, maxdist = max(x, na.rm = TRUE), n = 200L),
      error = function(e) NULL
    )
    if (!is.null(model_line)) {
      lattice::panel.lines(model_line$dist, model_line$gamma, col = "blue", lwd = 2)
    }
    no_digits <- function(a) {
      if (a > 10) {
        return(0)
      } else {
        if (a < 1) {
          return(2)
        } else {
          return(1)
        }
      }
    }
    nugget <- sum(model[1, "psill"])
    sill <- sum(model[, "psill"])
    range <- sum(model[, "range"])
    txt <- paste("Model: ", as.character(model[2, "model"]), "\nNugget: ",
      round(nugget, digits = no_digits(nugget)), "\nSill: ",
      round(sill, digits = no_digits(sill)), "\nRange: ", round(range,
        digits = no_digits(range)
      ),
      sep = ""
    )
    if (model[2, "model"] %in% c("Mat", "Ste")) {
      kappa <- model[2, "kappa"]
      txt <- paste(txt, "\nKappa: ", round(kappa, digits = no_digits(kappa)),
        sep = ""
      )
    }
    lattice::ltext(max(x), 0.02 * max(y), txt, font = 2, cex = 0.7, adj = c(
      1,
      0
    ), col = grDevices::grey(0.3))
  }


#' Plot the automatically fitted variogram
#'
#' @importFrom graphics plot
#' @method plot autofitVariogram
#' @param x A result object of autofitVariogram.
#' @param plotit boolean. Print graph or not.
#' @param title character. the title of the plot.
#' @param ... passed to xyplot
#' @return A lattice::xyplot object.
#' @export

plot.autofitVariogram <-
  function(x, plotit = TRUE, title = "Experimental variogram and fitted variogram model", ...) {
    shift <- 0.03
    labels <- as.character(x$exp_var$np)
    vario <- lattice::xyplot(gamma ~ dist,
      data = x$exp_var, panel = autokrige.vgm.panel,
      labels = labels, shift = shift, model = x$var_model,
      direction = c(1, 0, 0),
      ylim = c(min(0, 1.04 * min(x$exp_var$gamma)), 1.04 *
        max(x$exp_var$gamma)), xlim = c(0, 1.04 * max(x$exp_var$dist)),
      xlab = "Distance", ylab = "Semi-variance", main = title,
      mode = "direct", ...
    )
    if (plotit) {
      print(vario)
    } else {
      vario
    }
  }
