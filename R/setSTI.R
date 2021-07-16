### setSTI.R
### AUTHOR: Insang Song (sigmafelix@hotmail.com)
### INPUT
### stf: ST*DF
### formula: formula (inherits the same parameter in variogramST)
### tlags: temporal lags to compute semivariance (inherits the same parameter in variogramST)
### cutoff: maximum bound of the set of spatial lags (inherits the same parameter in variogramST)
### width: spatial lag (inherits the same parameter in variogramST)
### logarithm: LOGICAL. log-transformation
### wireframe: LOGICAL. Whether you plot a StVariogram in wireframe or not. If not, the return will be in class of data.frame, not a list
### plot3d: LOGICAL. Wheter you make a three-dimensional graph with rgl package
#' A convenience function for the sample spatiotemporal variogram 
#'
#' @param stf A ST*DF object.
#' @param formula formula (inherits the same parameter in variogramST)
#' @param tlags temporal lags to compute semivariance (inherits the same parameter in variogramST)
#' @param cutoff the maximum bound of the set of spatial lags (inherits the same parameter in variogramST)
#' @param width integer (1). spatial lag (inherits the same parameter in variogramST)
#' @param assumeRegular Boolean. Assuming regular grid?
#' @param pseudo Boolean. See ?gstat::variogramST
#' @param logarithm Boolean. log-transformation
#' @param na.omit Boolean. Omit NA values.
#' @param wireframe Boolean. Whether you plot a StVariogram in wireframe or not. If not, the return will be in class of data.frame, not a list
#' @param plot3d Boolean. Wheter you make a three-dimensional graph with rgl package
#' @return Depends on the arguments wireframe (if TRUE, list of length 2) and plot3d (if TRUE, list of length 3), a StVariogram object otherwise.
#' @export
setSTI <-
  function(stf, 
           formula, 
           tlags = 0:6, 
           cutoff = 30000, 
           width = 1000,
           assumeRegular = TRUE, 
           pseudo = TRUE, 
           logarithm = FALSE, 
           na.omit = TRUE,
           wireframe = TRUE, 
           plot3d = FALSE, 
           cores = 1) {
    formula <- as.formula(formula)
    ncol.stf <- (cutoff / width) + 1
    nrow.stf <- max(tlags)

    if (logarithm){
      stf@data <- log(stf@data)
    }
    else {
      stf <- stf
    }

    apo.pmsub.stf <- variogramST(formula = formula,
                                 data = stf,
                                 tlags = tlags,
                                 assumeRegular = assumeRegular,
                                 pseudo = pseudo,
                                 na.omit = na.omit,
                                 cutoff = cutoff,
                                 width = width,
                                 cores = cores)
    if (!wireframe & !plot3d){
      plot.sti.set <- apo.pmsub.stf
    }

    if (wireframe) {
      wireframe.stf <- lattice::wireframe(gamma ~ spacelag * timelag,
                                          apo.pmsub.stf,
                                          drape=TRUE,
                                          col.regions = colorRampPalette(colors = c('white', 'red'))(100),
                                          zlim=c(0, max(apo.pmsub.stf$gamma)*1.02))
      dev.new()
      print(wireframe.stf)
      plot.sti.set <- list(apo.pmsub.stf, wireframe.stf)
    }
    if (wireframe * plot3d == 1) {
      apo.pmsub.stf.mat <- matrix(apo.pmsub.stf$gamma,
                                  byrow=FALSE,
                                  nrow=nrow.stf, ncol=ncol.stf)
      persp.stf <- rgl::persp3d(x=unique(apo.pmsub.stf$spacelag),
                           y=unique(apo.pmsub.stf$timelag),
                           z=apo.pmsub.stf.mat,
                           color='green3')
      persp.stf
      plot.sti.set <- list(apo.pmsub.stf, wireframe.stf, persp.stf)
    }
    return(plot.sti.set)
  }
