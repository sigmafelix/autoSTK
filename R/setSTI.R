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
#' @param na.omit Boolean. Omit NA values.
#' @param wireframe Boolean. Whether you plot a StVariogram in wireframe or not. If not, the return will be in class of data.frame, not a list
#' @param plot3d Boolean. Wheter you make a three-dimensional graph with rgl package
#' @param cores Integer. Number of threads passed to \code{gstat::variogramST}.
#' @return Depends on the arguments wireframe (if TRUE, list of length 2) and plot3d (if TRUE, list of length 3), a StVariogram object otherwise.
#' @export
setSTI <-
  function(stf,
           formula,
           tlags = 0:6,
           cutoff = 30000,
           width = 1000,
           assumeRegular = FALSE,
           pseudo = TRUE,
           # logarithm = FALSE,
           na.omit = TRUE,
           wireframe = FALSE,
           plot3d = FALSE,
           cores = 1) {
    formula <- as.formula(formula)

    # Keep temporal lags valid for the input series to avoid empty lag bins.
    max_tlag <- max(length(stf@time) - 1L, 0L)
    tlags <- unique(as.integer(tlags))
    tlags <- sort(tlags[is.finite(tlags) & tlags >= 0L & tlags <= max_tlag])
    if (length(tlags) == 0L) {
      tlags <- 0L
    }

    # if (logarithm){
    #   stf@data <- log(stf@data)
    # }
    # else {
    #   stf <- stf
    # }

    compute_stv <- function(cutoff_val, width_val) {
      variogramST(
        formula = formula,
        data = stf,
        tlags = tlags,
        assumeRegular = assumeRegular,
        pseudo = pseudo,
        na.omit = na.omit,
        cutoff = cutoff_val,
        width = width_val,
        cores = cores
      )
    }

    apo.pmsub.stf <- tryCatch(
      compute_stv(cutoff, width),
      error = function(e) {
        msg <- conditionMessage(e)
        if (!grepl("dim\\(X\\) must have a positive length", msg, fixed = FALSE)) {
          stop(e)
        }

        sp_coords <- sp::coordinates(stf@sp)
        if (nrow(sp_coords) < 2L) {
          stop("setSTI requires at least two spatial locations.")
        }

        d <- as.numeric(stats::dist(sp_coords))
        d <- d[is.finite(d) & d > 0]
        if (length(d) == 0L) {
          stop("setSTI could not compute positive pairwise spatial distances.")
        }

        cutoff_retry <- max(cutoff, stats::quantile(d, probs = 0.9, names = FALSE))
        width_retry <- width
        if (!is.finite(width_retry) || width_retry <= 0 || width_retry >= cutoff_retry) {
          width_retry <- cutoff_retry / 15
        }

        warning(
          "No valid ST lag bins for the given tlags/cutoff/width; retrying with adaptive spatial lags.",
          call. = FALSE
        )
        compute_stv(cutoff_retry, width_retry)
      }
    )

    ncol.stf <- length(unique(apo.pmsub.stf$spacelag))
    nrow.stf <- length(unique(apo.pmsub.stf$timelag))
    if (!wireframe & !plot3d) {
      plot.sti.set <- apo.pmsub.stf
    }

    if (wireframe) {
      label.tlag <- units(apo.pmsub.stf$timelag)
      apo.pmsub.stf$timelag <- as.numeric(apo.pmsub.stf$timelag)
      wireframe.stf <- lattice::wireframe(gamma ~ spacelag * timelag,
        apo.pmsub.stf,
        drape = TRUE,
        ylab = paste("Time lag (", label.tlag, ")", sep = ""),
        col.regions = colorRampPalette(colors = c("white", "red"))(100),
        zlim = c(0, max(apo.pmsub.stf$gamma) * 1.02)
      )
      dev.new()
      print(wireframe.stf)
      plot.sti.set <- list(apo.pmsub.stf, wireframe.stf)
    }
    if (wireframe * plot3d == 1) {
      apo.pmsub.stf.mat <- matrix(apo.pmsub.stf$gamma,
        byrow = FALSE,
        nrow = nrow.stf, ncol = ncol.stf
      )
      persp.stf <- rgl::persp3d(
        x = unique(apo.pmsub.stf$spacelag),
        y = unique(apo.pmsub.stf$timelag),
        z = apo.pmsub.stf.mat,
        color = "green3"
      )
      persp.stf
      plot.sti.set <- list(apo.pmsub.stf, wireframe.stf, persp.stf)
    }
    return(plot.sti.set)
  }
