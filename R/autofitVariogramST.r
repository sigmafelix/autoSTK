### autofitVariogramST.R
### Author: Insang Song (sigmafelix@hotmail.com)
### INPUT
#### stf: ST*DF
#### formula: formula
#### typestv: character. type of spatio-temporal covariance function.
####          should be one of c('separable','sumMetric','simplesumMetric','productSum','productSumOld','metric')
#### guess_nugget: spatio-temporal joint nugget
#### guess_psill: spatio-temporal joint partial sill
#### aniso_method: character. method of estimating spatio-temporal anisotropy ratio.
####               should be one of c('range','linear','vgm')
#### type_joint: model type of spatio-temporal joint process. only used when typestv is metric-variant type
#### prodsum_k: positive numeric. only used for productSum or productSumOld
#### surface: whether plot the empirical spatiotemporal variogram or not
#### measurement_error: measurement error variance component. Refer to ?gstat::vgm; the input vector includes spatial, temporal, and joint measurement error components.
#### cores: passed to variogramST. The number of cores to be used to compute semivariogram values
#' @title Automatic Spatio-Temporal Variogram Fitting
#' @description
#' Fits a spatio-temporal variogram model to a given spatio-temporal data frame using automatic parameter estimation.
#'
#' @param stf A spatio-temporal data frame (STFDF or similar) containing the data to be modeled.
#' @param formula A formula specifying the response and explanatory variables.
#' @param typestv Character string specifying the type of spatio-temporal variogram model to fit. Options include "sumMetric", "separable", "productSum", "productSumOld", "simpleSumMetric", and "metric".
#' @param candidate_model Character vector of candidate spatial/temporal variogram models to consider (e.g., c("Ste", "Exc", "Exp", "Wav")).
#' @param guess_nugget Optional initial guess for the nugget parameter. If NULL, it is estimated automatically.
#' @param guess_psill Optional initial guess for the partial sill parameter. If NULL, it is estimated automatically.
#' @param tlags Integer vector specifying the temporal lags to use in the empirical variogram calculation.
#' @param cutoff Numeric value specifying the maximum spatial lag distance to consider.
#' @param width Numeric value specifying the width of spatial lag bins.
#' @param aniso_method Character string specifying the method for estimating spatio-temporal anisotropy ("vgm" or "linear").
#' @param type_joint Character string specifying the model type for the joint variogram component (e.g., "Exp").
#' @param prodsum_k Optional numeric value for the product-sum model parameter k. If NULL, it is estimated automatically.
#' @param surface Logical; if TRUE, computes and returns the variogram surface.
#' @param measurement_error Numeric vector of length 3 specifying measurement error for spatial, temporal, and joint variograms, respectively.
#' @param cores Integer specifying the number of CPU cores to use for parallel computation.
#' @param verbose Logical; if TRUE, prints progress and diagnostic messages.
#'
#' @return A list containing:
#'   \item{jointSTV}{The fitted spatio-temporal variogram model.}
#'   \item{empSTV}{The empirical spatio-temporal variogram.}
#'   \item{SpV}{The fitted spatial variogram model.}
#'   \item{TV}{The fitted temporal variogram model.}
#'   \item{STVsurface}{(Optional) The computed variogram surface, if \code{surface = TRUE}.}
#'
#' @details
#' This function automates the process of fitting a spatio-temporal variogram by:
#' \itemize{
#'   \item Computing empirical spatial and temporal variograms.
#'   \item Automatically selecting and fitting marginal models.
#'   \item Estimating spatio-temporal anisotropy.
#'   \item Constructing and fitting a joint spatio-temporal variogram model.
#' }
#' The function supports several popular spatio-temporal variogram structures, including separable, product-sum, and sum-metric models.
#'
#' @seealso \code{\link{vgmST}}, \code{\link{fit.StVariogram}}, \code{\link{marginal.variogramST}}
#' @examples
#' library(spacetime)
#' library(sp)
#' library(stars)
#' library(sftime)
#' data(air)
#' stations <- st_as_sf(stations)
#' stations <- st_transform(stations, 'EPSG:3857')
#'
#' airdf <- data.frame(PM10 = as.vector(air))
#' stations_full <- do.call(c, rep(stations, length(dates)))
#' dates_full <- rep(dates, each = nrow(stations))
#'
#' rural <- cbind(airdf, time = dates_full, stations_full)
#'
#' rural <- st_as_sftime(rural, sf_column_name = 'geometry')
#' rr <- rural[match(rural$time, dates[3001:3060], nomatch = FALSE, incomparables = FALSE) > 0, ]
#' rr <- as(as(as(rr, "STIDF"), "STFDF"), "STSDF")
#' rrstv <- autofitVariogramST(
#'  stf = rr, formula = PM10 ~ 1, surface = TRUE)
#' rrstv
#' @export
autofitVariogramST <- function(stf,
                               formula,
                               typestv = "sumMetric",
                               candidate_model = c("Ste", "Exc", "Exp", "Wav"),
                               guess_nugget = NULL,
                               guess_psill = NULL,
                               tlags = 0:6,
                               cutoff = 2e4,
                               width = 5e2,
                               aniso_method = "vgm",
                               type_joint = "Exp",
                               prodsum_k = NULL,
                               surface = FALSE,
                               measurement_error = c(0, 0, 0),
                               cores = 1,
                               verbose = FALSE) {
  stva <- setSTI(
    stf = stf,
    formula = formula,
    tlags = tlags,
    cutoff = cutoff,
    width = width,
    wireframe = FALSE,
    cores = cores
  )
  if (is.null(candidate_model)) {
    stva.sp <- marginal.variogramST(stva,
      bound = cutoff
    )
    stva.ts <- marginal.variogramST(stva,
      spatial = FALSE
    )
    stva.sp.fit <- autofitVariogram_(
      formula = NULL, verbose = verbose,
      input_data = NULL,
      input_vgm = stva.sp,
      measurement_error = measurement_error[1],
      model = NULL,
    )
    stva.ts.fit <- autofitVariogram_(
      formula = NULL, verbose = verbose,
      input_data = NULL,
      input_vgm = stva.ts,
      measurement_error = measurement_error[2],
      model = NULL
    )
    aniso_method <- "linear"
  } else {
    stva.sp <- marginal.variogramST(stva,
      bound = cutoff
    )
    stva.ts <- marginal.variogramST(stva,
      spatial = FALSE
    )
    stva.sp$np <- as.numeric(stva.sp$np)
    stva.ts$np <- as.numeric(stva.ts$np)

    stva.sp.fit <- autofitVariogram_(
      formula = NULL, verbose = verbose,
      input_data = NULL,
      input_vgm = stva.sp,
      measurement_error = measurement_error[1],
      model = candidate_model
    )
    stva.ts.fit <- autofitVariogram_(
      formula = NULL, verbose = verbose,
      input_data = NULL,
      input_vgm = stva.ts,
      measurement_error = measurement_error[2],
      model = candidate_model
    )
  }

  # if (typestv == 'separable'){
  #  sill.init <- median(stva$gamma)
  #  vgm.spatial$psill <- vgm.spatial$psill / sill.init
  #  vgm.temporal$psill <- vgm.temporal$psill / sill.init
  # }

  stv.ani <- estiStAni(stva,
    interval = c(0.2, 2) * median(stva$spacelag),
    method = aniso_method,
    spatialVgm = stva.sp.fit$var_model,
    temporalVgm = stva.ts.fit$var_model
  )
  if (is.null(guess_nugget)) {
    guess_nugget <- max(min(stva$gamma), min(stva$gamma) - 0.5 * (min(stva.sp$gamma) + min(stva.ts$gamma)))
    guess_nugget <- ifelse(guess_nugget < 0, 0, guess_nugget)
  }
  if (is.null(guess_psill)) {
    guess_psill_c1 <- ifelse(0.5 * (max(stva$gamma) - max(stva.sp$gamma, stva.ts$gamma)) < 0,
      min(stva$gamma) * 2, 0.5 * (max(stva$gamma) - max(stva.sp$gamma, stva.ts$gamma))
    )
    guess_psill_c2 <- ifelse(0.5 * (stva$gamma[length(stva$gamma)] - max(stva.sp$gamma, stva.ts$gamma)) < 0,
      min(stva$gamma) * 2, 0.5 * (stva$gamma[length(stva$gamma)] - max(stva.sp$gamma, stva.ts$gamma))
    )
    if (typestv == "metric") {
      guess_psill <- 0.5 * max(stva.sp$gamma)
    } else {
      guess_psill <- max(0.1 * max(stva.sp$gamma), min(guess_psill_c1, guess_psill_c2))
    }
  }
  sill <- max(stva$gamma) * 0.6
  stv.jo <- vgm(
    model = type_joint,
    psill = guess_psill, nugget = guess_nugget,
    Err = measurement_error[3],
    range = 0.25 * sqrt((stv.ani)^2 + (stv.ani * max(stva$spacelag)^2))
  )

  if (is.null(prodsum_k)) {
    prodsum_k <- 4 / max(stva.sp$gamma)
  }

  ## sp phi, sigmasq, tausq - ts-joint - total sill and nugget,  stani
  variost.mod <- switch(typestv,
    separable = vgmST(
      stModel = typestv, space = stva.sp.fit$var_model,
      time = stva.ts.fit$var_model, sill = sill,
      nugget = guess_nugget
    ),
    productSum = vgmST(
      stModel = typestv,
      space = stva.sp.fit$var_model, time = stva.ts.fit$var_model,
      k = prodsum_k
    ),
    productSumOld = vgmST(
      stModel = typestv,
      space = stva.sp.fit$var_model, time = stva.ts.fit$var_model,
      sill = guess_psill * sqrt(2), nugget = guess_nugget
    ),
    sumMetric = vgmST(
      stModel = typestv, space = stva.sp.fit$var_model, time = stva.ts.fit$var_model,
      joint = stv.jo, stAni = stv.ani
    ),
    simpleSumMetric = vgmST(
      stModel = typestv,
      space = stva.sp.fit$var_model, time = stva.ts.fit$var_model,
      joint = stv.jo,
      nugget = guess_nugget, stAni = stv.ani
    ),
    metric = vgmST(stModel = typestv, joint = stv.jo, stAni = stv.ani),
    stop(paste("model", typest, "unknown"))
  )

  joint.lower <- extractPar(variost.mod) * 0.25
  joint.upper <- extractPar(variost.mod) * 1.5
  if (any(is.infinite(joint.lower))) joint.lower[which(is.infinite(joint.lower))] <- 0
  if (any(is.infinite(joint.upper))) joint.upper[which(is.infinite(joint.upper))] <- 1e4
  if (any(is.na(joint.lower))) joint.lower[which(is.na(joint.lower))] <- 0
  if (any(is.na(joint.upper))) joint.upper[which(is.na(joint.upper))] <- 1e4
  if (any(joint.lower < 0)) joint.lower[which(joint.lower < 0)] <- 0
  if (any(joint.upper < 0)) joint.upper[which(joint.upper < 0)] <- 1e2

  stva.joint <- fit.StVariogram(
    object = stva, model = variost.mod, stAni = stv.ani,
    method = "L-BFGS-B",
    lower = joint.lower, upper = joint.upper,
    control = list(maxit = 2.5e3, REPORT = 1)
  )

  if (surface) {
    STVS <- variogramSurface(stva.joint, stva[, c("timelag", "spacelag")])
    autofitSTV <- list(
      jointSTV = stva.joint,
      empSTV = stva,
      SpV = stva.sp.fit,
      TV = stva.ts.fit,
      STVsurface = STVS
    )
  } else {
    autofitSTV <- list(
      jointSTV = stva.joint,
      empSTV = stva,
      SpV = stva.sp.fit,
      TV = stva.ts.fit
    )
  }

  return(autofitSTV)
}
