### autofitVariogramST.R
### Author: Insang Song (sigmafelix@hotmail.com)

#' @title Automatic Spatio-Temporal Variogram Fitting
#' @description
#' Fits a spatio-temporal variogram model to a spatio-temporal data frame using
#' automatic parameter estimation.  Supports multiple optimisers and objective
#' functions (v2.0+).
#'
#' @param stf A spatio-temporal data frame (STFDF, STSDF, STIDF, or sftime).
#' @param formula A formula specifying the response and explanatory variables.
#' @param typestv Character. ST variogram model type. One of \code{"sumMetric"},
#'   \code{"separable"}, \code{"productSum"}, \code{"productSumOld"},
#'   \code{"simpleSumMetric"}, \code{"metric"}.
#' @param candidate_model Character vector of candidate spatial/temporal
#'   variogram model families (e.g. \code{c("Ste", "Exc", "Exp", "Wav")}).
#' @param guess_nugget Optional initial guess for the nugget. Estimated if NULL.
#' @param guess_psill Optional initial guess for the partial sill. Estimated if NULL.
#' @param tlags Integer vector of temporal lags.
#' @param cutoff Numeric. Maximum spatial lag distance.
#' @param width Numeric. Width of spatial lag bins.
#' @param aniso_method Character. Method for estimating ST anisotropy ratio
#'   (\code{"vgm"} or \code{"linear"}).
#' @param type_joint Character. Model family for the joint variogram component.
#' @param prodsum_k Numeric. k parameter for productSum models.
#' @param surface Logical. If TRUE, also compute and return the variogram surface.
#' @param measurement_error Numeric vector of length 3: spatial, temporal, joint.
#' @param cores Integer. Number of cores for \code{variogramST}.
#' @param verbose Logical. If TRUE, print diagnostic messages.
#' @param optimizer Character. Optimisation strategy: \code{"lbfgsb"} (default,
#'   current behaviour) or \code{"grid"} (LHS grid search + L-BFGS-B).
#' @param objective Character. Fitting criterion: \code{"WLS"} (weighted least
#'   squares, default) or \code{"MLE"} (Gaussian log-likelihood; only for
#'   n < 500).
#' @param n_restart Integer. Number of optimisation restarts for
#'   \code{optimizer = "lbfgsb"}.  Default 1 matches v1.x behaviour.
#' @param optimizer_control Named list of extra arguments forwarded to the
#'   optimiser (e.g. \code{list(n_coarse = 100L)} for \code{"grid"}).
#'
#' @return A \code{STVariogramFit} object (list) with elements:
#'   \describe{
#'     \item{jointSTV}{Fitted ST variogram model.}
#'     \item{empSTV}{Empirical ST variogram.}
#'     \item{SpV}{Fitted spatial marginal variogram.}
#'     \item{TV}{Fitted temporal marginal variogram.}
#'     \item{STVsurface}{(Optional) Variogram surface if \code{surface = TRUE}.}
#'     \item{optimizer}{Character. Optimiser used.}
#'     \item{objective}{Character. Objective function used.}
#'     \item{loglik}{Numeric. Log-likelihood (only when \code{objective = "MLE"}).}
#'     \item{n_obs}{Integer. Number of observations (for AIC/BIC).}
#'   }
#'
#' @seealso \code{\link{fitVariogramST}}, \code{\link{selectModelST}},
#'   \code{\link{test_separability}}
#' @examples
#' library(spacetime)
#' library(sp)
#' library(sftime)
#' data(air)
#' stations <- st_as_sf(stations)
#' stations <- st_transform(stations, 'EPSG:3857')
#' airdf <- data.frame(PM10 = as.vector(air))
#' stations_full <- do.call(c, rep(stations, length(dates)))
#' dates_full <- rep(dates, each = nrow(stations))
#' rural <- cbind(airdf, time = dates_full, stations_full)
#' rural <- st_as_sftime(rural, sf_column_name = 'geometry')
#' rr <- rural[match(rural$time, dates[3001:3060], nomatch = FALSE) > 0, ]
#' rr <- as(as(as(rr, "STIDF"), "STFDF"), "STSDF")
#' rrstv <- autofitVariogramST(stf = rr, formula = PM10 ~ 1, surface = TRUE)
#' rrstv
#' @export
autofitVariogramST <- function(stf,
                               formula,
                               typestv           = "sumMetric",
                               candidate_model   = c("Ste", "Exc", "Exp", "Wav"),
                               guess_nugget      = NULL,
                               guess_psill       = NULL,
                               tlags             = 0:6,
                               cutoff            = 2e4,
                               width             = 5e2,
                               aniso_method      = "vgm",
                               type_joint        = "Exp",
                               prodsum_k         = NULL,
                               surface           = FALSE,
                               measurement_error = c(0, 0, 0),
                               cores             = 1L,
                               verbose           = FALSE,
                               optimizer         = c("lbfgsb", "grid"),
                               objective         = c("WLS", "MLE"),
                               n_restart         = 1L,
                               optimizer_control = list()) {
  optimizer <- match.arg(optimizer)
  objective <- match.arg(objective)

  # ---- Coerce input to STSDF -----------------------------------------------
  stf <- .to_stsdf(stf)

  # ---- Empirical ST variogram ----------------------------------------------
  stva <- setSTI(
    stf       = stf,
    formula   = formula,
    tlags     = tlags,
    cutoff    = cutoff,
    width     = width,
    wireframe = FALSE,
    cores     = cores
  )

  # ---- Marginal variograms -------------------------------------------------
  if (is.null(candidate_model)) {
    stva.sp <- marginal.variogramST(stva, bound = cutoff)
    stva.ts <- marginal.variogramST(stva, spatial = FALSE)
    stva.sp.fit <- autofitVariogram_(
      formula = NULL, verbose = verbose,
      input_data = NULL, input_vgm = stva.sp,
      measurement_error = measurement_error[1], model = NULL
    )
    stva.ts.fit <- autofitVariogram_(
      formula = NULL, verbose = verbose,
      input_data = NULL, input_vgm = stva.ts,
      measurement_error = measurement_error[2], model = NULL
    )
    aniso_method <- "linear"
  } else {
    stva.sp <- marginal.variogramST(stva, bound = cutoff)
    stva.ts <- marginal.variogramST(stva, spatial = FALSE)
    stva.sp$np <- as.numeric(stva.sp$np)
    stva.ts$np <- as.numeric(stva.ts$np)
    stva.sp.fit <- autofitVariogram_(
      formula = NULL, verbose = verbose,
      input_data = NULL, input_vgm = stva.sp,
      measurement_error = measurement_error[1], model = candidate_model
    )
    stva.ts.fit <- autofitVariogram_(
      formula = NULL, verbose = verbose,
      input_data = NULL, input_vgm = stva.ts,
      measurement_error = measurement_error[2], model = candidate_model
    )
  }

  # ---- Anisotropy ratio ---------------------------------------------------
  stv.ani <- estiStAni(
    stva,
    interval    = c(0.2, 2) * median(stva$spacelag),
    method      = aniso_method,
    spatialVgm  = stva.sp.fit$var_model,
    temporalVgm = stva.ts.fit$var_model
  )

  # ---- Initial parameter guesses ------------------------------------------
  if (is.null(guess_nugget)) {
    guess_nugget <- max(
      min(stva$gamma),
      min(stva$gamma) - 0.5 * (min(stva.sp$gamma) + min(stva.ts$gamma))
    )
    guess_nugget <- ifelse(guess_nugget < 0, 0, guess_nugget)
  }
  if (is.null(guess_psill)) {
    guess_psill_c1 <- ifelse(
      0.5 * (max(stva$gamma) - max(stva.sp$gamma, stva.ts$gamma)) < 0,
      min(stva$gamma) * 2,
      0.5 * (max(stva$gamma) - max(stva.sp$gamma, stva.ts$gamma))
    )
    guess_psill_c2 <- ifelse(
      0.5 * (stva$gamma[length(stva$gamma)] -
               max(stva.sp$gamma, stva.ts$gamma)) < 0,
      min(stva$gamma) * 2,
      0.5 * (stva$gamma[length(stva$gamma)] -
               max(stva.sp$gamma, stva.ts$gamma))
    )
    if (typestv == "metric") {
      guess_psill <- 0.5 * max(stva.sp$gamma)
    } else {
      guess_psill <- max(0.1 * max(stva.sp$gamma),
                         min(guess_psill_c1, guess_psill_c2))
    }
  }
  if (is.null(prodsum_k)) {
    prodsum_k <- 4 / max(stva.sp$gamma)
  }

  sill    <- max(stva$gamma) * 0.6
  stv.jo  <- vgm(
    model  = type_joint,
    psill  = guess_psill, nugget = guess_nugget,
    Err    = measurement_error[3],
    range  = 0.25 * sqrt(stv.ani^2 + (stv.ani * max(stva$spacelag))^2)
  )

  # ---- Initial vgmST model ------------------------------------------------
  variost.mod <- switch(typestv,
    separable = vgmST(
      stModel = typestv, space = stva.sp.fit$var_model,
      time = stva.ts.fit$var_model, sill = sill, nugget = guess_nugget
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
      stModel = typestv, space = stva.sp.fit$var_model,
      time = stva.ts.fit$var_model, joint = stv.jo, stAni = stv.ani
    ),
    simpleSumMetric = vgmST(
      stModel = typestv,
      space = stva.sp.fit$var_model, time = stva.ts.fit$var_model,
      joint = stv.jo, nugget = guess_nugget, stAni = stv.ani
    ),
    metric = vgmST(stModel = typestv, joint = stv.jo, stAni = stv.ani),
    stop(paste("model", typestv, "unknown"))
  )

  # ---- Bounds (improved semantics via st_param_bounds) --------------------
  bounds <- st_param_bounds(variost.mod, stva)

  # Apply legacy bound safety checks on top of st_param_bounds
  bounds$lower[is.infinite(bounds$lower) | is.na(bounds$lower)] <- 0
  bounds$upper[is.infinite(bounds$upper) | is.na(bounds$upper)] <- 1e4
  bounds$lower[bounds$lower < 0] <- 0
  bounds$upper[bounds$upper < 0] <- 1e2

  # ---- Optimise -----------------------------------------------------------
  ctrl <- modifyList(
    list(n_coarse = 50L, n_refine = 30L, maxit = 2500L),
    optimizer_control
  )

  stva.joint <- switch(optimizer,
    lbfgsb = .optimize_stv_lbfgsb(
      stva_emp       = stva,
      model_template = variost.mod,
      bounds         = bounds,
      n_restart      = n_restart,
      maxit          = ctrl$maxit
    ),
    grid = .optimize_stv_grid(
      stva_emp       = stva,
      model_template = variost.mod,
      bounds         = bounds,
      n_coarse       = ctrl$n_coarse,
      n_refine       = ctrl$n_refine,
      maxit          = ctrl$maxit
    )
  )

  # ---- Optional: compute log-likelihood -----------------------------------
  loglik <- NULL
  n_obs  <- nrow(stf@data)

  if (objective == "MLE") {
    response_col <- as.character(formula)[2L]
    loglik <- .loglik_stv(stva.joint, stf, response_col)
    if (verbose)
      message("MLE log-likelihood: ", signif(loglik, 6L))
  }

  # ---- Assemble output ----------------------------------------------------
  if (surface) {
    STVS <- variogramSurface(stva.joint, stva[, c("timelag", "spacelag")])
  }

  result <- list(
    jointSTV  = stva.joint,
    empSTV    = stva,
    SpV       = stva.sp.fit,
    TV        = stva.ts.fit,
    optimizer = optimizer,
    objective = objective,
    loglik    = loglik,
    n_obs     = n_obs
  )
  if (surface) result$STVsurface <- STVS

  class(result) <- c("STVariogramFit", "list")
  result
}
