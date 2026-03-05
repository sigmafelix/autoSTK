# last revision: autoSTK v2.0.0

#' Fit a Spatio-Temporal Variogram
#'
#' A named wrapper around \code{\link{autofitVariogramST}} that makes the
#' fit/predict separation explicit.  Returns a \code{STVariogramFit} object
#' which can be passed directly to \code{\link{predictKrigeST}}.
#'
#' @param formula A formula (e.g. \code{PM10 ~ 1}).
#' @param data A spatio-temporal data object (STFDF, STSDF, STIDF, or sftime).
#' @param ... Arguments forwarded to \code{\link{autofitVariogramST}}.
#' @return A \code{STVariogramFit} object.
#' @seealso \code{\link{predictKrigeST}}, \code{\link{autoKrigeST}}
#' @export
fitVariogramST <- function(formula, data, ...) {
  autofitVariogramST(stf = data, formula = formula, ...)
}

#' Predict Using a Fitted Spatio-Temporal Kriging Model
#'
#' Takes a \code{STVariogramFit} object (from \code{\link{fitVariogramST}} or
#' \code{\link{autofitVariogramST}}) and performs spatiotemporal Kriging at the
#' locations in \code{newdata}.
#'
#' @param fit A \code{STVariogramFit} object.
#' @param data The original data used to fit the variogram (STFDF, STSDF, STIDF,
#'   or sftime).
#' @param newdata Prediction grid (STFDF or sftime).
#' @param formula Formula used during fitting.
#' @param nmax Passed to \code{gstat::krigeST}.
#' @param predict_chunk Optional integer; splits \code{newdata} into chunks for
#'   large grids to avoid memory issues.
#' @param ... Additional arguments forwarded to \code{gstat::krigeST}.
#' @return An \code{autoKrigeST} object (list) with elements \code{krige_output}
#'   (STFDF) and \code{var_model}.
#' @seealso \code{\link{fitVariogramST}}, \code{\link{autoKrigeST}}
#' @export
predictKrigeST <- function(fit, data, newdata, formula,
                            nmax          = Inf,
                            predict_chunk = NULL,
                            ...) {
  stopifnot(inherits(fit, "STVariogramFit"))

  input_data <- .to_stsdf(data)
  new_data   <- .to_stfdf(newdata)

  if (!is.null(predict_chunk)) {
    len_new_data <- dim(new_data)[1L]
    len_chunks   <- ceiling(len_new_data / predict_chunk)
    idx_start    <- (predict_chunk * seq(0, len_chunks - 1L, 1L)) + 1L
    idx_end      <- idx_start + predict_chunk - 1L
    idx_end[length(idx_end)] <- len_new_data

    pb <- txtProgressBar(style = 3L, max = len_chunks)
    krige_results_l <- vector("list", length = len_chunks)
    for (i in seq_len(len_chunks)) {
      krige_results_l[[i]] <- krigeST(
        formula   = formula,
        data      = input_data,
        newdata   = new_data[idx_start[i]:idx_end[i], ],
        nmax      = nmax,
        computeVar = TRUE,
        bufferNmax = 2L,
        modelList  = fit$jointSTV,
        ...
      )
      setTxtProgressBar(pb, i)
    }
    close(pb)
    krige_result <- do.call("rbind", krige_results_l)
  } else {
    krige_result <- krigeST(
      formula    = formula,
      data       = input_data,
      newdata    = new_data,
      nmax       = nmax,
      computeVar = TRUE,
      bufferNmax = 2L,
      modelList  = fit$jointSTV,
      ...
    )
  }

  krige_result <- as(krige_result, "STFDF")

  result <- list(
    krige_output = krige_result,
    var_model    = fit$jointSTV
  )
  class(result) <- c("autoKrigeST", "list")
  result
}

#' Automatic Spatio-Temporal Kriging
#'
#' Convenience wrapper that fits the variogram and predicts in one call.
#' Internally calls \code{\link{fitVariogramST}} then
#' \code{\link{predictKrigeST}}.
#'
#' @param formula A formula (e.g. \code{PM10 ~ 1}).
#' @param input_data Training data (STFDF, STSDF, STIDF, or sftime).
#' @param new_data Prediction grid. If missing, generated automatically via
#'   \code{\link{create_new_data.ST}}.
#' @param type_stv Character. ST variogram model type. Passed as \code{typestv}
#'   to \code{\link{autofitVariogramST}}.
#' @param block Numeric. Block size for block Kriging.
#' @param model Character vector. Candidate marginal variogram model families.
#' @param kappa Numeric vector. Kappa values for Matern-family models.
#' @param fix.values Numeric vector (3). Fixed initial values for nugget,
#'   range, sill.
#' @param newdata_mode Character. \code{"rect"} or \code{"chull"}.
#' @param newdata_npoints Integer. Number of prediction grid points.
#' @param GLS.model Variogram model for GLS sample variogram (rarely needed).
#' @param tlags Integer vector. Temporal lags.
#' @param cutoff Numeric. Maximum spatial lag.
#' @param width Numeric. Spatial lag bin width.
#' @param forward Integer. Number of future time steps to predict.
#' @param predict_chunk Integer or NULL. Chunk size for large prediction grids.
#' @param nmax Numeric. Maximum number of ST neighbours.
#' @param aniso_method Character. Anisotropy estimation method.
#' @param type_joint Character. Joint variogram component model family.
#' @param prodsum_k Numeric. k parameter for productSum models.
#' @param surface Logical. Return variogram surface?
#' @param start_vals Numeric vector (3). Not used in v2.0 (kept for API compat).
#' @param miscFitOptions List. Not used in v2.0 (kept for API compat).
#' @param measurement_error Numeric vector (3). Measurement error components.
#' @param cores Integer. Cores for variogramST computation.
#' @param verbose Logical. Print progress messages.
#' @param optimizer Character. \code{"lbfgsb"} or \code{"grid"}.
#' @param objective Character. \code{"WLS"} or \code{"MLE"}.
#' @param n_restart Integer. Number of optimisation restarts.
#' @param optimizer_control List. Extra control arguments for the optimiser.
#' @return An \code{autoKrigeST} object.
#' @examples
#' data(air)
#' deair <- STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
#' deair_rs <- deair[, 3751:3800]
#' ## Not run:
#' # akst <- autoKrigeST(formula = PM10 ~ 1, input_data = deair_rs,
#' #                     cutoff = 300000, width = 30000, tlags = 0:7, cores = 4)
#' @export
autoKrigeST <- function(formula,
                        input_data,
                        new_data,
                        type_stv          = "sumMetric",
                        block             = 0,
                        model             = c("Sph", "Exp", "Gau", "Ste"),
                        kappa             = c(0.05, seq(0.2, 2, 0.1), 5, 10),
                        fix.values        = c(NA, NA, NA),
                        newdata_mode      = "rect",
                        newdata_npoints   = 3e3,
                        GLS.model         = NA,
                        tlags             = 0:6,
                        cutoff            = 2e4,
                        width             = 5e2,
                        forward           = 6,
                        predict_chunk     = NULL,
                        nmax              = Inf,
                        aniso_method      = "vgm",
                        type_joint        = "Exp",
                        prodsum_k         = 0.25,
                        surface           = FALSE,
                        start_vals        = c(NA, NA, NA),
                        miscFitOptions    = list(),
                        measurement_error = c(0, 0, 0),
                        cores             = 1L,
                        verbose           = TRUE,
                        optimizer         = "lbfgsb",
                        objective         = "WLS",
                        n_restart         = 1L,
                        optimizer_control = list()) {

  # Universal Kriging requires new_data for predictors
  if (as.character(formula)[3L] != "1" && missing(new_data))
    stop("Universal Kriging requires new_data to supply covariate values.")

  # Identical-value guard
  col_name <- as.character(formula)[2L]
  if (length(unique(input_data[[col_name]])) == 1L)
    stop(sprintf(
      "All values of '%s' are identical (%s). Cannot interpolate.",
      col_name, unique(input_data[[col_name]])[1L]
    ))

  # Build prediction grid if not supplied
  if (missing(new_data)) {
    new_data <- create_new_data.ST(input_data,
      form     = formula,
      gen_mode = newdata_mode,
      npoints  = newdata_npoints,
      forward  = forward
    )
  }

  if (inherits(input_data, "ST")) {
    input_data <- sftime::st_as_sftime(input_data)
  }
  if (inherits(new_data, "ST")) {
    new_data <- sftime::st_as_sftime(new_data)
  }

  # CRS compatibility check
  crs_in  <- sf::st_crs(sftime::st_as_sftime(input_data))
  crs_new <- sf::st_crs(new_data)
  if (!all(is.na(c(crs_in$wkt, crs_new$wkt)))) {
    if (is.na(crs_in$wkt) && !is.na(crs_new$wkt)) {
      sf::st_crs(input_data) <- crs_new
    } else if (!is.na(crs_in$wkt) && is.na(crs_new$wkt)) {
      sf::st_crs(new_data) <- crs_in
    } else if (!identical(crs_in$wkt, crs_new$wkt)) {
      stop("Projections of input_data and new_data do not match:\n",
           "  input_data: ", crs_in$input, "\n",
           "  new_data:   ", crs_new$input)
    }
  }

  # Coerce using single unified helper (replaces the 3-step chain)
  input_stsdf <- .to_stsdf(input_data)
  new_stfdf   <- .to_stfdf(new_data)

  message("Fitting the optimal spatiotemporal variogram model...")
  variogram_object <- autofitVariogramST(
    stf               = input_stsdf,
    formula           = formula,
    typestv           = type_stv,
    candidate_model   = model,
    tlags             = tlags,
    cutoff            = cutoff,
    width             = width,
    aniso_method      = aniso_method,
    type_joint        = type_joint,
    prodsum_k         = if (prodsum_k == 0.25) NULL else prodsum_k,
    surface           = surface,
    measurement_error = measurement_error,
    cores             = cores,
    verbose           = verbose,
    optimizer         = optimizer,
    objective         = objective,
    n_restart         = n_restart,
    optimizer_control = optimizer_control
  )

  message("Predicting ", forward, " time step(s)...")
  predictKrigeST(
    fit           = variogram_object,
    data          = input_stsdf,
    newdata       = new_stfdf,
    formula       = formula,
    nmax          = nmax,
    predict_chunk = predict_chunk
  )
}
