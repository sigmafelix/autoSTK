# last revision: 08/12/2021
#' Cross-validation of spatiotemporal Kriging
#' @description Performs an automatic Kriging on the data with ST*DF or sftime object
#' @param data a `sftime`-class or `ST*DF`-class object. `ST*DF`-class object will be converted to sftime automatically.
#' @param fold_dim character. the dimension at which you want to cross-validate (spatial, temporal, and random)
#' @param nfold integer. the number of folds. 10 as the default.
#' @return The cross-validated spatiotemporal Kriging results.
#' @examples
#' data(air)
#' deair <- STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
#' deair_sf <- st_as_stars(deair) %>%
#'   st_transform("+proj=longlat +ellps=sphere")
#' deair_sf <- st_transform(deair_sf, 3857)
#' deair_r <- as(deair_sf, "STFDF")
#' deair_rs <- deair_r[, 3751:3800]

#' ## autoKrigeST.cv test
#' akst = autoKrigeST(formula = PM10~1, data = deair_rs,
#'                          cutoff = 300000, width = 30000, tlags = 0:7, cores = 8)
autoKrigeST <- function(formula,
                        input_data, new_data,
                        type_stv = "sumMetric",
                        # data_variogram = input_data,
                        block = 0,
                        model = c("Sph", "Exp", "Gau", "Ste"),
                        kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
                        fix.values = c(NA, NA, NA),
                        # remove_duplicates = TRUE,
                        newdata_mode = "rect",
                        newdata_npoints = 3e3,
                        GLS.model = NA,
                        tlags = 0:6,
                        cutoff = 2e4,
                        width = 5e2,
                        forward = 6,
                        predict_chunk = NULL,
                        nmax = Inf,
                        aniso_method = "vgm",
                        type_joint = "Exp",
                        prodsum_k = 0.25,
                        surface = FALSE,
                        start_vals = c(NA, NA, NA),
                        miscFitOptions = list(),
                        measurement_error = c(0, 0, 0),
                        cores = 1,
                        verbose = TRUE) {
  if (inherits(input_data, "STIDF") | inherits(input_data, "STFDF") | inherits(input_data, "STSDF")) {
    input_data <- sftime::st_as_sftime(input_data)
    # input_data = formula
    # formula = as.formula(paste(names(input_data)[1], "~ 1"))
  }

  # if((!inherits(input_data,"STFDF") & !inherits(data_variogram,"STFDF") & !inherits(input_data,"STIDF") & !inherits(data_variogram,"STIDF") & !inherits(input_data,"STSDF") & !inherits(data_variogram,"STSDF")))
  # {
  #     stop(paste("\nInvalid input objects: input_data or data_variogram not of class 'ST*DF'.\n\tClass input_data: '",
  #                   class(input_data),"'",
  #                   "\n\tClass data_variogram: '",
  #                   class(data_variogram),"'",sep=''))
  # }

  if (as.character(formula)[3] != 1 & missing(new_data)) stop("If you want to use Universal Kriging, new_data needs to be specified \n  because the predictors are also required on the prediction of spatio-temporal data elements.")
  # TODO: remove duplicate objects in ST*DF (spatially)
  # Check if there are points or gridcells on the exact same coordinate and provide a more informative error message.
  # Points on the same spot causes the interpolation to crash.
  # if(remove_duplicates)
  # {
  #    zd = zerodist(input_data)
  #    if(length(zd) != 0)
  #    {
  #        warning("Removed ", length(zd) / 2, " duplicate observation(s) in input_data:", immediate. = TRUE)
  # 	print(input_data[c(zd), ])
  #      	input_data = input_data[-zd[, 2], ]

  #    }
  # }

  # If all the values return an informative error
  col_name <- as.character(formula)[2]
  if (length(unique(input_data[[col_name]])) == 1) {
    stop(sprintf("All data in attribute \'%s\' is identical and equal to %s\n   Can not interpolate this data", col_name, unique(input_data[[col_name]])[1]))
  }

  if (missing(new_data)) {
    new_data <- create_new_data.ST(input_data,
      form = formula,
      gen_mode = newdata_mode,
      npoints = newdata_npoints,
      forward = forward
    )
  }

  ## Perform some checks on the projection systems of input_data and new_data
  crs_in <- sf::st_crs(input_data)
  crs_new <- sf::st_crs(new_data)
  if (!all(is.na(c(crs_in, crs_new)))) {
    if (is.na(crs_in) & !is.na(crs_new)) {
      sf::st_crs(input_data) <- sf::st_crs(new_data)
    } else if (!is.na(crs_in) & is.na(crs_new)) {
      sf::st_crs(new_data) <- sf::st_crs(input_data)
    } else if (crs_in[["wkt"]] != crs_new[["wkt"]]) {
      stop(paste(
        "Projections of input_data and new_data do not match:\n",
        "  input_data: ", crs_in, "\n",
        "  new_data:    ", crs_new, "\n"
      ))
    } else if (any(!grepl("PROJ", crs_in[["wkt"]]), !grepl("PROJ", crs_new[["wkt"]]))) {
      stop(paste(
        "One of input_data and new_data is in LongLat, please reproject.\n",
        "  input_data: ", crs_in, "\n",
        "  new_data:   ", crs_new, "\n"
      ))
    }
  } else {
    stop("One of input_data and new_data has no coordinate systems.\n")
  }

  # Fit the variogram model, first check which model is used
  # will be deprecated after v.2.0
  input_data <- as(as(as(input_data, "STIDF"), "STFDF"), "STSDF")
  new_data <- as(as(new_data, "STIDF"), "STFDF")
  data_variogram <- input_data

  message("Fitting the optimal spatiotemporal variogram model...")
  variogram_object <- autofitVariogramST(
    formula = formula,
    stf = data_variogram,
    typestv = type_stv,
    candidate_model = model,
    tlags = tlags,
    cutoff = cutoff,
    width = width,
    aniso_method = aniso_method,
    type_joint = type_joint,
    prodsum_k = prodsum_k,
    surface = surface,
    measurement_error = measurement_error,
    cores = cores,
    verbose = verbose
  )

  message("Predicting the output for ", forward, " time steps...")
  if (!is.null(predict_chunk)) {
    ## Perform the interpolation by chunk
    # new_data_sti = as(new_data, 'STIDF')
    len_new_data <- dim(new_data)[1]
    len_chunks <- ceiling(len_new_data / predict_chunk)
    len_indices_start <- (predict_chunk * seq(0, len_chunks - 1, 1)) + 1
    len_indices_end <- len_indices_start + predict_chunk - 1
    len_indices_end[length(len_indices_end)] <- len_new_data

    krige_results_l <- vector("list", length = len_chunks)
    pb <- txtProgressBar(style = 3, max = len_chunks)

    for (i in 1:len_chunks) {
      krige_results_l[[i]] <- krigeST(
        formula = formula,
        data = input_data,
        newdata = new_data[len_indices_start[i]:len_indices_end[i], ],
        nmax = nmax,
        computeVar = TRUE,
        bufferNmax = 2,
        modelList = variogram_object$jointSTV
      )
      setTxtProgressBar(pb, i)
    }
    close(pb)

    krige_result <- do.call("rbind", krige_results_l)
  } else {
    krige_result <- krigeST(
      formula = formula,
      data = input_data,
      newdata = new_data,
      nmax = nmax,
      computeVar = TRUE,
      bufferNmax = 2,
      modelList = variogram_object$jointSTV
    )
  }

  krige_result <- as(krige_result, "STFDF")

  # Aggregate the results into an autoKrige object
  result <- list(
    krige_output = krige_result,
    var_model = variogram_object$jointSTV
  )
  class(result) <- c("autoKrigeST", "list")

  return(result)
}
