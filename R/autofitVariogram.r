autofitVariogram_ <- function(formula, input_data, input_vgm = NULL, model = c("Sph", "Exp", "Gau", "Ste"),
                             kappa = c(0.05, seq(0.1, 2, 0.1), 5, 10), fix.values = c(NA, NA, NA),
                             verbose = FALSE, GLS.model = NA, start_vals = c(NA, NA, NA), measurement_error = 0,
                             boundaries = c(2, 4, 6, 9, 12, 15, 25, 35, 50, 65, 80, 100) * 50000,
                             miscFitOptions = list(), ...)
# This function automatically fits a variogram to input_data
{
  getModel <- function(psill, model, range, kappa, nugget, fit_range, fit_sill, fit_nugget, measurement_error = 0, verbose) {
    if (verbose) debug.level <- 1 else debug.level <- 0
    if (model == "Pow") {
      warning("Using the power model is at your own risk, read the docs of autofitVariogram for more details.")
      if (is.na(start_vals[1])) nugget <- 0
      if (is.na(start_vals[2])) range <- 1 # If a power mode, range == 1 is a better start value
      if (is.na(start_vals[3])) sill <- 1
    }
    vgm_try <- vgm(
      psill = psill, model = model, range = range,
      nugget = nugget, kappa = kappa, Err = measurement_error
    )
    obj <- try(
      fit.variogram(experimental_variogram,
        model = vgm_try,
        fit.ranges = c(fit_range),
        fit.sills = c(fit_nugget, fit_sill),
        debug.level = 0
      ),
      TRUE
    )
    if ("try-error" %in% class(obj)) {
      # print(traceback())
      warning(
        "An error has occured during variogram fitting. Used:\n",
        "\tnugget:\t", nugget,
        "\n\tmodel:\t", model,
        "\n\tpsill:\t", psill,
        "\n\trange:\t", range,
        "\n\tkappa:\t", ifelse(kappa == 0, NA, kappa),
        "\n  as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n", obj
      )
      return(NULL)
    } else {
      return(obj)
    }
  }

  # if there is an input data
  if (!is.null(input_data)) {
    # Create boundaries
    # if ST* or Spatial* object
    if (sum(grepl("(Spatial|ST).*", class(input_data))) > 0) {
      longlat <- !is.projected(input_data)
      if (is.na(longlat)) {
        longlat <- FALSE
      }
      diagonal <- spDists(t(bbox(input_data)), longlat = longlat)[2]
    } else {
      longlat <- st_is_longlat(input_data)
      diagonal <- st_as_sfc(st_bbox(input_data))
      diagonal <- st_cast(diagonal, "POINT")
      diagonal <- st_distance(diagonal[c(1, 3)], longlat = longlat)[2, 1]
      diagonal <- as.vector(diagonal)
    }
    if (is.null(boundaries)) {
      boundaries <- c(2, 4, 6, 9, 12, 15, 25, 35, 50, 65, 80, 100) * diagonal * 0.35 / 100 # # 0.35 times the length of the central axis through the area, Boundaries for the bins in km
    }
  }

  if (!is.null(input_data) & is.null(input_vgm)) {
    # Check for anisotropy parameters
    if ("alpha" %in% names(list(...))) warning("Anisotropic variogram model fitting not supported, see the documentation of autofitVariogram for more details.")

    # Take the misc fit options and overwrite the defaults by the user specified ones
    miscFitOptionsDefaults <- list(merge.small.bins = TRUE, min.np.bin = 5)
    miscFitOptions <- modifyList(miscFitOptionsDefaults, miscFitOptions)

    # If you specifiy a variogram model in GLS.model the Generelised least squares sample variogram is constructed
    if (!is(GLS.model, "variogramModel")) {
      experimental_variogram <- variogram(formula, input_data, boundaries = boundaries, ...)
    } else {
      if (verbose) cat("Calculating GLS sample variogram\n")
      g <- gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls = 1))
      experimental_variogram <- variogram(g, boundaries = boundaries, ...)
    }

    # request by Jon Skoien
    if (miscFitOptions[["merge.small.bins"]]) {
      if (verbose) cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
      while (TRUE) {
        if (length(experimental_variogram$np[experimental_variogram$np < miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) == 1) break
        boundaries <- boundaries[2:length(boundaries)]
        if (!is(GLS.model, "variogramModel")) {
          experimental_variogram <- variogram(formula, input_data, boundaries = boundaries, ...)
        } else {
          experimental_variogram <- variogram(g, boundaries = boundaries, ...)
        }
      }
    }
  } else {
    experimental_variogram <- input_vgm
    diagonal <- max(input_vgm$dist)
  }


  if (!is.null(model) | !is.null(input_vgm)) {
    # set initial values
    if (is.na(start_vals[1])) { # Nugget
      initial_nugget <- min(experimental_variogram$gamma) * 0.5
    } else {
      initial_nugget <- start_vals[1]
    }
    if (is.na(start_vals[2])) { # Range
      initial_range <- 0.1 * diagonal # 0.10 times the length of the central axis through the area
    } else {
      initial_range <- start_vals[2]
    }
    if (is.na(start_vals[3])) { # Sill
      initial_sill <- max(experimental_variogram$gamma)
    } else {
      initial_sill <- start_vals[3]
    }

    # Determine what should be automatically fitted and what should be fixed
    # Nugget
    if (!is.na(fix.values[1])) {
      fit_nugget <- FALSE
      initial_nugget <- fix.values[1]
    } else {
      fit_nugget <- TRUE
    }
    # Range
    if (!is.na(fix.values[2])) {
      fit_range <- FALSE
      initial_range <- fix.values[2]
    } else {
      fit_range <- TRUE
    }
    # Partial sill
    if (!is.na(fix.values[3])) {
      fit_sill <- FALSE
      initial_sill <- fix.values[3]
    } else {
      fit_sill <- TRUE
    }

    if (measurement_error != 0) {
      # initial_nugget = min(initial_nugget - measurement_error, measurement_error)
      initial_sill <- max(initial_sill - measurement_error, measurement_error)
    }


    # Automatically testing different models, the one with the smallest sums-of-squares is chosen
    test_models <- model
    SSerr_list <- c()
    vgm_list <- list()
    counter <- 1

    for (m in test_models) {
      if (m != "Mat" && m != "Ste") { # If not Matern and not Stein
        model_fit <- getModel(initial_sill - initial_nugget, m, initial_range, kappa = 0, initial_nugget, fit_range, fit_sill, fit_nugget, verbose = verbose, measurement_error = measurement_error)
        if (!is.null(model_fit)) { # skip models that failed
          vgm_list[[counter]] <- model_fit
          SSerr_list <- c(SSerr_list, attr(model_fit, "SSErr"))
        }
        counter <- counter + 1
      } else { # Else loop also over kappa values
        for (k in kappa) {
          model_fit <- getModel(initial_sill - initial_nugget, m, initial_range, k, initial_nugget, fit_range, fit_sill, fit_nugget, verbose = verbose, measurement_error = measurement_error)
          if (!is.null(model_fit)) {
            vgm_list[[counter]] <- model_fit
            SSerr_list <- c(SSerr_list, attr(model_fit, "SSErr"))
          }
          counter <- counter + 1
        }
      }
    }

    # Check for negative values in sill or range coming from fit.variogram
    # and NULL values in vgm_list, and remove those with a warning
    strange_entries <- sapply(vgm_list, function(v) any(c(v$psill, v$range) < 0) | is.null(v))
    if (any(strange_entries)) {
      if (verbose) {
        print(vgm_list[strange_entries])
        cat("^^^ ABOVE MODELS WERE REMOVED ^^^\n\n")
      }
      warning("Some models where removed for being either NULL or having a negative sill/range/nugget, \n\tset verbose == TRUE for more information")
      SSerr_list <- SSerr_list[!strange_entries]
      vgm_list <- vgm_list[!strange_entries]
    }

    if (verbose) {
      cat("Selected:\n")
      print(vgm_list[[which.min(SSerr_list)]])
      cat("\nTested models, best first:\n")
      tested <- data.frame(
        "Tested models" = sapply(vgm_list, function(x) as.character(x[2, 1])),
        kappa = sapply(vgm_list, function(x) as.character(x[2, 4])),
        "SSerror" = SSerr_list
      )
      tested <- tested[order(tested$SSerror), ]
      print(tested)
    }

    result <- list(exp_var = experimental_variogram, var_model = vgm_list[[which.min(SSerr_list)]], sserr = min(SSerr_list))
  }

  # if (is.null(model) & )
  if (is.null(model) & !is.null(input_data)) {
    svar.af <- svariso(input_data, as.character(formula)[2], maxlag = boundaries[length(boundaries)], nlags = length(boundaries))
    svar.af <- fitsvar.sb.iso(svar.af, dk = 10)
    result <- list(exp_var = svar.af, var_model = vgm.tab.svarmod(svar.af, seq(0, svar.af$range, length = 5000)))
  }
  if (is.null(model) & !is.null(input_vgm)) {
    svar.af <- .as.svariso.variogram(input_vgm)
    svar.af <- fitsvar.sb.iso(svar.af, dk = 10)
    result <- list(exp_var = svar.af, var_model = vgm.tab.svarmod(svar.af, seq(0, svar.af$range, length = 5000)))
  }


  class(result) <- c("autofitVariogram", "list")
  return(result)
}
