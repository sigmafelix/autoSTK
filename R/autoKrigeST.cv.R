## autoKrigeST.cv

#' Cross-validation of spatiotemporal Kriging
#'
#' @param data a `STFDF`-class object
#' @param fold_dim character. the dimension at which you want to cross-validate (spatial, temporal, and random)
#' @param nfold integer. the number of folds. 10 as the default
#' @param formula formula. e.g., y~1
#' @param type_stv character. One of 'sumMetric', 'metric', 'productSum', and 'separable'
#' @param block numeric. passed to conduct block spatiotemporal Kriging.
#' @param model character vector. Default is c("Sph", "Exp", "Gau", "Ste"), but users can specify the list of theoretical variograms by referring gstat::vgm.
#' @param kappa numeric vector. Kappa values tested for Matern-family variogram models.
#' @param fix.values numeric vector. Initial values in order of nugget, range, and sill, respectively.
#' @param newdata_mode character. One of 'rect' (rectangular grid) and 'chull' (convex hull)
#' @param newdata_npoints integer. The number of points that will be generated in the range of geometry the user specified (one of rectangle or convex hull)
#' @param GLS.model a variogram model. The default value is NA. If a variogram model is passed, a Generalized Lease Squares sample variogram will be calculated.
#' @param tlags integer vector (increasing, preferably to be consecutive). temporal lags.
#' @param cutoff numeric. The maximum distance at which the sample variogram will be computed.
#' @param width numeric. The interval at which the variogram cloud will be summarized.
#' @param predict_chunk integer. The number of data points per chunk in the new data for the large data. It should be meticulously chosen according to the user's machine specification.
#' @param nmax integer or positive infinite. The maximum number of spatiotemporal neighbors to conduct the local spatiotemporal Kriging.
#' @param aniso_method character. One of 'vgm', 'linear', 'range', and 'metric'. Please refer to ?gstat::estiStAni.
#' @param type_joint character. The model form of joint spatiotemporal variogram.
#' @param prodsum_k numeric. The parameter for the case when 'productSum' is chosen for type_stv.
#' @param start_vals numeric vector (3). The initial values to optimize the spatiotemporal variogram model.
#' @param miscFitOptions list. See ?automap::autofitVariogram.
#' @param cores integer. The number of threads that will be used to compute the sample spatiotemporal variogram.
#' @return The cross-validated spatiotemporal Kriging results.
#' @examples
#' library(sp)
#' library(gstat)
#' library(spacetime)
#' library(stars)
#' library(dplyr)
#' data(air)
#' deair <- STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
#' deair_sf <- st_as_stars(deair, crs = "+proj=longlat +ellps=sphere")
#' deair_sf <- st_transform(deair_sf, 3857)
#' deair_r <- as(deair_sf, "STFDF")
#' deair_r@sp@proj4string <- CRS("+init=epsg:3857")
#' deair_rs <- deair_r[, 3751:3800]
#' ## autoKrigeST.cv test
#' akst_cv_t <- autoKrigeST.cv(
#'   formula = PM10 ~ 1, data = deair_rs, nfold = 3, fold_dim = "temporal",
#'   cutoff = 300000, width = 30000, tlags = 0:7, cores = 8
#' )
#' akst_cv_s <- autoKrigeST.cv(
#'   formula = PM10 ~ 1, data = deair_rs, nfold = 3, fold_dim = "spatial",
#'   cutoff = 300000, width = 30000, tlags = 0:7, cores = 8
#' )
#' # akst_cv_r = autoKrigeST.cv(formula = PM10~1, data = deair_rs,  nfold = 3, fold_dim = 'random',
#' #                          cutoff = 300000, width = 30000, tlags = 0:7, cores = 8)
#' akst_cv_spt <- autoKrigeST.cv(
#'   formula = PM10 ~ 1, data = deair_rs, nfold = 4, fold_dim = "spacetime",
#'   cutoff = 300000, width = 30000, tlags = 0:7, cores = 8
#' )
#' @export
autoKrigeST.cv <- function(data, fold_dim,
                           nfold = 10L,
                           formula,
                           type_stv = "sumMetric",
                           block = 0,
                           model = c("Sph", "Exp", "Gau", "Ste"),
                           kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
                           fix.values = c(NA, NA, NA),
                           tlags = 0:6,
                           cutoff = 2e4,
                           width = 5e2,
                           nmax = Inf,
                           aniso_method = "vgm",
                           type_joint = "Exp",
                           prodsum_k = 0.25,
                           surface = FALSE,
                           start_vals = c(NA, NA, NA),
                           miscFitOptions = list(),
                           measurement_error = c(0, 0, 0),
                           cores = 1,
                           seed = 130425L) {
  if (!grepl("^(spa|temp|tim|rand|)*", fold_dim)) {
    stop("The argument fold_dim is not valid. Enter one of spatial, temporal, and random.")
  }

  len_time <- length(data@time)
  len_space <- length(data@sp)
  set.seed(seed)
  get_data_fold <- function(data, dimension = fold_dim, nfold) {
    data_fold <- vector("list", length = nfold)
    data_validation <- vector("list", length = nfold)

    if (grepl("^spatial$|^space$", dimension)) {
      sp_coords <- coordinates(data@sp)
      sp_coords_km <- kmeans(sp_coords, nfold)
      indices <- sp_coords_km$cluster

      vv <- split(1:len_space, indices)
      for (i in seq_len(length(vv))) {
        idata <- data[-vv[[i]], ]
        vdata <- data[vv[[i]], ]
        data_fold[[i]] <- as(idata, "STSDF")
        data_validation[[i]] <- as(vdata, "STSDF")
      }
    }
    if (grepl("^(temp|tim)", dimension)) {
      q <- len_time %% nfold
      if (q != 0) {
        targ <- floor(len_time / nfold)
        v <- rep(targ, nfold)
        vindex <- sample(nfold, q, replace = FALSE)
        v[vindex] <- v[vindex] + 1
      } else {
        targ <- len_time / nfold
        v <- rep(targ, nfold)
      }
      vv <- split(c(1:len_time), rep(1:nfold, v))

      for (i in seq_len(length(vv))) {
        idata <- data[, -vv[[i]]]
        vdata <- data[, vv[[i]]]
        data_fold[[i]] <- as(idata, "STSDF")
        data_validation[[i]] <- as(vdata, "STSDF")
      }
    }
    if (grepl("^(spatiotemp|spacetime)", dimension)) {
      if (sqrt(nfold) %% 1 > 0) {
        stop("Spatiotemporal CV is only available for integer sqrt(nfold).")
      }
      q_t <- len_time %% sqrt(nfold)
      if (q_t != 0) {
        targ <- ceiling(len_time / nfold)
        v_t <- rep(targ, sqrt(nfold))
        v_tindex <- sample(sqrt(nfold), q_t, replace = FALSE)
        v_t[v_tindex] <- v[v_tindex] + 1
      } else {
        targ <- len_time / sqrt(nfold)
        v_t <- rep(targ, sqrt(nfold))
      }

      sp_coords <- coordinates(data@sp)
      sp_coords_km <- kmeans(sp_coords, sqrt(nfold))
      indices <- sp_coords_km$cluster

      vv_sp <- split(1:len_space, indices)
      vv_t <- split(c(1:len_time), rep(1:sqrt(nfold), v_t))
      # vv_sp = split(c(1:len_space), rep(1:sqrt(nfold), v_sp))

      for (i in seq_len(length(vv_sp))) {
        for (j in seq_len(length(vv_t))) {
          idata <- data[-vv_sp[[i]], -vv_t[[j]]]
          vdata <- data[vv_sp[[i]], vv_t[[j]]]
          data_fold[[sqrt(nfold) * (i - 1) + j]] <- as(idata, "STSDF")
          data_validation[[sqrt(nfold) * (i - 1) + j]] <- as(vdata, "STSDF")
        }
      }
    }

    if (grepl("^rand", dimension)) {
      for (i in 1:nfold) {
        data_sti <- as(data, "STIDF")
        indices <- sample(nrow(data_sti@data), ceiling(nrow(data_sti@data) / nfold), replace = FALSE)
        indices <- sort(indices, decreasing = TRUE)
        data_fold[[i]] <- data_sti[-indices, ]
        data_validation[[i]] <- data_sti[indices, ]
      }
    }
    return(list(data_fold, data_validation))
  }

  folded <- get_data_fold(data = data, dimension = fold_dim, nfold = nfold)
  data_fold <- folded[[1]]
  data_validation <- folded[[2]]

  cvresult <- mapply(function(id, vd) {
    kriged <- autoKrigeST(
      formula = formula,
      input_data = id,
      new_data = vd,
      type_stv = type_stv,
      model = model,
      block = block,
      kappa = kappa,
      fix.values = fix.values,
      tlags = tlags,
      cutoff = cutoff,
      width = width,
      nmax = nmax,
      prodsum_k = prodsum_k,
      start_vals = start_vals,
      miscFitOptions = miscFitOptions,
      measurement_error = measurement_error,
      cores = cores
    )$krige_output$var1.pred
    form_char <- as.character(formula)
    actual <- as.vector(vd@data[, form_char[2]])

    RMSE <- sqrt(mean((kriged - actual)^2))
    MAE <- mean(abs(kriged - actual))
    BIAS <- mean(kriged) - mean(actual)
    return(list(RMSE = RMSE, MAE = MAE, BIAS = BIAS))
  }, data_fold, data_validation, SIMPLIFY = FALSE)

  res <- do.call(rbind, cvresult)
  res <- data.frame(cbind(CVFold = 1:nfold, res))
  return(res)
}
