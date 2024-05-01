#' Create new spatial data for interpolation
#'
#' @param obj A Spatial*DataFrame.
#' @param gen_mode character. One of 'rect' (rectangular) and 'chull' (convex hull).
#' @param npoints integer. the number of points that will be generated
#' @param return_class character(1). One of 'sf' and 'sp'
#' @return A sf (return_class == 'sf') or SpatialPointsDataFrame (return_class == 'sp').
#' @export
create_new_data <- function(obj, gen_mode = "chull", npoints = 1e4, return_class = "sf") {
  # Function that creates a new_data object if one is missing
  distinguishsf <- function(infeat) {
    if (any(!"sf" %in% class(infeat)) | !grepl("^Spatial", class(infeat))) {
      stop("The class is neither sf nor Spatial*. Please check your input.")
    }
    if (any("sf" %in% class(infeat))) {
      return(infeat)
    } else {
      st_as_sf(infeat)
      # as(infeat, "Spatial")
    }
  }

  if (gen_mode == "rect") {
    obj.b <- st_as_sf(obj)
    genextent <- st_as_sfc(st_bbox(obj.b))
    # d <- as(obj.rect, 'Spatial')
  } else {
    genextent <- st_convex_hull(obj)
    #   convex_hull = chull(coordinates(obj)[,1],coordinates(obj)[,2])
    #   convex_hull = c(convex_hull, convex_hull[1]) # Close the polygon
    #   d = Polygon(coordinates(obj)[convex_hull, ])
    # d = convex_hull
  }
  new_data <- st_sample(genextent, npoints, type = "regular")
  # gridded(new_data) = TRUE
  if (return_class == "sp") {
    new_data <- as(new_data, "Spatial")
  }
  return(new_data)
}




#' Autodetect the temporal unit
#'
#' @param temporal xts object.
#' @return A character that indicates the temporal unit of the input xts object.
#' @export
detect_temporal_unit <- function(temporal) {
  if (!any(grepl("(zoo|xts)", class(temporal)))) {
    detect_temporal_unit.generic(temporal)
  } else {
    detect_temporal_unit.xts(temporal)
  }
}

#' Autodetect the temporal unit in a xts object
#'
#' @param temporal a xts object.
#' @return A character that indicates the temporal unit of the input xts object.
#' @export
detect_temporal_unit.xts <- function(temporal) {
  gcdt <- unique(diff(sort(unique(attributes(temporal)$index))))
  if (length(gcdt) > 1) stop("The data has the incompatible irregular temporal difference")
  if (gcdt == 86400) {
    tunit <- "days"
  } else if (gcdt == 3600) {
    tunit <- "hours"
  } else if (gcdt == 60) {
    tunit <- "minutes"
  } else {
    tunit <- "unknown"
  }
  return(tunit)
}


#' Autodetect the temporal unit in a lubridate/POSIXct object
#'
#' @param temporal a lubridate/POSIXct object.
#' @return A character that indicates the temporal unit of the input xts object.
#' @export
detect_temporal_unit.generic <- function(temporal) {
  cldt <- class(temporal)
  gcdt <- diff(sort(unique(temporal)))
  unique_gcdt <- unique(gcdt)
  len_gcdt <- length(unique_gcdt)

  if (is.null(gcdt)) {
    stop("No time difference is detected.\n")
  } else {
    if (len_gcdt > 1) {
      stop("The data has the incompatible irregular temporal difference.\n")
    } else if (len_gcdt == 0) {
      stop("No temporal difference is detected.\n")
    }
    if (any(grepl("Date", cldt))) {
      tunit <- "days"
    } else {
      tunit <- attr(gcdt, "units")
    }
  }
  return(tunit)
}



#' Generate a new spatiotemporal points for the spatiotemporal prediction and interpolation (legacy sp)
#'
#' @param obj a ST*DF object.
#' @param form formula.
#' @param gen_mode character. One of 'rect' (rectangular) and 'chull' (convex hull).
#' @param npoints integer. the number of points that will be generated
#' @param forward integer. the time length of the data will generate ahead of the last time point of the input data. If NULL is passed, the spatiotemporal interpolation mode in obj will be conducted.
#' @return A STFDF object.
#' @export
create_new_data.ST.legacy <- function(obj, form, gen_mode = "chull", npoints = 1e4, forward = 6) {
  sp.base <- obj@sp
  tunit <- detect_temporal_unit(obj@time)
  if (is.null(tunit)) {
    ts.base <- sort(unique(index(obj@time)))
    ts.base <- ts.base[length(ts.base)]
  } else {
    if (tunit == "days") {
      ts.base <- sort(unique(index(obj@time)))
      ts.base <- as.Date(ts.base)[length(ts.base)]
    } else if (tunit != "days" & tunit != "unknown") {
      ts.base <- sort(unique(index(obj@time)))
      ts.base <- ts.base[length(ts.base)]
    }
  }

  sp.base <- create_new_data(sp.base, gen_mode = gen_mode, npoints = npoints)
  sp.base <- sp::spTransform(sp.base, proj4string(obj@sp))
  if (!is.null(forward)) {
    # TODO: temporal unit-dependent setting
    if (tunit == "hours") {
      ts.base <- tunit + seq(0, 3600 * forward, 3600)
    } else if (tunit == "minutes") {
      ts.base <- ts.base + seq(0, 60 * forward, 60)
    } else {
      ts.base <- ts.base + 1:forward
    }
  } else {
    ts.base <- as.Date(sort(unique(index(obj@time))))
  }
  dat <- data.frame(dat = rep(1, length(sp.base) * length(ts.base)))
  colnames(dat)[1] <- as.character(form)[2]
  new_data_ST <- STFDF(
    sp = sp.base,
    time = ts.base,
    data = dat
  )
  new_data_ST@sp@proj4string <- obj@sp@proj4string
  return(new_data_ST)
}


#' Generate a new spatiotemporal points for the spatiotemporal prediction and interpolation
#'
#' @param obj a sftime object.
#' @param form formula.
#' @param gen_mode character. The type of shape by the surface is generated. One of 'rect' (rectangular) and 'chull' (convex hull).
#' @param npoints integer. the number of points that will be generated
#' @param forward integer. the time length of the data will generate ahead of the last time point of the input data. If NULL is passed, the spatiotemporal interpolation mode in obj will be conducted.
#' @return A sftime object.
#' @export
create_new_data.ST <- function(obj,
                               form,
                               gen_mode = "chull",
                               npoints = 1e4,
                               forward = 6) {
  st.geom <- obj[[attr(obj, "sf_column")]]
  st.time <- obj[[attr(obj, "time_column")]]
  sp.base <- st.geom
  # sp.base <- st_as_sfc(unique(st.geom))
  tunit <- detect_temporal_unit(st.time)
  if (tunit == "days") {
    ts.base <- sort(unique(st.time))
    ts.base <- as.Date(ts.base)[length(ts.base)]
  } else if (tunit != "days" & tunit != "unknown") {
    ts.base <- sort(unique(st.time))
    ts.base <- ts.base[length(ts.base)]
  }
  sp.base <- create_new_data(sp.base, gen_mode = gen_mode, npoints = npoints, return_class = "sf")
  # st_crs(sp.base) = st_crs(obj)
  # sp.base <- sp::spTransform(sp.base, proj4string(obj@sp))
  if (!is.null(forward)) {
    # TODO: temporal unit-dependent setting
    ts.base <-
      switch(tunit,
        days = ts.base + seq(1, forward, 1),
        hours = ts.base + seq(3600, 3600 * forward, 3600),
        minutes = ts.base + seq(60, 60 * forward, 60),
        seconds = ts.base + seq(1, forward, 1)
      )

    #   if (tunit == 'hours') {
    #     ts.base <- tunit + seq(0, 3600 * forward, 3600)
    #   } else if (tunit == 'minutes') {
    #     ts.base <- ts.base + seq(0, 60 * forward, 60)
    #   } else {
    #     ts.base <- ts.base + 1:forward
    #   }
  } else {
    ts.base <- sort(unique(st.time))
  }
  # dat = data.frame(dat = rep(1, length(sp.base) * length(ts.base)))
  # new_data_ST <- STFDF(sp = sp.base,
  #                      time = ts.base,
  #                      data = dat)
  NT <- length(unique(ts.base))
  NS <- length(sp.base)
  sp.base <- rep(sp.base, NT)
  ts.base <- rep(ts.base, each = NS)
  new_data_ST <- data.frame(
    outcome = NA,
    time = ts.base,
    geometry = sp.base
  )
  colnames(new_data_ST)[1] <- as.character(form)[2]
  new_data_ST <- st_as_sftime(new_data_ST,
    sf_column_name = "geometry",
    time_column_name = "time"
  )
  st_crs(new_data_ST) <- st_crs(obj)
  # new_data_ST@sp@proj4string = obj@sp@proj4string
  return(new_data_ST)
}
