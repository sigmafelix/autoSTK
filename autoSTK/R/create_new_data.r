create_new_data = function(obj, gen_mode = 'chull', npoints = 1e4) {
# Function that creates a new_data object if one is missing

	if (gen_mode == 'rect') {
		obj.b <- st_as_sf(obj)
		obj.rect <- st_as_sfc(st_bbox(obj.b))
		d <- as(obj.rect, 'Spatial')
	} else {
	  convex_hull = chull(coordinates(obj)[,1],coordinates(obj)[,2])
	  convex_hull = c(convex_hull, convex_hull[1]) # Close the polygon
	  d = Polygon(coordinates(obj)[convex_hull, ])
	}
	new_data = spsample(d, npoints, type = "regular")
	gridded(new_data) = TRUE
	return(new_data)
}


detect_temporal_unit <- function(temporal){
  if (sum(class(temporal) %in% c('zoo', 'xts')) < 1) stop('Please check the class of the input data')

  gcdt <- unique(diff(sort(unique(attributes(temporal)$index))))
  if (length(gcdt) > 1) stop('The data has the incompatible irregular temporal difference')
  if (gcdt == 86400){
    tunit <- 'days'
  } else if (gcdt == 3600){
    tunit <- 'hours'
  } else if (gcdt == 60){
    tunit <- 'minutes'
  } else {
    tunit <- 'unknown'
  }
  return(tunit)
}


create_new_data.ST <- function(obj, form = NULL, gen_mode = 'chull', npoints = 1e4, forward=6){
	sp.base <- obj@sp
	tunit <- detect_temporal_unit(obj@time)
	if (tunit == 'days'){
	  ts.base <- sort(unique(index(obj@time)))
	  ts.base <- as.Date(ts.base)[length(ts.base)]
	} else if (tunit != 'days' & tunit != 'unknown'){
	  ts.base <- sort(unique(index(obj@time)))
	  ts.base <- ts.base[length(ts.base)]
	}
	sp.base <- create_new_data(sp.base, gen_mode = gen_mode, npoints = npoints)
	sp.base <- spTransform(sp.base, proj4string(obj@sp))
	if (!is.null(forward)){
	  # TODO: temporal unit-dependent setting
	  if (tunit == 'hours'){
	    ts.base <- tunit + seq(0, 3600 * forward, 3600)
	  } else if (tunit == 'minutes') {
	    ts.base <- ts.base + seq(0, 60 * forward, 60)
	  } else {
	    ts.base <- ts.base + 1:forward
	  }
	}
	dat = data.frame(dat = rep(NA, length(sp.base) * length(ts.base)))
	colnames(dat) <- as.character(form)[2]
	new_data_ST <- STFDF(sp = sp.base,
	                     time = ts.base,
	                     data = dat)
	return(new_data_ST)
}
