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

create_new_data.ST <- function(obj, form = NULL, gen_mode = 'chull', npoints = 1e4, forward=NULL){
	sp.base <- obj@sp
	ts.base <- sort(as.Date(as.POSIXlt(obj@time)))
	sp.base <- create_new_data(sp.base, gen_mode = gen_mode, npoints = npoints)
	sp.base <- spTransform(sp.base, proj4string(obj@sp))
	if (!is.null(forward)){
	  # TODO: temporal unit-dependent setting
	  ts.base <- as.Date(as.POSIXlt(obj@time))[length(obj@time)] + 1:forward
	}
	dat = data.frame(dat = rep(NA, length(sp.base) * length(ts.base)))
	colnames(dat) <- as.character(form)[2]
	new_data_ST <- STFDF(sp = sp.base,
	                     time = ts.base,
	                     data = dat)
	return(new_data_ST)
}
