create_new_data = function(obj, mode = NULL, npoints = 1e4) {
# Function that creates a new_data object if one is missing
	convex_hull = chull(coordinates(obj)[,1],coordinates(obj)[,2])
	convex_hull = c(convex_hull, convex_hull[1]) # Close the polygon
	d = Polygon(coordinates(obj)[convex_hull, ]) 

	if (mode == 'rect'){
		obj.b <- st_as_sf(obj)
		obj.rect <- st_as_sfc(st_bbox(obj.b))
		d <- as(obj.rect, 'Spatial')
	}
	new_data = spsample(d, npoints, type = "regular")
	gridded(new_data) = TRUE
	return(new_data)
}

create_new_data.ST <- function(obj, mode = NULL){
	sp.base <- obj@sp
	ts.base <- unique(obj@time)
	sp.base <- create_new_data(sp.base, mode = mode)
	dat = data.frame(dat = rep(NA, length(sp.base) * length(ts.base)))
	new_data_ST <- STFDF(sp = sp.base,
						 time = ts.base,
						 data = dat)
	return(new_data_ST)
}