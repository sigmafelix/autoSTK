## Misc for npsp integration

svariso <- function(input, vars, maxlag = 30000, nlags = 10, estimator = 'cressie'){
    if (sum(grepl('Spatial.*', class(input))) > 1){
        .svariso.sp(input, vars, maxlag, nlags, estimator)
    } else if (sum(grepl('sf.*', class(input))) > 1){
        .svariso.sf(input, vars, maxlag, nlags, estimator)
    } else {
        stop('The input is not sf or sp-compatible dataset')
    }
}

.svariso.sp <- function(Sp, vars, ml, nlags, estimator = 'cressie'){
  coord.n <- coordinates(Sp)

  y <- Sp@data[,vars]
  ssp <- svar.bin(coord.n, y, estimator=estimator, maxlag=ml, nlags = nlags)
  return(ssp)
}

# sf implementation
.svariso.sf <- function(sf, vars, ml, nlags, estimator = 'cressie'){
  coord.n <- st_coordinates(sf)
  y <- st_set_geometry(sf, NULL)[,vars] %>% unlist
  ssp <- svar.bin(coord.n, y, estimator=estimator, maxlag = ml, nlags = nlags)
  return(ssp)
}
