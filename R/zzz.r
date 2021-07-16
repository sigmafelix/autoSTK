## Misc for npsp integration: subject to retire

svariso <- function(input, vars, maxlag = 30000, nlags = 10, estimator = 'modulus'){
    if (sum(grepl('Spatial.*', class(input))) != 0){
        .svariso.sp(input, vars, maxlag, nlags, estimator)
    } else if (sum(grepl('sf.*', class(input))) != 0){
        .svariso.sf(input, vars, maxlag, nlags, estimator)
    } else {
        stop('The input is not sf or sp-compatible dataset')
    }
}

.as.svariso.variogram <- function(vg, estimator = 'modulus'){
  sv <- list(biny = vg$gamma,
       binw = vg$np,
       grid = npsp::grid.par(n = nrow(vg),
                         min = min(vg$dist, na.rm = T),
                         max = max(vg$dist, na.rm = T),
                         lag = (vg$spacelag[2] - vg$spacelag[1])),
       data = list(x = 2, y = 1, med = weighted.mean(vg$gamma, vg$np)),
       svar = list(type = 'isotropic',
                   estimator = estimator))
  attr(sv, 'class') <- c('svar.bin', 'bin.data', 'bin.den', 'data.grid')
  return(sv)
}

.svariso.sp <- function(Sp, vars, ml, nlags, estimator = 'modulus'){
  coord.n <- coordinates(Sp)

  y <- Sp@data[,vars]
  ssp <- svar.bin(coord.n, y, estimator=estimator, maxlag=ml, nlags = nlags)
  return(ssp)
}

# sf implementation
.svariso.sf <- function(sf, vars, ml, nlags, estimator = 'modulus'){
  coord.n <- st_coordinates(sf)
  y <- st_set_geometry(sf, NULL)[,vars] %>% unlist
  ssp <- svar.bin(coord.n, y, estimator=estimator, maxlag = ml, nlags = nlags)
  return(ssp)
}
