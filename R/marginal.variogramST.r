#' Compute the marginal spatial or temporal sample variogram
#'
#' @param stv A StVariogram and data.frame object.
#' @param bound numeric. The maximum distance that will be used to compute a spatial variogram.
#' @param spatial boolean. if TRUE, the spatial marginal variogram will be obtained, temporal otherwise.
#' @return A gstatVariogram object.
#' @export

### marginal.variogramST.R
## input: stv. StVariogram and data.frame.
## bound: numeric. maximum distance to function fit
## Subset empirical variogram by lag distance
marginal.variogramST <- function(stv, bound, spatial = TRUE) {
  # if (!grepl('^timelag.*', colnames(stv))){
  #    stop('You may want to check the input stv is valid.')
  # }
  # print(stv)
  if (spatial) {
    vg <- stv[(1:nrow(stv)) * ((stv$timelag == min(stv$timelag)) * (stv$spacelag <= bound) * !is.na(stv$gamma)), ]
    vg$np <- as.numeric(vg$np)
  } else {
    vg <- stv[(1:nrow(stv)) * ((stv$spacelag == min(stv$spacelag)) * !is.na(stv$gamma)), ]
    vg$dist <- vg$timelag
    vg$id <- 0
    vg$np <- as.numeric(vg$np)
  }
  class(vg) <- c("gstatVariogram", "data.frame")
  return(vg)
}
