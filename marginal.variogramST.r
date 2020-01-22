### marginal.variogramST.R
## input: stv. StVariogram and data.frame.
## bound: numeric. maximum distance to function fit
## Subset empirical variogram by lag distance
marginal.variogramST <- function(stv, bound, spatial=TRUE){
  if (!grepl('^timelag*', colnames(stv))){
      stop('You may want to check the input stv is valid.')
  }
  if (spatial) {
    vg = stv %>% dplyr::filter(timelag == 0 & spacelag <= bound & !is.na(gamma))
  }
  else {
    vg = stv %>% dplyr::filter(spacelag == 0 & !is.na(gamma)) %>% mutate(dist = timelag, id = 0)
  }
  class(vg) = c('gstatVariogram', 'data.frame')
  return(vg)
}