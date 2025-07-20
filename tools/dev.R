library(sp)
library(spacetime)
sumMetricVgm <- vgmST("sumMetric",
                      space = vgm( 4.4, "Lin", 196.6,  3),
                      time  = vgm( 2.2, "Lin",   1.1,  2),
                      joint = vgm(34.6, "Exp", 136.6, 12),
                      stAni = 51.7)
data(air)
suppressWarnings(proj4string(stations) <- CRS(proj4string(stations)))
rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
rr <- rural[,"2005-06-01/2005-06-03"]
rr <- as(rr,"STSDF")
x1 <- seq(from=6,to=15,by=1)
x2 <- seq(from=48,to=55,by=1)
DE_gridded <- SpatialPoints(cbind(rep(x1,length(x2)), rep(x2,each=length(x1))), 
                            proj4string=CRS(proj4string(rr@sp)))
gridded(DE_gridded) <- TRUE
DE_pred <- STF(sp=as(DE_gridded,"SpatialPoints"), time=rr@time)

rr@sp <- sp::spTransform(rr@sp, CRS("EPSG:3857"))
air_stk0 <- autoKrigeST(formula = PM10~1,
                       input_data = rr,
                       type_stv = 'sumMetric',
                       tlags = 0:7,
                       model = c('Mat', 'Ste', 'Wav', 'Exp', 'Exc'),
                       cutoff = 3e5,
                       width = 3e4,
                       cores = 4)


data(air)
deair <- STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
deair_sf <- st_as_stars(deair) %>%
  st_transform("+proj=longlat +ellps=sphere")
deair_sf <- st_transform(deair_sf, 3857)
deair_r <- as(deair_sf, "STFDF")
deair_rs <- deair_r[, 3751:3800]

## autoKrigeST.cv test
akst = autoKrigeST(formula = PM10~1, input_data = deair_rs,
                         cutoff = 300000, width = 30000, tlags = 0:7, cores = 8)



DE_kriged <- krigeST(PM10~1, data=rr, newdata=DE_pred,
                     modelList=sumMetricVgm)
gridded(DE_kriged@sp) <- TRUE
stplot(DE_kriged)





###
library(spacetime)
library(sp)
library(stars)
library(sftime)
data(air)
stations <- st_as_sf(stations)
stations <- st_transform(stations, 'EPSG:3857')

airdf <- data.frame(PM10 = as.vector(air))
stations_full <- do.call(c, rep(stations, length(dates)))
dates_full <- rep(dates, each = nrow(stations))

rural <- cbind(airdf, time = dates_full, stations_full)

#ruralsft = st_as_stars(airdf, time = dates)
rural <- st_as_sftime(rural, sf_column_name = 'geometry')
rr <- rural[match(rural$time, dates[3001:3060], nomatch = FALSE, incomparables = FALSE) > 0, ]
# rr[is.na(rr$PM10), 'PM10'] <- mean(rr$PM10, na.rm = TRUE)

air_stk <- autoKrigeST(formula = PM10~1,
                       input_data = rr,
                       type_stv = 'metric',
                       tlags = 0:7,
                       model = c('Mat', 'Ste', 'Wav', 'Exp', 'Exc'),
                       cutoff = 3e5,
                       width = 3e4,
                       cores = 4)
