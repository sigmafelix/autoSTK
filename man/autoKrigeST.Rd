\name{autoKrigeST}
\alias{autoKrigeST}
\title{Performs an automatic spatiotemporal interpolation}
\description{This function performs automatic spatiotemporal kriging on the given dataset.
                The variogram is generated automatically using \link{autofitVariogramST}.}
\usage{autoKrigeST(formula,
                   input_data,
                   new_data,
                   type_stv = 'sumMetric',
                   data_variogram = input_data,
                   block = 0,
                   model = c("Sph", "Exp", "Gau", "Ste"),
                   kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
				   fix.values = c(NA, NA, NA),
                   newdata_mode = 'rect',
				   newdata_npoints = 3e3,
                   GLS.model = NA,
          		   tlags = 0:6,
          		   cutoff = 2e4,
          		   width = 5e2,
                   predict_chunk = NULL,
                   nmax = Inf,
          		   aniso_method = 'vgm',
          		   type_joint = 'Exp',
          		   prodsum_k = 0.25,
          		   start_vals = c(NA, NA, NA),
                   miscFitOptions = list(),
                   measurement_error = c(0,0,0),
                   cores = 1,
                   verbose = TRUE)
                    }
\arguments{
	\item{formula}{formula that defines the dependent variable as a linear model
				of independent variables; suppose the dependent variable has
				name 'z', for ordinary and simple kriging use the formula
				'z~1'; for simple kriging also define 'beta' (see below); for
				universal kriging, suppose 'z' is linearly dependent on 'x'
				and 'y', use the formula 'z~x+y'.}
    \item{input_data}{A 
	            \link[sftime]{sftime}, \link[spacetime]{STI}, \link[spacetime]{STS} or \link[spacetime]{STF} object containing the data to be interpolated.}
    \item{new_data}{A \link[sftime]{sftime}, \link[spacetime]{STI} or \link[spacetime]{STF} object containing the prediction locations.
                Must not contain NA's. If this object is not provided
				a default is calculated. This is done by taking the convex hull of \code{input_data} and
				placing around 3000 gridcells in that convex hull.}
    \item{data_variogram}{An optional way to provide a different dataset for
                the building of the variogram then for the spatial
                interpolation. }
    \item{block}{Use this parameter to pass on a specification for the
                block size. e.g. c(1000,1000) }
	\item{model}{List of models that will be tested during automatic variogram fitting. Default values are "Sph", "Exp", "Gau", and "Ste"}
    \item{kappa}{List of values for the smoothing parameter of the Matern model that will be tested during automatic variogram fitting.}
    \item{fix.values}{Can be used to fix a variogram parameter to a certain value. It
                 consists of a list with a length of three. The items describe the
                 fixed value for the nugget, range and sill respectively. Setting
                 the value to NA means that the value is not fixed. Is passed on to autofitVariogram.}
    \item{newdata_mode}{How the new data will be generated in shape. One of "rect" (rectangle) and "hull" (convex hull)}
    \item{newdata_npoints}{The number of points will be generated in the extent of the new data}
	\item{GLS.model}{If a variogram model is passed on through this parameter a Generalized Least Squares
				 sample variogram is calculated.}
    \item{tlags}{The range of time lags for fitting STVariogram}
    \item{cutoff}{The maximum range of spatial lags for fitting STVariogram}
    \item{width}{The interval for fitting spatial part of the STVariogram}
    \item{predict_chunk}{The number of chunks to predict values in `new_data`. If this value is not `NULL`, the new data will be split into chunks in size of `predict_chunk` and the prediction will be done per chunk.}
    \item{nmax}{The maximum number of spatiotemporal neighborhood to make predictions (not stable)}
    \item{aniso_method}{The method to estimate the spatiotemporal anisotropy (one of \code{linear}, \code{range}, \code{vgm}, or \code{metric})}
    \item{type_joint}{The type of theoretical model of the joint spatiotemporal variogram. Only applied when \code{joint} type of spatiotemporal variogram is chosen.}
    \item{prodsum_k}{k value for \code{stModel} in \link[gstat]{vgmST} is \code{productSum}. See \link[gstat]{vgmST} for more.}
    \item{start_vals}{Can be used to give the starting values for the variogram fitting. The items describe the
                 fixed value for the nugget, range and sill respectively. They need to be given in that order.
                 Setting the value to NA means that the value will be automatically chosen.}
    \item{miscFitOptions}{Additional options to set the behavior of \link{autofitVariogram}. For details see the
                 documentation of \link{autofitVariogram}.}
    \item{measurement_error}{integer vector (3). Adds measurement error components for spatial, temporal, and joint spatiotemporal variogram models, respectively. IT IS HIGHLY EXPERIMENTAL. MAY RESULT IN ERRORS.}
    \item{cores}{The number of cores to be used for estimating \link[gstat]{variogramST}. See \link[gstat]{variogramST} for more detail.}
  	\item{verbose}{logical, if TRUE autoKrige will give extra information on the fitting process. Default is TRUE.}
}
\details{
\code{autoKrigeST} calls the function \code{autofitVariogramST} that fits a spatiotemporal variogram model to the
given dataset. This variogram model and the data are used to make predictions on the spatiotemporal locations
in \code{new_data}. If \code{new_data} is not specified, an internal function will automatically generate a new \code{ST*DF} data to perform the spatiotemporal interpolation.
}

\value{This function returns an \code{autoKrige} object containing the results of the interpolation
(prediction, variance and standard deviation), the sample variogram, the variogram model that
was fitted by \code{autofitVariogram} and the sums of squares between the sample variogram and the
fitted variogram model. The attribute names are \code{krige_output}, \code{exp_var}, \code{var_model}
and \code{sserr} respectively.
}
\note{\code{cutoff} and \code{width} should be respectfully chosen to the scale of the data.
      The \code{width} should be larger than the minimum distance between the points in the dataset.
      The \code{cutoff} should be larger than the \code{width}}
\author{Insang Song, \email{sigmafelix@hotmail.com}}
\seealso{\code{\link{autofitVariogramST}}, \code{\link[gstat]{krigeST}}}
\examples{
# The first part of the example is from the example of krigeST
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
                       type_stv = 'sumMetric',
                       tlags = 0:7,
                       model = c('Mat', 'Ste', 'Wav', 'Exp', 'Exc'),
                       cutoff = 2e5,
                       width = 2e4,
                       cores = 4)

plot(st_as_stars(air_stk[[1]])[1,], pch = 19, cex = 0.5)

}