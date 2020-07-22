\name{autofitVariogramST}
\alias{autofitVariogramST}
\title{
Automatically fit a spatiotemporal variogram from ST*DF}
\description{

}
\usage{
autofitVariogramST(
                   stf,
                   formula,
                   typestv ='sumMetric',
                   candidate_model = c('Ste','Exc','Exp','Wav'),
                   guess_nugget = NULL,
                   guess_psill = NULL,
                   tlags = 0:6,
                   cutoff = 2e4,
                   width = 5e2,
                   aniso_method = 'vgm',
                   type_joint = 'Exp',
                   prodsum_k = NULL,
                   surface = FALSE,
                   cores = 1
                       )}

\arguments{
    \item{stf}{An object of the
	            \link[spacetime]{STI-class} or \link[spacetime]{STF-class} containing the data from which the spatiotemporal variogram to be estimated.}
	\item{formula}{formula that defines the dependent variable as a linear model
				of independent variables; suppose the dependent variable has
				name 'z', for ordinary and simple kriging use the formula
				'z~1'; for simple kriging also define 'beta' (see below); for
				universal kriging, suppose 'z' is linearly dependent on 'x'
				and 'y', use the formula 'z~x+y'.}
    \item{typestv}{Model for the spatiotemporal variogram}
	\item{candidate_model}{List of models that will be tested during automatic variogram fitting. Default values are "Sph", "Exp", "Gau", and "Ste"}
    \item{guess_nugget}{A user-defined value for the nugget of the spatiotemporal variogram.}
    \item{guess_psill}{A user-defined value for the partial sill of the spatiotemporal variogram}
    \item{tlags}{The range of time lags for fitting STVariogram}
    \item{cutoff}{The maximum range of spatial lags for fitting STVariogram}
    \item{width}{The interval for fitting spatial part of the STVariogram}
    \item{aniso_method}{The method to estimate the spatiotemporal anisotropy (one of \code{linear}, \code{range}, \code{vgm}, or \code{metric})}
    \item{type_joint}{The type of theoretical model of the joint spatiotemporal variogram. Only applied when \code{joint} type of spatiotemporal variogram is chosen.}
    \item{prodsum_k}{k value for \code{stModel} in \link[gstat]{vgmST} is \code{productSum}. See \link[gstat]{vgmST} for more.}
    #\item{surface}{\code{logical}. Set \code{TRUE} if you want to get the wireframe plot of the fitted spatiotemporal variogram}
    \item{cores}{The number of cores to be used for estimating \link[gstat]{variogramST}. See \link[gstat]{variogramST} for detail.}
}
}
\details{
}
\value{
\item{jointSTV}{The fitted theoretical spatiotemporal variogram}
\item{empSTV}{The empirical spatiotemporal variogram estimated from the input data}
\item{SpV}{The spatial component of \code{jointSTV}}
\item{TV}{The temporal component of \code{jointSTV}}
\item{STVsurface}{\code{data.frame}. The variogram surface of \code{jointSTV}}
}
\author{Insang Song, \email{sigmafelix@hotmail.com}}
\examples{
}
