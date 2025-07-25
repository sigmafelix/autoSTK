% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setSTI.R
\name{setSTI}
\alias{setSTI}
\title{A convenience function for the sample spatiotemporal variogram}
\usage{
setSTI(
  stf,
  formula,
  tlags = 0:6,
  cutoff = 30000,
  width = 1000,
  assumeRegular = TRUE,
  pseudo = TRUE,
  na.omit = TRUE,
  wireframe = FALSE,
  plot3d = FALSE,
  cores = 1
)
}
\arguments{
\item{stf}{A ST*DF object.}

\item{formula}{formula (inherits the same parameter in variogramST)}

\item{tlags}{temporal lags to compute semivariance (inherits the same parameter in variogramST)}

\item{cutoff}{the maximum bound of the set of spatial lags (inherits the same parameter in variogramST)}

\item{width}{integer (1). spatial lag (inherits the same parameter in variogramST)}

\item{assumeRegular}{Boolean. Assuming regular grid?}

\item{pseudo}{Boolean. See ?gstat::variogramST}

\item{na.omit}{Boolean. Omit NA values.}

\item{wireframe}{Boolean. Whether you plot a StVariogram in wireframe or not. If not, the return will be in class of data.frame, not a list}

\item{plot3d}{Boolean. Wheter you make a three-dimensional graph with rgl package}

\item{logarithm}{Boolean. log-transformation}
}
\value{
Depends on the arguments wireframe (if TRUE, list of length 2) and plot3d (if TRUE, list of length 3), a StVariogram object otherwise.
}
\description{
A convenience function for the sample spatiotemporal variogram
}
