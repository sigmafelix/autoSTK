% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.autofitVariogram.R
\name{plot.autofitVariogram}
\alias{plot.autofitVariogram}
\title{Plot the automatically fitted variogram}
\usage{
plot.autofitVariogram(
  x,
  plotit = TRUE,
  title = "Experimental variogram and fitted variogram model",
  ...
)
}
\arguments{
\item{x}{A result object of autofitVariogram.}

\item{plotit}{boolean. Print graph or not.}

\item{title}{character. the title of the plot.}

\item{`...`}{passed to xyplot}
}
\value{
A lattice::xyplot object.
}
\description{
Plot the automatically fitted variogram
}
