% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nkde_execution_functions.R
\name{split_by_grid}
\alias{split_by_grid}
\title{Split data with a grid}
\usage{
split_by_grid(grid, samples, events, lines, bw, tol, digits, split_all = TRUE)
}
\arguments{
\item{grid}{A spatial grid to split the data within}

\item{samples}{A spatialPointsDataFrame of the samples points}

\item{events}{A spatialPointsDataFrame of the events points}

\item{lines}{A SpatialLinesDataFrame representing the network}

\item{bw}{The kernel bandwidth (used to avoid edge effect)}

\item{tol}{A float indicating the spatial tolerance when snapping events on
lines}

\item{digits}{The number of digits to keep}

\item{split_all}{A boolean indicating if we must split the lines at each vertex
(TRUE) or only at event vertices (FALSE)}
}
\value{
A list with the split dataset
}
\description{
Function to split the dataset according to a grid.
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}
