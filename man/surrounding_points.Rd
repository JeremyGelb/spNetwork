% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometrical_functions.R
\name{surrounding_points}
\alias{surrounding_points}
\title{Points along polygon boundary}
\usage{
surrounding_points(polygons, dist)
}
\arguments{
\item{polygons}{A SpatialPolygonsDataFrame}

\item{dist}{The distance between the points}
}
\value{
A SpatialPolygonsDataFrame representing the grid
}
\description{
Generate a SpatialPointsDataFrame by placing points along the border of
polygons of a SpatialPolygonsDataFrame.
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}