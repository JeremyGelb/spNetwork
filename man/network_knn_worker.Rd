% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn_sf.R
\name{network_knn_worker}
\alias{network_knn_worker}
\title{worker function for K-nearest points on network}
\usage{
network_knn_worker(
  points,
  lines,
  k,
  direction = NULL,
  use_dest = FALSE,
  verbose = verbose,
  digits = digits,
  tol = tol
)
}
\arguments{
\item{points}{A feature collection of points, for each point, its k nearest
neighbours will be found on the network.}

\item{lines}{A feature collection of lines representing the network}

\item{k}{An integer indicating the number of neighbours to find..}

\item{direction}{Indicates a field providing information about authorized
travelling direction on lines. if NULL, then all lines can be used in both
directions. Must be the name of a column otherwise. The values of the
column must be "FT" (From - To), "TF" (To - From) or "Both".}

\item{use_dest}{A boolean indicating if the origins and separations are
separated (TRUE), FALSE if only origins are used.}

\item{verbose}{A Boolean indicating if the function should print its
progress}

\item{digits}{The number of digits to retain in the spatial coordinates (
simplification used to reduce risk of topological error)}

\item{tol}{A float indicating the spatial tolerance when points are
added as vertices to lines.}
}
\value{
A list with two matrices, one with the index of the neighbours and
one with the distances.
}
\description{
The worker the K-nearest points for a set of points on a network.
}
\examples{
#no example provided, this is an internal function
}
\keyword{internal}
