% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_functions_sf.R
\name{dist_mat_dupl}
\alias{dist_mat_dupl}
\title{Distance matrix with dupicated}
\usage{
dist_mat_dupl(graph, start, end, ...)
}
\arguments{
\item{graph}{The Graph to use}

\item{start}{The vertices to use as starting points}

\item{end}{The vertices to use as ending points}

\item{...}{parameters passed to the function igraph::distances}
}
\value{
A matrix with the distances between the vertices
}
\description{
Function to Create a distance matrix when some vertices are duplicated.
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}
