% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometrical_functions.R
\name{split_at_vertices}
\alias{split_at_vertices}
\title{Split a line at vertices in a SpatialLinesDataFrame}
\usage{
split_at_vertices(line, points, i, mindist)
}
\arguments{
\item{line}{The SpatialLine to split}

\item{points}{The SpatialPoints to add to as vertex to the lines}

\item{mindist}{The minimum distance between one point and the extremity of
the line to add the point as a vertex.}

\item{nearest_lines_idx}{For each point, the index of the nearest line}
}
\value{
An object of class SpatialLinesDataFrame (package sp)
}
\description{
Split a line (SpatialLine) at their nearest vertices
(SpatialPoints), may fail if the lines geometries are self intersecting.
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}
