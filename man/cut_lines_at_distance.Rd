% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometrical_functions_sf.R
\name{cut_lines_at_distance}
\alias{cut_lines_at_distance}
\title{Cut lines at a specified distance}
\usage{
cut_lines_at_distance(lines, dists)
}
\arguments{
\item{lines}{The feature collection of linestrings to cut}

\item{dists}{A vector of distances, if only one value is given,
each line will be cut at that distance.}
}
\value{
A feature collection of linestrings
}
\description{
Cut lines in a feature collection of linestrings at a specified distance from the
begining of the lines.
}
\examples{
# This is an interal function, no example provided
}
\keyword{internal}
