% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometrical_functions.R
\name{add_center_lines.mc}
\alias{add_center_lines.mc}
\title{Add center vertex to lines (multicore)}
\usage{
add_center_lines.mc(lines, show_progress = TRUE, chunk_size = 100)
}
\arguments{
\item{lines}{The SpatialLinesDataframe to use}

\item{show_progress}{A Boolean indicating if a progress bar must be displayed}

\item{chunk_size}{The size of a chunk used for multiprocessing. Default is
100.}
}
\value{
An object of class SpatialLinesDataframe (package sp)
}
\description{
Add to each line of a SpatialLinesDataFrame an additional vertex
 at its center with multicore support.
}
\examples{
#This is an internal function, no example provided
\dontshow{
   ## R CMD check: make sure any open connections are closed afterward
   if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
}
}
\keyword{internal}
