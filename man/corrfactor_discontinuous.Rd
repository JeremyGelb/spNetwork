% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{corrfactor_discontinuous}
\alias{corrfactor_discontinuous}
\title{A function to calculate the necessary informations to apply the
Diggle correction factor with a discontinuous method}
\usage{
corrfactor_discontinuous(neighbour_list, events, line_list, bws, max_depth)
}
\arguments{
\item{neighbour_list}{a list of the neighbours of each node}

\item{events}{a numeric vector of the node id of each event}

\item{line_list}{a DataFrame representing the lines of the graph}

\item{bws}{the kernel bandwidth for each event}

\item{max_depth}{the maximum recursion depth (after which recursion is stopped)}
}
\value{
a list of dataframes, used to calculate the Diggel correction factor
}
\description{
A function to calculate the necessary informations to apply the
Diggle correction factor with a discontinuous method
}
\keyword{internal}
