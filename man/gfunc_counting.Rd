% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gfunc_counting}
\alias{gfunc_counting}
\title{c++ g function counting worker}
\usage{
gfunc_counting(dist_mat, wc, wr, breaks, width)
}
\arguments{
\item{dist_mat}{A matrix with the distances between points}

\item{wc}{The weight of the points represented by the columns (destinations)}

\item{wr}{The weight of the points represented by the rows (origins)}

\item{breaks}{A numeric vector with the distance to consider}

\item{width}{The width of each donut}
}
\value{
A numeric matrix with the countings of the g function evaluated at the required distances
}
\description{
c++ k function counting (INTERNAL)
}
\keyword{internal}
