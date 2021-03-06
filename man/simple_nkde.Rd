% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_functions.R
\name{simple_nkde}
\alias{simple_nkde}
\title{Simple NKDE algorithm}
\usage{
simple_nkde(graph, events, samples, bws, kernel_func, nodes, edges)
}
\arguments{
\item{graph}{a graph object from igraph representing the network}

\item{events}{a SpatialPointsDataFrame representing the events. It must be
snapped on the network, and be nodes of the network. A column vertex_id
must indicate for each event its corresponding node}

\item{samples}{a SpatialPointsDataFrame representing the sampling points.
The samples must be snapped on the network. A column edge_id must indicate
for each sample on which edge it is snapped.}

\item{bws}{a vector indicating the kernel bandwidth (in meters) for each
event}

\item{kernel_func}{a function obtained with the function select_kernel}

\item{nodes}{a SpatialPointsDataFrame representing the nodes of the network}

\item{edges}{a SpatialLinesDataFrame representing the edges of the network}
}
\value{
a dataframe with two columns. sum_k is the sum for each sample point
 of the kernel values. n is the number of events influencing each sample
point
}
\description{
Function to perform the simple nkde.
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}
