% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NEW_k_functions_sf.R
\name{kfunctions}
\alias{kfunctions}
\title{Network k and g functions (maturing)}
\usage{
kfunctions(
  lines,
  points,
  start,
  end,
  step,
  width,
  nsim,
  conf_int = 0.05,
  digits = 2,
  tol = 0.1,
  agg = NULL,
  verbose = TRUE,
  return_sims = FALSE,
  calc_g_func = TRUE,
  resolution = NULL
)
}
\arguments{
\item{lines}{A feature collection of linestrings representing the underlying network. The
geometries must be simple Linestrings (may crash if some geometries
are invalid) without MultiLineSring}

\item{points}{A feature collection of points representing the points on the
network. These points will be snapped on their nearest line}

\item{start}{A double, the lowest distance used to evaluate the k and g functions}

\item{end}{A double, the highest distance used to evaluate  the k and g functions}

\item{step}{A double, the step between two evaluations of the k and g
function. start, end and step are used to create a vector of distances with the function seq}

\item{width}{The width of each donut for the g-function. Half of the width is applied on
both sides of the considered distance}

\item{nsim}{An integer indicating the number of Monte Carlo simulations
to perform for inference}

\item{conf_int}{A double indicating the width confidence interval (default =
0.05) calculated on the Monte Carlo simulations}

\item{digits}{An integer indicating the number of digits to retain from the
spatial coordinates}

\item{tol}{When adding the points to the network, specify the minimum
distance between these points and the lines' extremities. When points are
closer, they are added at the extremity of the lines}

\item{agg}{A double indicating if the events must be aggregated within a distance.
If NULL, the events are aggregated only by rounding the coordinates}

\item{verbose}{A Boolean indicating if progress messages should be displayed}

\item{return_sims}{a boolean indicating if the simulated k and g values must also
be returned.}

\item{calc_g_func}{A Boolean indicating if the G function must also be calculated (TRUE by default). If FALSE, then only the K function is calculated}

\item{resolution}{When simulating random points on the network, selecting a
resolution will reduce greatly the calculation time. When resolution is null
the random points can occur everywhere on the graph. If a value is specified,
the edges are split according to this value and the random points can only be
vertices on the new network}
}
\value{
A list with the following values :
\itemize{
\item plotk: A ggplot2 object representing the values of the k-function
\item plotg: A ggplot2 object representing the values of the g-function
\item values: A DataFrame with the values used to build the plots
}
}
\description{
Calculate the k and g functions for a set of points on a
network (maturing).
}
\details{
The k-function is a method to characterize the dispersion of a set
of points. For each point, the numbers of other points in subsequent radii
are calculated. This empirical k-function can be more or less clustered
than a k-function obtained if the points were randomly located in space. In
a network, the network distance is used instead of the Euclidean distance.
This function uses Monte Carlo simulations to assess if the points are
clustered or dispersed, and gives the results as a line plot. If the line
of the observed k-function is higher than the shaded area representing the
values of the simulations, then the points are more clustered than what we
can expect from randomness and vice-versa. The function also calculates the
g-function, a modified version of the k-function using rings instead of
disks. The width of the ring must be chosen. The main interest is to avoid
the cumulative effect of the classical k-function. This function is maturing,
it works as expected (unit tests) but will probably be modified in the
future releases (gain speed, advanced features, etc.).
}
\examples{
\donttest{
data(main_network_mtl)
data(mtl_libraries)
result <- kfunctions(main_network_mtl, mtl_libraries,
     start = 0, end = 2500, step = 100,
     width = 200, nsim = 50,
     conf_int = 0.05, tol = 0.1, agg = NULL,
     calc_g_func = TRUE,
     verbose = FALSE)
}
}
