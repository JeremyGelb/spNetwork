#' @param start A double, the lowest distance used to evaluate the k and g functions
#' @param end A double, the highest distance used to evaluate  the k and g functions
#' @param step A double, the step between two evaluations of the k and g
#'   function. start, end and step are used to create a vector of distances with the function seq
#' @param width The width of each donut for the g-function. Half of the width is applied on
#' both sides of the considered distance
#' @param nsim An integer indicating the number of Monte Carlo simulations
#'   to perform for inference
#' @param conf_int A double indicating the width confidence interval (default =
#'   0.05) calculated on the Monte Carlo simulations
#' @param digits An integer indicating the number of digits to retain from the
#'   spatial coordinates
#' @param tol When adding the points to the network, specify the minimum
#'   distance between these points and the lines' extremities. When points are
#'   closer, they are added at the extremity of the lines
#' @param resolution When simulating random points on the network, selecting a
#'   resolution will reduce greatly the calculation time. When resolution is null
#'   the random points can occur everywhere on the graph. If a value is specified,
#'   the edges are split according to this value and the random points can only be
#'   vertices on the new network
#' @param agg A double indicating if the events must be aggregated within a distance.
#' If NULL, the events are aggregated only by rounding the coordinates
#' @param verbose A Boolean indicating if progress messages should be displayed
#' @param calc_g_func A Boolean indicating if the G function must also be calculated (TRUE by default). If FALSE, then only the K function is calculated
#' @param return_sims a boolean indicating if the simulated k and g values must also
#' be returned.