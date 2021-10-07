#' @param digits The number of digits to retain from the spatial coordinates. It
#'   ensures that topology is good when building the network. Default is 3. A
#'   too high precision (high number of digits) might break some connections
#' @param tol A float indicating the minimum distance between the events and the
#'   lines' extremities when adding the point to the network. When points are
#'   closer, they are added at the extremity of the lines.
#' @param agg A double indicating if the events must be aggregated within a
#'   distance. If NULL, the events are aggregated only by rounding the
#'   coordinates.
