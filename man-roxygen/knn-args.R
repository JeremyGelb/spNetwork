#' @param origins A SpatialPointsDataFrame, for each point, its k nearest
#' neighbours will be found on the network.
#' @param lines A SpatialLinesDataFrame representing the underlying network
#' @param k An integer indicating the number of neighbours to find.
#' @param destinations A SpatialPointsDataFrame, might be used if the neighbours
#' must be found in a separate set of points NULL if the neighbours must be found in
#' origins.
#' @param maxdistance The maximum distance between two observations to
#' consider them as neighbours. It is useful only if a grid is used, a
#' lower value will reduce calculating time, but one must be sure that the
#' k nearest neighbours are within this radius. Otherwise NAs will be present
#' in the results.
#' @param snap_dist The maximum distance to snap the start and end points on
#' the network.
#' @param line_weight The weighting to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional to the geographical length of the lines.
#' @param direction The name of a column indicating authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(1,1), so all the calculation is
#' done in one go. It might be necessary to split it if the dataset is large.
#' @param verbose A Boolean indicating if the function should print its
#' progress
#' @param digits The number of digits to retain from the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the minimum distance between the points and the
#'   lines' extremities when adding the point to the network. When points are
#'   closer, they are added at the extremity of the lines.
