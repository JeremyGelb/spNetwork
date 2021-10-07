#' @param lines A SpatialLinesDataFrame representing the underlying network. The
#' geometries must be a SpatialLinesDataFrame (may crash if some geometries
#'  are invalid) without MultiLineSring
#' @param points A SpatialPointsDataFrame representing the points on the
#'   network. These points will be snapped on their nearest line
