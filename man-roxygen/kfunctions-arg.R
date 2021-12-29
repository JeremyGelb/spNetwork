#' @param lines A feature collection of linestrings representing the underlying network. The
#' geometries must be simple Linestrings (may crash if some geometries
#'  are invalid) without MultiLineSring
#' @param points A feature collection of points representing the points on the
#'   network. These points will be snapped on their nearest line
