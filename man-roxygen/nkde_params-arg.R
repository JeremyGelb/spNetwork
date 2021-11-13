#' @param lines A SpatialLinesDataFrame representing the underlying network. The
#' geometries must be a SpatialLinesDataFrame (may crash if some geometries
#'  are invalid) without MultiLineSring.
#' @param events events A SpatialPointsDataFrame representing the events on the
#' network. The points will be snapped on the network to their closest line.
#' @param w A vector representing the weight of each event
#' @param kernel_name The name of the kernel to use. Must be one of triangle,
#' gaussian, tricube, cosine ,triweight, quartic, epanechnikov or uniform.
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see nkde details for more information)
#' @param max_depth when using the continuous and discontinuous methods, the
#' calculation time and memory use can go wild  if the network has many
#' small edges (area with many of intersections and many events). To
#' avoid it, it is possible to set here a maximum depth. Considering that the
#' kernel is divided at intersections, a value of 10 should yield good
#' estimates in most cases. A larger value can be used without a problem for the
#' discontinuous method. For the continuous method, a larger value will
#' strongly impact calculation speed.
