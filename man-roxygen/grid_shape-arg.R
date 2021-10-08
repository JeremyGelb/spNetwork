#' @param grid_shape A vector of two values indicating how the study area
#' must be split when performing the calculus. Default is c(1,1) (no split). A finer grid could
#' reduce memory usage and increase speed when a large dataset is used. When using
#' multiprocessing, the work in each grid is dispatched between the workers.
