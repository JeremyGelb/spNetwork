#' @name spatial_index
#' @title An object (c++ pointer) to do some spatial query of rectanlges (internal)
#' @description An object (c++ pointer) that uses the geometry
#' index library from boost to perform spatial queries
#' @param x - A numeric matrix with 4 columns (minX, minY, maxX, maxY)
#' @return a new instance of the spatial_index class (c++ pointer)
#' @export spatial_index
#' @examples
#' # no example provided, this is an internal class

loadModule(module = "spatial_index_cpp", TRUE)
