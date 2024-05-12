#' @useDynLib spNetwork, .registration = TRUE
#' @import methods Rcpp
"_PACKAGE"


Rcpp::loadModule(module = "spatial_index_cpp", TRUE)

#' @name spatial_index
#' @title An object (c++ pointer) to do some spatial query of rectanlges (internal)
#' @description An object (c++ pointer) that uses the geometry
#' index library from boost to perform spatial queries
#' @param x - A numeric matrix with 4 columns (minX, minY, maxX, maxY)
#' @return a new instance of the spatial_index class (c++ pointer)
#' @export
#' @examples
#' # no example provided, this is an internal class
NULL
