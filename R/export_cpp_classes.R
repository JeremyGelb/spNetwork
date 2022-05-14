#' @useDynLib spNetwork, .registration = TRUE
#' @import methods Rcpp
"_PACKAGE"


# Exporting spatial_index class
#' @export spatial_index

Rcpp::loadModule(module = "spatial_index_cpp", TRUE)

#' @name spatial_index
#' @title An object (c++ pointer) to do some spatial query of rectanlges (internal)
#' @description An object (c++ pointer) that uses the geometry
#' index library from boost to perform spatial queries
#' @param x - A numeric matrix with 4 columns (minX, minY, maxX, maxY)
#' @return a new instance of the spatial_index class (c++ pointer)
NULL

#' @name spatial_index$new
#' @title Constructor method for a spatial_index object
#' @description An object (c++ pointer) that uses the geometry
#' index library from boost to perform spatial queries
#' @param x - A numeric matrix with 4 columns (minX, minY, maxX, maxY)
#' @return a new instance of the spatial_index class (c++ pointer)
NULL

#' @name spatial_index$tree_request
#' @title spatial request on rtree index (internal)
#' @description A method to get the boxes in the rtree intersecting another box
#' @param reqBbox - A numeric vector with 4 values (minX, minY, maxX, maxY)
#' @return an IntegerVector with the indices of the intersected boxes
NULL
