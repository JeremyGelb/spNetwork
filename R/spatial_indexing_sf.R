# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### spatial indexing ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Obtain all the bounding boxes of a feature collection
#'
#' @description Obtain all the bounding boxes of a feature collection (INTERNAL).
#'
#' @param x a feature collection
#' @return a matrix (xmin, ymin, xmax, ymax)
#' @importFrom sf st_geometry st_as_sfc st_bbox
#' @examples
#' #This is an internal function, no example provided
st_bbox_by_feature = function(x) {
  t(sapply(st_geometry(x), function(y){st_bbox(y)}))
}


#' @title Build a quadtree
#'
#' @description Generate a quadtree object from package SearchTrees, useful to speed up
#' spatial requesting (INTERNAL).
#'
#' @param data a feature collection of linestrings or a feature collection of points
#' @return a spatial_index object (pointer to a c++ instance)
#' @export
#' @examples
#' data(mtl_network)
#' tree <- build_quadtree(mtl_network)
#' buff <- sf::st_buffer(mtl_network[55,], 50)
#' selection <- spatial_request(buff, tree, mtl_network)
build_quadtree <- function(data){
  #step1 : extracting the bbox of the geometrie
  if(class(data)[[1]] != "sf"){
    stop("The data argument must be a feature collection from the package sf")
  }

  boxes <- st_bbox_by_feature(data)
  spIndex <- new(spatial_index, boxes)

  return(spIndex)
}


# build_quadtree <- function(data){
#   #step1 : extracting the bbox of the geometrie
#   if(class(data)[[1]] != "sf"){
#     stop("The data argument must be a feature collection from the package sf")
#   }
#   geom_type <- unique(st_geometry_type(data))
#
#   if(geom_type == "LINESTRING"){
#
#     coords <- st_coordinates(data)
#     ids <- coords[,3]
#     coords_list <- split(coords[,1:2], f = ids)
#
#     bbox_coords <- lapply(coords_list,function(i){
#       line_coords <- matrix(i, ncol = 2, byrow = FALSE)
#       maxs <- apply(line_coords,MARGIN=2,FUN = max)
#       mins <- apply(line_coords,MARGIN=2,FUN = min)
#       row <- c(maxs,mins)
#       return(row)
#     })
#
#     bbox_coords <- do.call(rbind,bbox_coords)
#     #step2 : generate the spatial index
#     spIndex <- SearchTrees::createTree(bbox_coords,dataType = "rect")
#
#   }else if (geom_type == "POINT"){
#     coords <- st_coordinates(data)
#     spIndex <- SearchTrees::createTree(coords,dataType = "point")
#   }else {
#     stop("The supported geometry types are POINT and LINESTRING")
#   }
#
#   return(spIndex)
# }


#' @title Spatial request
#'
#' @description Use a quadtree index to perform spatial request.
#'
#' @param geometry sf like object (feature collection or simple geometry)
#' @param tree a tree object from package SearchTrees
#' @param data the original data used to build the tree object
#' @return a subset of data, intersecting geometry
#' @importFrom sf st_intersects
#' @export
#' @examples
#' data(mtl_network)
#' tree <- build_quadtree(mtl_network)
#' buff <- sf::st_buffer(mtl_network[55,], 50)
#' selection <- spatial_request(buff, tree, mtl_network)
spatial_request <- function(geometry,tree,data){
  ## step1 : find candidates
  #box <- t(raster::bbox(geometry))
  box <- st_bbox(geometry)
  idx <- tree$tree_request(box)
  candidates <- data[idx,]
  if(nrow(candidates) > 0){
    ## step2 : find real intersection
    final_vector <- st_intersects(candidates, geometry, sparse = FALSE)[,1]
    final_data <- subset(candidates,final_vector)
    return(final_data)
  }else{
    return(candidates)
  }
}


#' @title Find closest points
#'
#' @description Solve the nearest neighbour problem for two feature collections of points
#' This is a simple wrap-up of the dbscan::kNN function
#'
#' @param origins a feature collection of points
#' @param targets a feature collection of points
#' @return for each origin point, the index of the nearest target point
#' @export
#' @examples
#' data(mtl_libraries)
#' data(mtl_theatres)
#' close_libs <- closest_points(mtl_theatres, mtl_libraries)
closest_points <- function(origins, targets){

  xy_origins <- st_coordinates(origins)
  xy_targets <- st_coordinates(targets)

  # idx <- FNN::knnx.index(data = xy_targets,
  #                       query = xy_origins,
  #                       k = 1)
  if(nrow(xy_targets) > 1){
    idx <- dbscan::kNN(x = xy_targets, query = xy_origins, k = 1)$id
  }else if(nrow(xy_targets) == 1){
    return(rep(1, nrow(xy_origins)))
  }else{
    stop("Error in the function closest_points, less than one target...")
  }

  return(as.vector(idx[,1]))
}




