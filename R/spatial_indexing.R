# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### spatial indexing ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Build a quadtree
#'
#' @description Generate quadtree object from package SearchTrees, useful to speed up
#' spatial requesting.
#'
#' @param data a SpatialLinesDataFrame or a SpatialPointsDataFrame
#' @return quadtree object from package SearchTrees
#' @examples
#' #This is an internal function, no example provided
#' @keywords internal
build_quadtree <- function(data){
  #step1 : extracting the bbox of the geometrie
  if(class(data)=="SpatialLinesDataFrame"){
    coords <- sp::coordinates(data)
    bbox_coords <- lapply(coords,function(i){
      line_coords <- do.call(rbind,i)
      maxs <- apply(line_coords,MARGIN=2,FUN = max)
      mins <- apply(line_coords,MARGIN=2,FUN = min)
      row <- c(maxs,mins)
      return(row)
    })
    bbox_coords <- do.call(rbind,bbox_coords)
    #step2 : generate the spatial index
    spIndex <- SearchTrees::createTree(bbox_coords,dataType = "rect")
  }else if (class(data)=="SpatialPointsDataFrame"){
    coords <- sp::coordinates(data)
    spIndex <- SearchTrees::createTree(coords,dataType = "point")
  }else {
    stop("The class of the data argument must be SpatialLinesDataFrame or SpatialPointsDataFrame")
  }

  return(spIndex)
}


#' @title Spatial request
#'
#' @description Use a quadtree index to perform spatial request.
#'
#' @param geometry objects such as SpatialLine, SpatialPolygon or SpatialPoint
#' @param tree a tree object from package SearchTrees
#' @param data the original data used to build the tree object
#' @return a subset of data, intersecting geometry
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
spatial_request <- function(geometry,tree,data){
  ## step1 : find candidates
  box <- t(raster::bbox(geometry))
  idx <- SearchTrees::rectLookup(tree,box[1,],box[2,])
  candidates <- data[idx,]
  ## step2 : find real intersection
  final_vector <- as.vector(rgeos::gIntersects(candidates,
                                               geometry, byid = TRUE))
  final_data <- subset(candidates,final_vector)
  return(final_data)
}


#' @title Find closest points
#'
#' @description Solve the nearest neighbour problem for two SpatialPointsDataFrame.
#' This is a simple wrapup of the FNN::knnx.index function
#'
#' @param origins a SpatialPointsDataFrame
#' @param targets a SpatialPointsDataFrame
#' @return for each origin point, the index of the nearest target point
#' @export
#' @examples
#' #This is an internal function, no example provided
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_libraries <- rgdal::readOGR(eventsgpkg,layer="mtl_libraries", verbose=FALSE)
#' mtl_theatres <- rgdal::readOGR(eventsgpkg,layer="mtl_theatres", verbose=FALSE)
#' close_libs <- closest_points(mtl_theatres, mtl_libraries)
closest_points <- function(origins, targets){

  xy_origins <- sp::coordinates(origins)
  xy_targets <- sp::coordinates(targets)

  idx <- FNN::knnx.index(data = xy_targets,
                         query = xy_origins,
                         k = 1)

  return(idx[,1])
}




