# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### spatial indexing ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Generate quad tree object from package SearchTrees, usefull to speed up
#' spatial requesting
#'
#' @param data a SpatialLinesDataFrame or a SpatialPointsDataFrame
#' @return quad tree object from package SearchTrees
#' @examples
#' #This is an internal function, no example provided
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


#' Use a quad tree index to perform spatial request
#'
#' @param geometry object like SpatialLine, SpatialPolygon or SpatialPoint
#' @param tree a tree object from package SearchTrees
#' @param data the original data used to build the tree object
#' @return a subset of data, intersecting geometry
#' @examples
#' #This is an internal function, no example provided
spatial_request <- function(geometry,tree,data){
  ## step1 : find candidates
  box <- t(raster::bbox(geometry))
  idx <- SearchTrees::rectLookup(tree,box[1,],box[2,])
  candidates <- data[idx,]
  ## step2 : find real intersection
  final_vector <- as.vector(rgeos::gIntersects(candidates, geometry, byid = T))
  final_data <- subset(candidates,final_vector)
  return(final_data)
}


#' build a quad tree index and solve the nearest neighbour problem for two
#' SpatialPointsDataFrame
#'
#' @param origins a SpatialPointsDataFrame
#' @param targets a SpatialPointsDataFrame
#' @return for each origin point, the index of the nearest target point
#' @examples
#' #This is an internal function, no example provided
closest_points <- function(origins, targets){
  ## step1 : create the spatial index for the target points
  original_coords <- sp::coordinates(targets)
  original_coords <- cbind(original_coords,1:nrow(targets))
  tree <- SearchTrees::createTree(original_coords[,1:2],dataType = "point")
  ## step2 : performe the spatial request
  pts <- sp::coordinates(origins)
  k1 <- SearchTrees::knnLookup(tree,newdat=pts,k = 1)
  extract <- original_coords[k1,]
  if(nrow(origins)==1){
    idx <- extract[3]
  }else{
    idx <- extract[,3]
  }

  return(idx)
}

