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
  #etape 1 : extraire les bbox des geometries
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
    #etape 2 : generer l'index spatial
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
  ##step1 : find candidates
  box <- t(raster::bbox(geometry))
  idx <- SearchTrees::rectLookup(tree,box[1,],box[2,])
  candidates <- data[idx,]
  ##step2 : find real intersection
  final_vector <- as.vector(rgeos::gIntersects(candidates, geometry, byid = T))
  final_data <- subset(candidates,final_vector)
  return(final_data)
}


#' build a quad tree index and solve the nearest neighbour problem for two
#' SpatialPointsDataFrame
#'
#' @param points a SpatialPointsDataFrame
#' @param features a SpatialPointsDataFrame
#' @return for each point, the index of the nearest line (row-wise)
#' @examples
#' #This is an internal function, no example provided
closest_points <- function(points, features){
  ## step1 : creer l'index spatial pour les li
  original_coords <- sp::coordinates(features)
  original_coords <- cbind(original_coords,1:nrow(features))
  tree <- SearchTrees::createTree(original_coords[,1:2],dataType = "point")
  ## realiser la requete a partir des points
  pts <- sp::coordinates(points)
  k1 <- SearchTrees::knnLookup(tree,newdat=pts,k = 1)
  extract <- original_coords[k1,]
  idx <- extract[,3]
  return(idx)
}

