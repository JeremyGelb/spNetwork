#' @title worker function for K-nearest points on network
#'
#' @description The worker the K-nearest points for a set of points on a network.
#'
#' @param points A SpatialPointsDataFrame, for each point, its k nearest
#' neighbours will be found on the network.
#' @param lines A SpatialLinesDataFrame representing the network
#' @param k An integer indicating the number of neighbours to find..
#' @param direction Indicates a field providing information about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @param use_dest A boolean indicating if the origins and separations are
#' separated (TRUE), FALSE if only origins are used.
#' @param verbose A Boolean indicating if the function should print its
#' progress
#' @param digits The number of digits to retain in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A list with two matrices, one with the index of the neighbours and
#' one with the distances.
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom rgeos gLength
#' @examples
#' #no example provided, this is an internal function
network_knn_worker <- function(points, lines, k, direction = NULL, use_dest = FALSE, verbose = verbose, digits = digits, tol=tol){

  ## step1 : adding the points to the lines
  lines$worker_id <- 1:nrow(lines)
  points$nearest_line_id <- as.numeric(as.character(points$nearest_line_id))
  joined <- data.table(points@data)
  B <- data.table(lines@data)
  joined[B,  on = c("nearest_line_id" = "tmpid"),
         names(B) := mget(paste0("i.", names(B)))]

  ## step2 : adding the points to the lines
  if(verbose){
    print("adding the points as vertices to nearest lines")
  }
  graph_lines <- split_lines_at_vertex(lines, points, joined$worker_id, 1)
  graph_lines$lx_length <- gLength(graph_lines,byid = TRUE)

  ## step4 : building the graph
  if(verbose){
    print("build the local graph")
  }
  graph_lines$lx_weight <- (graph_lines$lx_length / graph_lines$line_length) * graph_lines$line_weight

  if (is.null(direction)){
    result_graph <- build_graph(graph_lines, digits = digits,
                                attrs = TRUE, line_weight = "lx_weight")
  }else{
    #dir <- ifelse(graph_lines[[direction]]=="Both",0,1)
    #graph_lines$direction <- graph_lines[[direction]]

    result_graph <- build_graph_directed(graph_lines, digits = digits,
                                         attrs = TRUE, line_weight='lx_weight',
                                         direction = direction)
  }

  points$vertex <- closest_points(points,result_graph$spvertices)

  if(use_dest) {
    endV <- unique(subset(points, points$type == "destination")$vertex)
    startV <- unique(subset(points, points$pttype == "start" & points$type == "origin")$vertex)
  }else{
    endV <- unique(points$vertex)
    startV <- unique(subset(points, points$pttype == "start")$vertex)
  }


  ## step4 : calculating the distances between each vertex
  if(verbose){
    print("calculating the distances on the graph")
  }
  distmat <- igraph::distances(result_graph$graph,
                               mode = "out", v = startV, to = endV)

  rownames(distmat) <- startV
  colnames(distmat) <- endV
  points$ch_vertex <- as.character(points$vertex)

  ## step5 find for each observation its n nearest neighbours
  ok_points <- subset(points, points$type == "origin")
  values <- lapply(1:nrow(ok_points), function(i){
    row <- ok_points@data[i,]
    vert <- row$ch_vertex
    dists <- distmat[vert,]
    bests <- sort(dists)[1:(k+1)]
    sub <- subset(points, points$ch_vertex %in% names(bests) & points$oids != row$oids)
    distsf <- bests[sub$ch_vertex]
    fids <- sub$oids
    fids <- fids[order(distsf)]
    distsf <- distsf[order(distsf)]
    n <- length(fids)
    if(n < k){
      fids <- c(fids, rep(NA,(k-n-1)))
      distsf <- c(distsf, rep(NA,(k-n-1)))
    }

    return(list(fids, distsf))
  })

  ## creating matrices
  matdists <- do.call(rbind,lapply(values, function(l){
    return(l[[2]])
  }))
  matoids <- do.call(rbind,lapply(values, function(l){
    return(l[[1]])
  }))

  rownames(matdists) <- ok_points$oids
  rownames(matoids) <- ok_points$oids

  return (list(matdists, matoids))
}



#' @title K-nearest points on network
#'
#' @description Calculate the K-nearest points for a set of points on a network.
#'
#' @param origins A SpatialPointsDataFrame, for each point, its k nearest
#' neighbours will be found on the network.
#' @param lines A SpatialLinesDataFrame representing the network
#' @param k An integer indicating the number of neighbours to find.
#' @param destinations A SpatialPointsDataFrame, might be used if the neighbours
#' must be found in a separate dataset. NULL if the neighbours must be found in
#' origins.
#' @param maxdistance The maximum distance between two observations to
#' consider them as neighbours. It is usefull only if a grid is used, a
#' lower value will reduce calculating time, but one must be sure that the
#' k nearest neighbours are within this radius. Otherwise NAs will be present
#' in the final matrices.
#' @param snap_dist The maximum distance to snap the start and end points on
#' the network.
#' @param line_weight The weighting to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional to the geographical length of the lines.
#' @param direction Indicates a field providing information about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(1,1), so all the calculation is
#' done in one go. It might be necessary to split it if the dataset is large.
#' @param verbose A Boolean indicating if the function should print its
#' progress
#' @param digits The number of digits to retain in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A list with two matrices, one with the index of the neighbours and
#' one with the distances.
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gLength gBuffer gIntersects gPointOnSurface
#' @export
#' @examples
#' \donttest{
#'     networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#'     eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#'     main_network_mtl <- rgdal::readOGR(networkgpkg,layer="main_network_mtl", verbose=FALSE)
#'     mtl_libraries <- rgdal::readOGR(eventsgpkg,layer="mtl_libraries", verbose=FALSE)
#'     results <- network_knn(mtl_libraries, main_network_mtl,
#'         k = 3, maxdistance = 1000, line_weight = "length",
#'         grid_shape=c(1,1), verbose = FALSE)
#' }
network_knn <- function(origins, lines, k, destinations = NULL, maxdistance = 0, snap_dist=Inf, line_weight = "length", direction=NULL, grid_shape=c(1,1), verbose = FALSE, digits = 3, tol=0.1){

  ## quick sanity check before starting
  sanity_check_knn(origins = origins, destinations = destinations,
                   lines = lines, k = k, maxdistance = maxdistance,
                   snap_dist = snap_dist, line_weight = line_weight,
                   direction = direction, grid_shape = grid_shape,
                   verbose = verbose, digits = digits, tol = tol)

  ## step0 adjusting the weights of the lines
  lines$tmpid <- 1:nrow(lines)
  lines$line_length <- gLength(lines,byid = TRUE)
  if(line_weight=="length"){
    lines$line_weight <- gLength(lines,byid = TRUE)
  }else {
    lines$line_weight <- lines[[line_weight]]
  }
  if(min(lines$line_weight)<=0){
    stop("the weights of the lines must be superior to 0")
  }

  ## step1 adjusting the directions of the lines
  if(is.null(direction) == FALSE){
    lines <- lines_direction(lines,direction)
  }

  ## step2 snap points on the lines
  if(is.null(destinations)){
    use_dest <- FALSE
    origins$type <- "origin"
    origins$base_oid <- 1:nrow(origins)
    comb_pts <- origins[c("base_oid","type")]
  }else{
    use_dest <- TRUE
    destinations$base_oid <- 1:nrow(destinations)
    destinations$type <- "destination"
    origins$base_oid <- 1:nrow(origins)
    origins$type <- "origin"
    comb_pts <- rbind(origins[c("base_oid","type")], destinations[c("base_oid","type")])
  }

  if(verbose){
    print("snapping the points to the lines (only once)")
  }
  snapped_points <- snapPointsToLines2(comb_pts,lines, idField="tmpid")

  ## step 3 building grid
  grid <- build_grid(grid_shape,list(comb_pts,lines))
  if(verbose){
    print("preparing the data in the grid")
  }
  ids <- 1:length(grid)
  list_elements <- prepare_elements_netlistw(ids,grid,snapped_points,lines,maxdistance)

  ## step7 iterating over the grid
  listvalues <- lapply(1:length(grid),function(i){
    quadra <- grid[i,]
    if(verbose){
      print(paste("working on quadra : ",i,"/",length(grid),sep=""))
    }
    elements <- list_elements[[i]]
    if(length(elements)==0){
      return()
    }else {
      all_pts <- elements[[1]]
      selected_lines <- elements[[2]]
      #calculating the elements
      values <- network_knn_worker(all_pts, selected_lines, k, direction=direction,
                                    use_dest = use_dest,
                                    verbose = verbose, digits = digits, tol=tol)
      return(values)
    }
  })

  ## step8 combining the results in two global matrices
  okvalues <- listvalues[lengths(listvalues) != 0]
  matdists <- do.call(rbind, lapply(okvalues, function(l){
    return(l[[1]])
  }))
  matoids <- do.call(rbind, lapply(okvalues, function(l){
    return(l[[2]])
  }))
  matoids <- matoids[order(as.numeric(rownames(matoids))),]
  matdists <- matdists[order(as.numeric(rownames(matdists))),]

  ## dealing with special cases
  if(k == 1){
    matoids <- matrix(matoids, ncol = 1)
    matdists <- matrix(matdists, ncol = 1)
  }
  if (use_dest){
    matoids <- matoids - (nrow(origins))
  }

  return(list("distances" = matdists,
              "ids" = matoids))
}



#' @title K-nearest points on network (multicore version)
#'
#' @description Calculate the K-nearest points for a set of points on a network with multicore support.
#'
#' @param origins A SpatialPointsDataFrame, for each point, its k nearest
#' neighbours will be found on the network.
#' @param lines A SpatialLinesDataFrame representing the network
#' @param k An integer indicating the number of neighbours to find.
#' @param destinations A SpatialPointsDataFrame, might be used if the neighbours
#' must be found in a separate dataset. NULL if the neighbours must be found in
#' origins.
#' @param maxdistance The maximum distance between two observations to
#' consider them as neighbours. It is usefull only if a grid is used, a
#' lower value will reduce calculating time, but one must be sure that the
#' k nearest neighbours are within this radius. Otherwise NAs will be present
#' in the final matrices.
#' @param snap_dist The maximum distance to snap the start and end points on
#' the network.
#' @param line_weight The weighting to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional to the geographical length of the lines.
#' @param direction Indicates a field providing information about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(1,1), so all the calculation is
#' done in one go. It might be necessary to split it if the dataset is large.
#' @param verbose A Boolean indicating if the function should print its
#' progress
#' @param digits The number of digits to retain in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A list with two matrices, one with the index of the neighbours and
#' one with the distances.
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gLength gBuffer gIntersects gPointOnSurface
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' main_network_mtl <- rgdal::readOGR(networkgpkg,layer="main_network_mtl", verbose=FALSE)
#' mtl_libraries <- rgdal::readOGR(eventsgpkg,layer="mtl_libraries", verbose=FALSE)
#' future::plan(future::multisession(workers=2))
#' results <- network_knn.mc(mtl_libraries, main_network_mtl,
#'     k = 3, maxdistance = 1000, line_weight = "length",
#'     grid_shape=c(1,1), verbose = FALSE)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
network_knn.mc <- function(origins, lines, k, destinations = NULL, maxdistance = 0, snap_dist = Inf, line_weight = "length", direction=NULL, grid_shape=c(1,1), verbose = FALSE, digits = 3, tol=0.1){

  ## quick sanity check before starting
  sanity_check_knn(origins, destinations,
                   lines, k, maxdistance, snap_dist,
                   line_weight, direction, grid_shape,
                   verbose, digits, tol)

  ## step0 adjusting the weights of the lines
  lines$tmpid <- 1:nrow(lines)
  lines$line_length <- gLength(lines,byid = TRUE)
  if(line_weight=="length"){
    lines$line_weight <- gLength(lines,byid = TRUE)
  }else {
    lines$line_weight <- lines[[line_weight]]
  }
  if(min(lines$line_weight)<=0){
    stop("the weights of the lines must be superior to 0")
  }

  ## step1 adjusting the directions of the lines
  if(is.null(direction) == FALSE){
    lines <- lines_direction(lines,direction)
  }

  ## step2 snap points on the lines
  if(is.null(destinations)){
    use_dest <- FALSE
    origins$type <- "origin"
    origins$base_oid <- 1:nrow(origins)
    comb_pts <- origins[c("base_oid","type")]
  }else{
    use_dest <- TRUE
    destinations$base_oid <- 1:nrow(destinations)
    destinations$type <- "destination"
    origins$base_oid <- 1:nrow(origins)
    origins$type <- "origin"
    comb_pts <- rbind(origins[c("base_oid","type")], destinations[c("base_oid","type")])
  }

  if(verbose){
    print("snapping the points to the lines (only once)")
  }
  snapped_points <- snapPointsToLines2(comb_pts,lines, idField="tmpid")

  ## step 3 building grid
  grid <- build_grid(grid_shape,list(comb_pts,lines))
  if(verbose){
    print("preparing the data in the grid")
  }

  all_is <- 1:length(grid)
  iseq <- list()
  cnt <- 0
  for(i in 1:grid_shape[[1]]){
    start <- cnt*grid_shape[[2]]+1
    iseq[[length(iseq)+1]] <- list(cnt+1,all_is[start:(start+grid_shape[[2]]-1)])
    cnt<-cnt+1
  }
  listelements <- future.apply::future_lapply(iseq,function(ii){
    elements <- prepare_elements_netlistw(ii[[2]],grid,snapped_points,lines,maxdistance)
    return(elements)
  })

  listelements <- unlist(listelements,recursive = FALSE)

  ## step7 iterating over the grid
  listvalues <- future.apply::future_lapply(listelements,function(elements){
    ##step1 : preparing elements
    if(is.null(elements)){
      return()
    }else {
      all_pts <- elements[[1]]
      selected_lines <- elements[[2]]
      #calculating the elements
      values <- network_knn_worker(all_pts, selected_lines, k, direction = direction,
                                   use_dest = use_dest,
                                   verbose = verbose, digits = digits, tol=tol)
      return(values)
    }

  })

  ## step8 combining the results in two global matrices
  okvalues <- listvalues[lengths(listvalues) != 0]
  matdists <- do.call(rbind, lapply(okvalues, function(l){
    return(l[[1]])
  }))
  matoids <- do.call(rbind, lapply(okvalues, function(l){
    return(l[[2]])
  }))
  matoids <- matoids[order(as.numeric(rownames(matoids))),]
  matdists <- matdists[order(as.numeric(rownames(matdists))),]

  ## dealing with special cases
  if(k == 1){
    matoids <- matrix(matoids, ncol = 1)
    matdists <- matrix(matdists, ncol = 1)
  }
  if (use_dest){
    matoids <- matoids - (nrow(origins))
  }

  return(list("distances" = matdists,
              "ids" = matoids))
}


#' @title Sanity check for the knn functions
#'
#' @description Check if all the parameters are valid for the knn functions
#'
#' @param origins A SpatialPointsDataFrame, for each point, its k nearest
#' neighbours will be found on the network.
#' @param destinations A SpatialPointsDataFrame, might be used if the neighbours
#' must be found in a separate dataset. NULL if the neighbours must be found in
#' origins.
#' @param lines A SpatialLinesDataFrame representing the network
#' @param k An integer indicating the number of neighbours to find..
#' @param maxdistance The maximum distance between two observations to
#' consider them as neighbours. It is usefull only if a grid is used, a
#' lower value will reduce calculating time, but one must be sure that the
#' k nearest neighbours are within this radius. Otherwise NAs will be present
#' in the final matrices.
#' @param snap_dist The maximum distance to snap the start and end points on
#' the network.
#' @param line_weight The weighting to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional to the geographical length of the lines.
#' @param direction Indicates a field providing information about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(1,1), so all the calculation is
#' done in one go. It might be necessary to split it if the dataset is large.
#' @param verbose A Boolean indicating if the function should print its
#' progress
#' @param digits The number of digits to retain in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A list with two matrices, one with the index of the neighbours and
#' one with the distances.
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gLength gBuffer gIntersects gPointOnSurface
#' @examples
#' #no example provided, this is an internal function
sanity_check_knn <- function(origins, destinations, lines, k, maxdistance, snap_dist, line_weight, direction, grid_shape, verbose , digits, tol){ # nocov start

  ## check des types destinations, lines et origins
  if(class(origins)[[1]] != "SpatialPointsDataFrame"){
    stop("The origins must be a SpatialPointsDataFrame")
  }
  if(is.null(destinations) == FALSE){
    if(class(destinations)[[1]] != "SpatialPointsDataFrame"){
      stop("The destinations must be NULL or a SpatialPointsDataFrame")
    }
  }
  if(class(lines)[[1]] != "SpatialLinesDataFrame"){
    stop("The lines must be a SpatialLinesDataFrame")
  }

  ## check des types numerics
  if((all.equal(k, as.integer(k))) == FALSE | k < 1){
    stop("The k parameter must be an integer >=1 ")
  }

  if(is.numeric(maxdistance) == FALSE | maxdistance < 0){
    stop("The maxdistance parameter must be a postive numeric")
  }

  if(is.numeric(snap_dist) == FALSE | snap_dist < 0){
    stop("The snap_dist parameter must be a postive numeric")
  }

  if(line_weight != "length" & line_weight %in% names(lines) == FALSE){
    stop("line_weight must be 'length' or a column in lines")
  }
  if(line_weight != "length"){
    if(is.numeric(lines[[line_weight]]) == FALSE){
      stop("line_weight must a numeric column in lines")
    }
  }

  ## check de directions
  if(is.null(direction) == FALSE){
    if(is.null(lines[[direction]])){
      stop("direction must be the name of a column in lines")
    }else{
      vals <-lines[[direction]]
      diffs <- union(unique(vals), c("TF","FT","Both"))
      if(length(diffs) != 3){
        stop("the values in the column direction must be in c('TF', 'FT', 'Both')")
      }
    }
  }

  ## check grid shape
  if(all(grid_shape == floor(grid_shape)) == FALSE){
    stop("the values in grid shape must all be integers")
  }
  if(length(grid_shape)!=2){
    stop("grid_shape must be a vector with a length of 2")
  }

  ## check the final parameters
  if(verbose %in% c(TRUE,FALSE) == FALSE){
    stop("verbose must be a boolean")
  }

  if((all.equal(digits, as.integer(digits))) == FALSE){
    stop("digits must be an integer")
  }

  if(is.numeric(tol) == FALSE){
    stop("tol must be a numeric")
  }

} # nocov end
