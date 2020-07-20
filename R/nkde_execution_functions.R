# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### helper functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function to check if the geometries given by the user are valide
#'
#' @param lines a SpatialLinesDataFrame
#' @param samples a SpatialPointsDataFrame of the samples
#' @param events a SpatialPointsDataFrame of the events
#' @param study_area a SpatialPointsDataFrame of the study_area
#' @return TRUE if all the checks are passed
#' @importFrom rgeos gIsSimple gIsValid
#' @importFrom sp is.projected
#' @examples
#' #This is an internal function, no example provided
check_geometries <- function(lines,samples,events, study_area){

  # checking if geometries are all valid, simple and planar
  obj_names <- c("lines","samples","events")
  objs <- list(lines,samples,events)
  for(i in 1:length(obj_names)){
    obj <- objs[[i]]
    obj_name <- obj_names[[i]]
    if(any(gIsSimple(obj,byid = TRUE)==FALSE)){
      stop(paste("the ",obj_name," must be simple geometries",sep=""))
    }
    if(any(gIsValid(obj,byid = TRUE)==FALSE)){
      stop(paste("the ",obj_name," must be valid geometries",sep=""))
    }
    if(is.projected(obj)==FALSE){
      stop(paste("the ",obj_name," must be projected (planar coordinates)",sep=""))
    }
  }
  # checking if the geometries type are good
  if(class(events)!="SpatialPointsDataFrame"){
    stop("the events must be given as a SpatialPointsDataFrame")
  }
  if(class(samples)!="SpatialPointsDataFrame"){
    stop("the samples must be given as a SpatialPointsDataFrame")
  }
  if(class(lines)!="SpatialLinesDataFrame"){
    stop("the lines must be given as a SpatialPointsDataFrame")
  }

  # checking if the CRS are good
  if(is.null(study_area)){
    comp <- c(raster::compareCRS(raster::crs(samples),raster::crs(events)),
              raster::compareCRS(raster::crs(lines),raster::crs(events)),
              raster::compareCRS(raster::crs(lines),raster::crs(samples)))
  }else{
    comp <- c(raster::compareCRS(raster::crs(study_area),raster::crs(events)),
              raster::compareCRS(raster::crs(study_area),raster::crs(samples)),
              raster::compareCRS(raster::crs(study_area),raster::crs(lines)),
              raster::compareCRS(raster::crs(samples),raster::crs(events)),
              raster::compareCRS(raster::crs(lines),raster::crs(events)),
              raster::compareCRS(raster::crs(lines),raster::crs(samples)))
  }

  if(any(comp==FALSE)){
    stop("the lines, events and samples must have the same Coordinates Reference System (crs)")
  }
  return(TRUE)
}

#defining some global variables (weird felx but ok)
utils::globalVariables(c("spid", "weight", "."))

#' Function to avoid having events at the same location
#'
#' @param events the SpatialPointsDataFrame to contract (must have a weight column)
#' @param digits the number of digits to keep
#' @param agg a double indicating if the points must be aggregated within a distance.
#' if NULL, then the points are aggregated by rouding the coordinates.
#' @return a new SpatialPointsDataFrame
#' @importFrom data.table tstrsplit setDF
#' @examples
#' #This is an internal function, no example provided
clean_events <- function(events,digits=5,agg=NULL){
  if(is.null(agg)){
    events$spid <- sp_char_index(raster::coordinates(events),digits)
    new_events <- data.table(events@data[c("weight","spid")])
    agg_events <- new_events[, .(sum(weight)), by = .(spid)]
    agg_events[,  c("X", "Y") := tstrsplit(spid, "_", fixed=TRUE)]
    agg_events$X <- as.numeric(agg_events$X)
    agg_events$Y <- as.numeric(agg_events$Y)
    agg_events$weight <- agg_events$V1
    new_events <- setDF(agg_events)
    new_events <- new_events[c("weight","spid","X","Y")]
    sp::coordinates(new_events) <- cbind(new_events$X,new_events$Y)
    raster::crs(new_events) <- raster::crs(events)
    return(new_events)
  }else{
    new_events <- aggregate_points(events,agg)
    new_events$spid <- sp_char_index(raster::coordinates(new_events),digits)
    return(new_events)
  }

}

#' Function to aggregate points within a radius
#'
#' @param points the SpatialPointsDataFrame to contract (must have a weight column)
#' @param maxdist the distance to use
#' @return a new SpatialPointsDataFrame
#' @importFrom rgeos gBuffer
#' @examples
#' #This is an internal function, no example provided
aggregate_points <- function(points, maxdist){
  Continue <- TRUE
  while(Continue){
    #generer l'arbre
    tree <- build_quadtree(points)
    #mettre a 0 les appartenances
    points$aggregated <- 0
    points$oid <- 1:nrow(points)
    #demarrer les iterations
    new_features <- lapply(1:nrow(points), function(i){
      row <- points[i,]
      buff <- gBuffer(row,width=maxdist)
      candidates <- spatial_request(buff,tree,points)
      ok_cand <- subset(candidates,candidates$aggregated==0)
      if(nrow(ok_cand)>0){
        coords <- sp::coordinates(ok_cand)
        points[ok_cand$oid,"aggregated"] <<- 1
        new_row <- c(sum(ok_cand$weight),mean(coords[,1]), mean(coords[,2]))
        return(new_row)
      }else{
        return(NULL)
      }
    })
    #let us remove all the empty quadra
    new_features <- new_features[lengths(new_features) != 0]
    new_points <- data.frame(do.call(rbind,new_features))
    names(new_points) <- c("weight","x","y")
    sp::coordinates(new_points) <- cbind(new_points$x,new_points$y)
    raster::crs(new_points) <- raster::crs(points)
    if(nrow(new_points) == nrow(points)){
      return(new_points)
    }else{
      points <- new_points
    }
  }
  return(points)
}



#' A simple function to prepare data before the NKDE calculation
#'
#' @param samples A spatialPointsDataFrame of the samples points
#' @param lines A SpatialLinesDataFrame representing the network
#' @param events A spatialPointsDataFrame of the events points
#' @param w A numeric vector reprsenting the weight of the events
#' @param digits the number of digits to keep
#' @param tol a float indicating the spatial tolerance when snapping events on
#' lines
#' @param agg a double indicating if the points must be aggregated within a distance.
#' if NULL, then the points are aggregated by rouding the coordinates.
#' @return the data prepared for the rest of the operations
#' @importFrom data.table tstrsplit setDF
#' @importFrom rgeos gLength
#' @importFrom maptools snapPointsToLines
#' @examples
#' #This is an internal function, no example provided
prepare_data <- function(samples,lines,events, w ,digits,tol, agg){

  ## step1 cleaning the events
  events$weight <- w
  events <- clean_events(events,digits,agg)

  ## step2 defining the global IDS
  events$goid <- 1:nrow(events)
  samples$goid <- 1:nrow(samples)
  samples <- samples[c("goid")]

  ## step3 remove lines with no length
  lines$length <- gLength(lines,byid=TRUE)
  lines <- subset(lines, lines$length>0)
  lines$oid <- 1:nrow(lines)

  return(list("samples" = samples,
              "lines" = lines[c("length","oid")],
              "events" = events))

}



#' Function to split the dataset according to a grid
#'
#' @param grid a spatial grid to split the data within
#' @param samples A spatialPointsDataFrame of the samples points
#' @param events A spatialPointsDataFrame of the events points
#' @param lines A SpatialLinesDataFrame representing the network
#' @param bw the kernel bandwidth (used to avoid egde effect)
#' @param digits the number of digits to keep
#' @param tol a float indicating the spatial tolerance when snapping events on
#' lines
#' @return a list with the splitted dataset
#' @importFrom rgeos gBuffer
#' @examples
#' #This is an internal function, no example provided
split_by_grid <- function(grid,samples,events,lines,bw,tol, digits){

  ## step1 : creating the spatial trees
  tree_samples <- build_quadtree(samples)
  tree_events <- build_quadtree(events)
  tree_lines <- build_quadtree(lines)

  ## step2 : split the datasets

  selections <- lapply(1:length(grid),function(i){
    square <- grid[i,]
    # selecting the samples in the grid
    sel_samples <- spatial_request(square,tree_samples,samples)
    # if there is no sampling points in the rectangle, then return NULL
    if(nrow(sel_samples)==0){
      return(NULL)
    }
    # selecting the events in a buffer
    buff <- gBuffer(square,width=bw)
    buff2 <- gBuffer(square,width=(bw+0.5*bw))
    sel_events <- spatial_request(buff,tree_events,events)
    # selecting the lines in a buffer
    sel_lines <- spatial_request(buff2,tree_lines,lines)
    sel_lines$oid <- 1:nrow(sel_lines)

    # snapping the events on the lines
    if(nrow(sel_events)==0){
      new_lines <- sel_lines
    }else{
      a <- nrow(sel_events)
      b <- nrow(sel_lines)
      x <-  a*b
      if(is.na(x)){
        stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      }
      if(x >= 2*10^9){
        stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      }
      snapped_events <- snapPointsToLines(sel_events,sel_lines,idField = "oid")
      sel_events <- cbind(snapped_events,sel_events)
      new_lines <- add_vertices_lines(sel_lines,sel_events,sel_events$nearest_line_id,tol)
    }

    # split lines at events
    new_lines <- simple_lines(new_lines)
    new_lines$length <- gLength(new_lines,byid = T)
    new_lines <- subset(new_lines,new_lines$length>0)
    new_lines$oid <- 1:nrow(new_lines)
    new_lines <- new_lines[c("length","oid")]


    return(list("samples" = sel_samples,
                "events" = sel_events,
                "lines" = new_lines))
  })
  #let us remove all the empty quadra
  selections <- selections[lengths(selections) != 0]

  for(i in 1:length(selections)){
    selections[[i]]$index <- i
  }

  return(selections)

}



#' Function to split the dataset according to a grid (multicore version)
#'
#' @param grid a spatial grid to split the data within
#' @param samples A spatialPointsDataFrame of the samples points
#' @param events A spatialPointsDataFrame of the events points
#' @param lines A SpatialLinesDataFrame representing the network
#' @param bw the kernel bandwidth (used to avoid egde effect)
#' @param digits the number of digits to keep
#' @param tol a float indicating the spatial tolerance when snapping events on
#' lines
#' @return a list with the splitted dataset
#' @importFrom rgeos gBuffer
#' @examples
#' #This is an internal function, no example provided
split_by_grid.mc <- function(grid,samples,events,lines,bw,tol,digits){

  ## step1 creating the spatial trees
  tree_samples <- build_quadtree(samples)
  tree_events <- build_quadtree(events)
  tree_lines <- build_quadtree(lines)

  ## step2 split the datasets
  #NB : because we can't send c++ pointer to child process, we must start with
  #splitting the spatial objects in a sequential way and only after
  #snapping and splitting in a multicore fashion

  sub_samples <- lapply(1:length(grid),function(i){
    square <- grid[i,]
    #selecting the samples in the grid
    sel_samples <- spatial_request(square,tree_samples,samples)
    ##return NULL if there is no sampling point in the rectangle
    if(nrow(sel_samples)==0){
      return(NULL)
    }
    #selecting the events in a buffer
    buff <- gBuffer(square,width=bw)
    buff2 <- gBuffer(square,width=(bw+0.5*bw))
    sel_events <- spatial_request(buff,tree_events,events)
    #selecting the lines in a buffer
    sel_lines <- spatial_request(buff2,tree_lines,lines)
    sel_lines$oid <- 1:nrow(sel_lines)
    return(list("sel_lines" = sel_lines,
                "sel_events"=sel_events,
                "sel_samples"=sel_samples))
    })

  selections <- future.apply::future_lapply(sub_samples,function(sub){
    if(is.null(sub)){
      return (NULL)
    }
    sel_lines <- sub$sel_lines
    sel_events <- sub$sel_events
    sel_samples <- sub$sel_samples

    #snapping the events on the lines
    if(nrow(sel_events)==0){
      new_lines <- sel_lines
    }else{
      a <- nrow(sel_events)
      b <- nrow(sel_lines)
      x <-  a*b
      if(is.na(x)){
        stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      }
      if(x >= 2*10^9){
        stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      }
      snapped_events <- snapPointsToLines(sel_events,sel_lines,idField = "oid")
      sel_events <- cbind(snapped_events,sel_events)
      invisible(capture.output(new_lines <- add_vertices_lines(sel_lines,sel_events,sel_events$nearest_line_id,tol)))
    }

    #split lines at events
    new_lines <- simple_lines(new_lines)
    new_lines$length <- gLength(new_lines,byid = T)
    new_lines <- subset(new_lines,new_lines$length>0)
    new_lines$oid <- 1:nrow(new_lines)
    new_lines <- new_lines[c("length","oid")]


    return(list("samples" = sel_samples,
                "events" = sel_events,
                "lines" = new_lines))
  })
  #let us remove the empty regions
  selections <- selections[lengths(selections) != 0]

  for(i in 1:length(selections)){
    selections[[i]]$index <- i
  }

  return(selections)

}


#' Function to calculate adaptative bandwidth according to Abramson’s smoothing regimen
#'
#' @param grid a spatial grid to split the data within
#' @param events A spatialPointsDataFrame of the events points
#' @param lines A SpatialLinesDataFrame representing the network
#' @param bw the fixed kernel bandwidth
#' @param trim_bw the maximum size of local bandiwidths
#' @param method The method to use when calculating the NKDE
#' @param kernel_name The name of the kernel to use
#' @param max_depth the maximum recursion depth
#' @param digits the number of digits to keep
#' @param tol a float indicating the spatial tolerance when snapping events on
#' lines
#' @param sparse a boolean indicating if sparse matrix should be used
#' @param verbose a boolean indicating if update messages should be printed
#' @return a vector with the local bandwidth
#' @examples
#' #This is an internal function, no example provided
adaptive_bw <- function(grid,events,lines,bw,trim_bw,method,kernel_name,max_depth,tol,digits,sparse,verbose){
  ##step 1 split the datas !
  selections <- split_by_grid(grid,events,events,lines,trim_bw, tol, digits)
  ## step 2 calculating the temp NKDE values
  if(verbose){
    print("start calculating the kernel values for the adaptive bandwidth...")
  }
  n_quadra <- length(selections)
  dfs <- lapply(1:n_quadra,function(i){
    sel <- selections[[i]]
    bws <- rep(bw,nrow(sel$events))
    if(verbose){
      print(paste("    quadra ",i,"/",n_quadra,sep=""))
    }
    values <- nkde_worker(sel$lines, sel$events,
                          sel$samples, kernel_name,bw,
                          bws, method, div = "none", digits,
                          tol,sparse, max_depth, verbose)
    if(any(is.na(values))){
      print(paste("NA values here : ",i,sep=""))
    }
    df <- data.frame("goid"=sel$samples$goid,
                     "k" = values)
    return(df)
  })

  ## step 3  combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]

  ## step 4 calculating the new bandwidth !
  delta <- calc_gamma(tot_df$k)
  new_bw <- bw * (tot_df$k**(-1/2) * delta**(-1))
  return(new_bw)
}


#' Function to calculate adaptative bandwidth according to Abramson's smoothing regimen (multicore version)
#'
#' @param grid a spatial grid to split the data within
#' @param events A spatialPointsDataFrame of the events points
#' @param lines A SpatialLinesDataFrame representing the network
#' @param bw the fixed kernel bandwidth
#' @param trim_bw the maximum size of local bandiwidths
#' @param method The method to use when calculating the NKDE
#' @param kernel_name The name of the kernel to use
#' @param max_depth the maximum recursion depth
#' @param digits the number of digits to keep
#' @param tol a float indicating the spatial tolerance when snapping events on
#' lines
#' @param sparse a boolean indicating if sparse matrix should be used
#' @param verbose a boolean indicating if update messages should be printed
#' @return a vector with the local bandwidth
#' @examples
#' #This is an internal function, no example provided
adaptive_bw.mc <- function(grid,events,lines,bw,trim_bw,method,kernel_name,max_depth,tol,digits,sparse,verbose){
  ##step 1 split the datas !
  selections <- split_by_grid.mc(grid,events,events,lines,trim_bw, tol, digits)
  ## step 2 calculating the temp NKDE values
  if(verbose){
    print("start calculating the kernel values for the adaptive bandwidth...")
  }

  if(verbose){
    print("starting the kernel calculations ...")
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {
        bws <- rep(bw,nrow(sel$events))
        invisible(capture.output(values <- nkde_worker(sel$lines, sel$events,
                                                       sel$samples, kernel_name,bw,
                                                       bws, method, div = "none", digits,
                                                       tol,sparse, max_depth, verbose)))

        df <- data.frame("goid"=sel$samples$goid,
                         "k" = values)
        p(sprintf("i=%g", sel$index))
        return(df)
      })
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {
      bws <- rep(bw,nrow(sel$events))
      values <- nkde_worker(sel$lines, sel$events,
                                  sel$samples, kernel_name,bw,
                                  bws, method, div = "none", digits,
                                  tol,sparse, max_depth, verbose)

      df <- data.frame("goid"=sel$samples$goid,
                       "k" = values)
      return(df)
    })
  }

  ## step 3  combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]

  ## step 4 calculating the new bandwidth !
  delta <- calc_gamma(tot_df$k)
  new_bw <- bw * (tot_df$k**(-1/2) * delta**(-1))
  return(new_bw)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### worker functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' The worker function for nkde and nkde.mc
#'
#' @param lines A SpatialLinesDataFrame with the sampling points. The
#' geoemtries must be a SpatialLinesDataFrame (may crash if some geometries
#'  are invalid)
#' @param events A SpatialPointsDataFrame representing the events on the
#' network. The points will be snapped on the network.
#' @param samples A SpatialPointsDataFrame representing the locations for
#' which the densities will be estimated.
#' @param kernel_name The name of the kernel to use
#' @param bw the global kernel bandwidth
#' @param bws The kernel bandwidth (in meters) for each event
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#' events within the radius around each sampling point), "bw" (the bandwith)
#' "none" (the simple sum).
#' @param digits The number of digits to keep in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines extremities. When
#' points are closer, they are added at the extermity of the lines.
#' @param sparse a boolean indicating if sparse or regular matrice should be
#' used by the Rcpp functions. Regular matrices are faster, but require more
#' memory and could lead to error, in particular with multiprocessing. Sparse
#' matrices are slower, but require much less memory.
#' @param max_depth when using the continuous and discontinuous methods, the
#' calculation time and memory use can go wild  if the network has a lot of
#' small edges (area with a lot of intersections and a lot of events). To
#' avoid it, it is possible to set here a maximum depth. Considering that the
#' kernel is divided at intersections, a value of 8 should yield good
#' estimates. A larger value can be used without problem for the discontinuous
#' method. For the continuous method, a larger value will strongly impact
#' calculation speed.
#' @param verbose A boolean, indicating if the function should print messages
#' about process.
#' @importFrom igraph adjacent_vertices get.edge.ids
#' @return A numerci vector with the nkde values
#' @examples
#' #This is an internal function, no example provided
nkde_worker <- function(lines, events, samples, kernel_name,bw, bws, method, div, digits, tol, sparse, max_depth, verbose = FALSE){

  # if we do not have event in that space, just return 0 values
  if(nrow(events)==0){
    values <- rep(0,nrow(samples))
    return(values)
  }

  ## step1 creating the graph
  if(verbose){
    print("    build graph ...")
  }
  graph_result <- build_graph(lines,digits = digits,line_weight = "length")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  edges <- graph_result$spedges

  ## step2 for each sample, find its belonging line
  a <- nrow(events)
  b <- nrow(lines)
  x <-  a*b
  if(is.na(x)){
    stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
  }
  if(x >= 2*10^9){
    stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
  }
  snapped_samples <- maptools::snapPointsToLines(samples,edges,idField = "edge_id")
  samples$edge_id <- snapped_samples$nearest_line_id

  ## step3 finding for each event, its corresponding node
  events$vertex_id <- closest_points(events, nodes)

  ## step4 adding the spatial coordinates to samples and nodes
  XY_nodes <- sp::coordinates(nodes)
  nodes$X_coords <- XY_nodes[,1]
  nodes$Y_coords <- XY_nodes[,2]

  XY_samples <- sp::coordinates(samples)
  samples$X_coords <- XY_samples[,1]
  samples$Y_coords <- XY_samples[,2]

  ## step5 adding a local oid for samples and events
  events$oid <- 1:nrow(events)
  samples$oid <- 1:nrow(samples)

  ## step6 starting the calculations !

  if(verbose){
    print("        calculating NKDE values ...")
  }

  if(method == "simple"){
    # the cas of the simple method, no c++ here
    kernel_func <- select_kernel(kernel_name)
    if(verbose){
      values <- simple_nkde(graph, events, samples, bws, kernel_func, nodes, edges)
    }else{
      invisible(capture.output(values <- simple_nkde(graph, events, samples, bws, kernel_func, nodes, edges)))
    }

  }else{
    # we have to call the package with the Rcpp functions
    # this is necessary because this function can be used in a multicore context

    ## step6.5 preparing the complementary object for the rcpp function
    neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")
    neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})

    if(method=="continuous"){
      ##and finally calculating the values
      if (sparse){
        values <- spNetworkCpp::continuous_nkde_cpp_arma_sparse(neighbour_list, events$vertex_id, events$weight,
                                                         samples@data, bws, kernel_name, nodes@data, graph_result$linelist, max_depth, verbose)
      }else{
        values <- spNetworkCpp::continuous_nkde_cpp_arma(neighbour_list, events$vertex_id, events$weight,
                                                                samples@data, bws, kernel_name, nodes@data, graph_result$linelist, max_depth, verbose)
      }

    }

    if(method == "discontinuous"){
        #let this commented here to debug and test sessions
        # invisible(capture.output(values <- discontinuous_nkde2(edge_list,neighbour_list, events$vertex_id, events$weight,
        #                                                       samples@data, bw, kernel_func, nodes@data, graph_result$linelist, max_depth, verbose)))
      if(sparse){
        values <- spNetworkCpp::discontinuous_nkde_cpp_arma_sparse(neighbour_list, events$vertex_id, events$weight,
                                                            samples@data, bws, kernel_name, nodes@data, graph_result$linelist, max_depth, verbose)
      }else{
        values <- spNetworkCpp::discontinuous_nkde_cpp_arma(neighbour_list, events$vertex_id, events$weight,
                                                           samples@data, bws, kernel_name, nodes@data, graph_result$linelist, max_depth, verbose)
      }


    }

  }

  ## step7 adjusting the kernel values !
  if(div == "n"){
    return(values$sum_k / values$n)
  }else if (div == "bw"){
    return(values$sum_k * (1/bw))
  }else if (div == "none"){
    return(values$sum_k)
  }
}




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### main functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Network Kernel density estimate
#'
#' Calculate the Network Kernel Density Estimate based on a network of lines,
#' sampling points, and events
#'
#' **The three NKDE methods**\cr
#' Estimating the density of a point process is commonly done by using an
#' ordinary two dimensional kernel density function. However, there is
#' numerous cases for which the events do not occur in a two dimensional
#' space but on a network (like car crash, outdoor crimes, leaks in pipelines,
#' etc.). New methods were developped to adapt the methodology to networks,
#' three of them are available in this package.
#' \itemize{
#'   \item{method="simple"}{this first method was presented by Xie et al.
#'   (2008) and proposes an intuitive solution. The distances between events
#'   and sampling points are replaced by network distances, and the formula of
#'   the kernel is adapted to calculate the density over a linear unit
#'   instead of an areal unit.}
#'   \item{method="discontinuous"}{the previous method has been critized by
#'   Okabe et al (2008), arguing that the estimator proposed is biaised,
#'   conducting to overestimation of density in hot-spots of events. More
#'   specifically, the simple method does not conserve mass and the induced
#'   kernel is not a probability density along the network. They thus
#'   proposed a discontinuous version of the kernel function on network, which
#'   "divide" equaly the mass density of an event at intersections}
#'   \item{method="continuous"}{if the discontinuous method is unbiased, it
#'   leads to a discontinuous kernel function which is a bit counter-intuitive.
#'   Okabe et al (2008) proposed another version of the kernel, that divide
#'   the mass of the density at intersection but adjust the density before the
#'   intersection to keep the function continuous.}
#' }
#' The three methods are available because, despite the fact that the simple
#' method is less exact statistically speaking, it might be more intuitive.
#' In a purely geographical view, it migh be seen as sort of distance decay
#' function like used in Geographically Weighted Regression.\cr
#' \cr\cr
#' **adaptive bandwidth**\cr
#' It is possible to use adaptive bandiwdth instead of fixed bandwidth. The
#' adaptive bandwidth are calculated using the Abramson’s smoothing regimen.
#' To do so, a orignal fixed bandiwdth must be specified (bw parameter), and
#' is used to estimate a priory densities at event locations. These densities
#' are then used to calculate local bandwidth. The maximum size of the local
#' bandwidth can be limited with the parameter trim_bw. For more details, look
#' at the vignette a_NKDE.
#' \cr\cr
#' **Optimization parameters**\cr
#' The grid_shape parameter allows to split the calculus of the NKDE according
#' to a grid dividing the study area. It might be necessary for big dataset
#' to reduce the memory used. If the grid_shape is c(1,1), then a full network
#' is build for the area. If the grid_shape is c(2,2), then the area is
#' split in 4 rectangles. For each rectangle, the sample points falling in the
#' rectangle are used, the events and the lines in a radius of the bandwidth
#' length are used The results are combined at the end and ordered to match
#' the original order of the samples.
#' \cr\cr
#' The geographical coordinates of the start and end of lines are used to build
#' the network. To avoid troubles with digits, we truncate the coordinates
#' according to the digit parameter. A minimal loss of precision is expected
#' but results in a fast construction of the network.
#' \cr\cr
#' To calculate the distances on the network, all the events are added as
#' vertices. To reduce the size of the network, it is possible to reduce the
#' number of vertices by adding the events at the extremity of
#' the lines if they are close to them. This is controled by the parameter tol.
#' \cr\cr
#' In the same way, it is possible to limit the number of vertices by
#' aggregating the events that are close to each other. In that case, the
#' weights of the aggregated events are summed. According to an aggregation
#' distance, a buffer is drawn around the fist event, each other event falling
#' in that buffer are aggregated to the first event, forming a new event. The
#' coordinates of this new event are the mean of the original events
#' coordinates. This procedure is repeated until no events are aggregated. The
#' aggregation distance can be fixed with the parameter agg.
#' \cr\cr
#' When using the continuous and discontinuous kernel, the density is reduced
#' at each intersection crossed. After 3 intersections with four directions
#' each, the density is divided by 27 (3x3x3), leading to very small values.
#' To reduce calculation time with a small precision loss, it is recommanded
#' to set a maximum depth value for the two methods. This is controled by the
#' max_depth parameter.
#' \cr\cr
#' When using the continuous and discontinuous kernel, the connexions between
#' graph nodes are stored in a matrix. This matrix is typically sparse, and
#' so a sparse matrix object is used to limit memory use. If the network is
#' small (typically when the grid used to split the data has small rectangles)
#' then a classical matrix could be used instead of a sparse one. It increases
#' singnificantly speed, but could lead to memory issues.
#'
#' @param lines A SpatialLinesDataFrame with the sampling points. The
#' geoemtries must be a SpatialLinesDataFrame (may crash if some geometries
#'  are invalid)
#' @param events A SpatialPointsDataFrame representing the events on the
#' network. The points will be snapped on the network.
#' @param w A vector representing the weight of each event
#' @param samples A SpatialPointsDataFrame representing the locations for
#' which the densities will be estimated.
#' @param kernel_name The name of the kernel to use. Must be one of triangle,
#' gaussian, tricube, cosine ,triweight, quartic, or epanechnikov.
#' @param bw The kernel bandwidth (in meters)
#' @param adaptive A boolean, indicating if an adaptive bandwidth must be
#' used
#' @param trim_bw A float, indicating the maximum value for the adaptive
#' bandwidth
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#' events within the radius around each sampling point), "bw" (the bandwith)
#' "none" (the simple sum).
#' @param diggle_correction A boolean indicating if the correction factor
#' for edge effect must be used.
#' @param study_area A SpatialPolygonsDataFrame or a SpatialPolygon
#' representing the limits of the study area.
#' @param max_depth when using the continuous and discontinuous methods, the
#' calculation time and memory use can go wild  if the network has a lot of
#' small edges (area with a lot of intersections and a lot of events). To
#' avoid it, it is possible to set here a maximum depth. Considering that the
#' kernel is divided at intersections, a value of 10 should yield good
#' estimates in most cases. A larger value can be used without problem for the
#' discontinuous method. For the continuous method, a larger value will
#' strongly impact calculation speed.
#' @param digits The number of digits to keep in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines extremities. When
#' points are closer, they are added at the extermity of the lines.
#' @param agg a double indicating if the events must be aggregated within a distance.
#' if NULL, then the events are aggregated by rounding the coordinates.
#' @param sparse a boolean indicating if sparse or regular matrice should be
#' used by the Rcpp functions. Regular matrices are faster, but require more
#' memory and could lead to error, in particular with multiprocessing. Sparse
#' matrices are slower, but require much less memory.
#' @param grid_shape A vector of two values indicating how the study area
#' must be splitted when performing the calculus (see details). Defaut is c(1,1)
#' @param verbose A boolean, indicating if the function should print messages
#' about process.
#' @param check A boolean indicating if the geometry checks must be run before
#' calculating the densities
#' @return A vector of values, they are the density estimates at samplings
#' points
#' @export
#' @examples
#' data(mtl_network)
#' data(bike_accidents)
#' lixels <- lixelize_lines(mtl_network,200,mindist = 50)
#' samples <- lines_center(lixels)
#' densities <- nkde(mtl_network,
#'                   events = bike_accidents,
#'                   w = rep(1,nrow(bike_accidents)),
#'                   samples = samples,
#'                   kernel_name = "quartic",
#'                   bw = 300, div= "bw",
#'                   adaptive = FALSE,
#'                   method = "discontinuous", digits = 1, tol = 1,
#'                   agg = 15,
#'                   grid_shape = c(1,1),
#'                   verbose=FALSE)
nkde <- function(lines, events, w, samples, kernel_name, bw, adaptive=FALSE, trim_bw=NULL, method, div="bw",diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE, grid_shape=c(1,1), verbose=TRUE, check=TRUE){

  ## step0 basic checks
  if(verbose){
    print("checking inputs ...")
  }
  if((kernel_name %in% c("triangle", "gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov'))==FALSE){
    stop('the name of the kernel function must be one of c("triangle", "gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov")')
  }

  if(method %in% c("simple","continuous","discontinuous") == FALSE){
    stop('the method must be one of c("simple","continuous","discontinuous"')
  }

  if(bw<=0){
    stop("the bandwidth for the kernel must be superior to 0")
  }

  if(adaptive & is.null(trim_bw)){
    stop("if adaptive is TRUE, a value for trim_bw must be supplied")
  }
  if(diggle_correction & is.null(study_area)){
    stop("the study_area must be defined if the Diggle correction factor is used")
  }
  if(check){
    check_geometries(lines,samples,events, study_area)
  }


  ## step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
  }
  data <- prepare_data(samples, lines, events,w,digits,tol,agg)
  lines <- data$lines
  samples <- data$samples
  events <- data$events

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))

  ## adaptive bandwidth !
  if(adaptive==FALSE){
    bws <- rep(bw,nrow(events))
  }else{
    ## we want to use an adaptive bw
    bws <- adaptive_bw(grid,events,lines,bw,trim_bw,method,kernel_name,max_depth,tol,digits,sparse,verbose)
  }

  ## calculating the correction factor
  if(diggle_correction){
    if(verbose){
      print("Calculating the correction factor")
    }
    corr_factor <- correction_factor(study_area,events,lines,method,bws, kernel_name, tol, digits, max_depth, sparse)
  }else{
    corr_factor <- rep(1,nrow(events))
  }
  events$weight <- events$weight * corr_factor

  events$bw <- bws
  max_bw <- max(bws)
  ## step3 splitting the dataset with each rectangle
  selections <- split_by_grid(grid,samples,events,lines,max_bw, tol, digits)

  ## step 4 calculating the values
  if(verbose){
    print("start calculating the kernel values ...")
  }
  n_quadra <- length(selections)
  dfs <- lapply(1:n_quadra,function(i){
    sel <- selections[[i]]

    if(verbose){
      print(paste("    quadra ",i,"/",n_quadra,sep=""))
    }

    values <- nkde_worker(sel$lines, sel$events,
                          sel$samples, kernel_name,bw,
                          sel$events$bw, method, div, digits,
                          tol,sparse, max_depth, verbose)

    df <- data.frame("goid"=sel$samples$goid,
                     "k" = values)
    return(df)
  })

  if(verbose){
    print("combining the results ...")
  }

  ## step5  combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]
  if(adaptive){
    return(list("events" = events,
                "k" = tot_df$k))
  }else{
    return(tot_df$k)
  }
}


#' Network Kernel density estimate
#'
#' Calculate the Network Kernel Density Estimate based on a network of lines,
#' sampling points, and events. This version can use the current defined plan
#' with the package future. This slightly increase the calculus to do, but can
#' increase speed a lot by dividing the burden on multiple cores.
#'
#' for details, please see the function nkde
#'
#' @param lines A SpatialLinesDataFrame with the sampling points. The
#' geoemtries must be a SpatialLinesDataFrame (may crash if some geometries
#'  are invalid)
#' @param events A SpatialPointsDataFrame representing the events on the
#' network. The points will be snapped on the network.
#' @param w A vector representing the weight of each event
#' @param samples A SpatialPointsDataFrame representing the locations for
#' which the densities will be estimated
#' @param kernel_name The name of the kernel to use. Must be one of triangle,
#' gaussian, tricube, cosine ,triweight, quartic, or epanechnikov.
#' @param bw The kernel bandwidth (in meters)
#' @param adaptive A boolean, indicating if an adaptive bandwidth must be
#' used
#' @param trim_bw A float, indicating the maximum value for the adaptive
#' bandwidth
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#' events within the radius around each sampling point), "bw" (the bandwith)
#' "none" (the simple sum).
#' @param diggle_correction A boolean indicating if the correction factor
#' for edge effect must be used.
#' @param study_area A SpatialPolygonsDataFrame or a SpatialPolygon
#' representing the limits of the study area.
#' @param max_depth when using the continuous and discontinuous methods, the
#' calculation time and memory use can go wild  if the network has a lot of
#' small edges (area with a lot of intersections and a lot of events). To
#' avoid it, it is possible to set here a maximum depth. Considering that the
#' kernel is divided at intersections, a value of 10 should yield good
#' estimates in most cases. A larger value can be used without problem for the
#' discontinuous method. For the continuous method, a larger value will
#' strongly impact calculation speed.
#' @param digits The number of digits to keep in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines extremities. When
#' points are closer, they are added at the extermity of the lines.
#' @param agg a double indicating if the events must be aggregated within a distance.
#' if NULL, then the events are aggregated by rounding the coordinates.
#' @param sparse a boolean indicating if sparse or regular matrice should be
#' used by the Rcpp functions. Regular matrices are faster, but require more
#' memory and could lead to error, in particular with multiprocessing. Sparse
#' matrices are slower, but require much less memory.
#' @param grid_shape A vector of two values indicating how the study area
#' must be splitted when performing the calculus (see details). Defaut is c(1,1)
#' @param verbose A boolean, indicating if the function should print messages
#' about process.
#' @param check A boolean indicating if the geometry checks must be run before
#' calculating the densities
#' @return A vector of values, they are the density estimates at samplings
#' points
#' @export
#' @examples
#' data(mtl_network)
#' data(bike_accidents)
#' future::plan(future::multiprocess(workers=2))
#' lixels <- lixelize_lines(mtl_network,200,mindist = 50)
#' samples <- lines_center(lixels)
#' densities <- nkde.mc(mtl_network,
#'                   events = bike_accidents,
#'                   w = rep(1,nrow(bike_accidents)),
#'                   samples = samples,
#'                   kernel_name = "quartic",
#'                   bw = 300, div= "bw",
#'                   adaptive = FALSE, agg = 15,
#'                   method = "discontinuous", digits = 1, tol = 1,
#'                   grid_shape = c(3,3),
#'                   verbose=FALSE)
#' \dontshow{
#'    ## R CMD check: make sure any open connections are closed afterward
#'    if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'    }
nkde.mc <- function(lines, events, w, samples, kernel_name, bw, adaptive=FALSE, trim_bw=NULL, method, div="bw", diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1,agg=NULL, sparse=TRUE, grid_shape=c(1,1), verbose=TRUE, check=TRUE){

  ## step0 basic checks
  if(verbose){
    print("checking inputs ...")
  }

  if((kernel_name %in% c("triangle", "gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov'))==FALSE){
    stop('the name of the kernel function must be one of c("triangle", "gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov")')
  }

  if(method %in% c("simple","continuous","discontinuous") == FALSE){
    stop('The method must be one of c("simple","continuous","discontinuous"')
  }

  if(bw<=0){
    stop("the bandwidth for the kernel must be superior to 0")
  }
  if(adaptive & is.null(trim_bw)){
    stop("if adaptive is TRUE, a value for trim_bw must be supplied")
  }
  if(diggle_correction & is.null(study_area)){
    stop("the study_area must be defined if the Diggle correction factor is used")
  }

  if(check){
    check_geometries(lines,samples,events,study_area)
  }

  ## step1 preparing the data
  if(verbose){
    print("prior data preparation ...")
  }
  data <- prepare_data(samples, lines, events,w,digits,tol,agg)
  lines <- data$lines
  samples <- data$samples
  events <- data$events

  ## step2 creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))

  ## adaptive bandwidth !
  if(adaptive==FALSE){
    bws <- rep(bw,nrow(events))
  }else{
    if(verbose){
      print("calculating the local bandwidth ...")
    }
    ## we want to use an adaptive bw
    bws <- adaptive_bw.mc(grid,events,lines,bw,trim_bw,method,kernel_name,max_depth,tol,digits,sparse,verbose)
  }

  ## calculating the correction factor
  if(diggle_correction){
    if(verbose){
      print("Calculating the correction factor")
    }
    corr_factor <- correction_factor(study_area,events,lines,method,bws, kernel_name, tol, digits, max_depth, sparse)
  }else{
    corr_factor <- rep(1,nrow(events))
  }
  events$weight <- events$weight * corr_factor

  events$bw <- bws
  max_bw <- max(bws)

  ## step3 splitting the dataset with each rectangle
  selections <- split_by_grid.mc(grid,samples,events,lines,max_bw, digits,tol)

  ## step4 calculating the values

  if(verbose){
    print("starting the kernel calculations ...")
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {

        invisible(capture.output(values <- nkde_worker(sel$lines, sel$events,
                              sel$samples, kernel_name,bw,
                              sel$events$bw, method, div, digits,
                              tol,sparse, max_depth, verbose)))

        df <- data.frame("goid"=sel$samples$goid,
                         "k" = values)
        p(sprintf("i=%g", sel$index))
        return(df)
      })
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {

      values <- nkde_worker(sel$lines, sel$events,
                            sel$samples, kernel_name,bw,
                            sel$events$bw, method, div, digits,
                            tol,sparse, max_depth, verbose)

      df <- data.frame("goid"=sel$samples$goid,
                       "k" = values)
      return(df)
    })
  }


  if(verbose){
    print("combining the results ...")
  }
  ## step5 combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]
  if(adaptive){
    return(list("events" = events,
                "k" = tot_df$k))
  }else{
    return(tot_df$k)
  }
}


