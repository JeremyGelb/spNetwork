# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### helper functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function to avoid having events at the same location
#'
#' @param events the SpatialPointsDataFrame to contract (must have a weight column)
#' @param digits the number of digits to keep
#' @return a new SpatialPointsDataFrame
#' @importFrom data.table tstrsplit setDF
#' @examples
#' #This is an internal function, no example provided
clean_events <- function(events,digits=5){
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
#' @return the data prepared for the rest of the operations
#' @importFrom data.table tstrsplit setDF
#' @importFrom rgeos gLength
#' @importFrom maptools snapPointsToLines
#' @examples
#' #This is an internal function, no example provided
prepare_data <- function(samples,lines,events, w ,digits,tol){

  ##step1 : cleaning the events
  events$weight <- w
  events <- clean_events(events,digits)

  ##step2 : defining the global IDS
  events$goid <- 1:nrow(events)
  samples$goid <- 1:nrow(samples)
  samples <- samples[c("goid")]

  ##step3 : remove lines with no length
  lines$length <- gLength(lines,byid=TRUE)
  lines <- subset(lines, lines$length>0)

  ##step4 : snapp events on lines
  lines$oid <- 1:nrow(lines)
  snapped_events <- snapPointsToLines(events,lines,idField = "oid")
  events <- cbind(snapped_events,events)

  ##step5 : split lines at events
  new_lines <- add_vertices_lines(lines,events,events$nearest_line_id,tol)
  new_lines <- simple_lines(new_lines)
  new_lines$length <- gLength(new_lines,byid = T)
  new_lines <- subset(new_lines,new_lines$length>0)
  new_lines$oid <- 1:nrow(new_lines)
  new_lines <- new_lines[c("length","oid")]

  return(list("samples" = samples,
              "lines" = new_lines,
              "events" = events))

}



#' Function to split the dataset according to a grid
#'
#' @param grid a spatial grid to split the data within
#' @param samples A spatialPointsDataFrame of the samples points
#' @param events A spatialPointsDataFrame of the events points
#' @param lines A SpatialLinesDataFrame representing the network
#' @param bw the kernel bandwidth (used to avoid egde effect)
#' @return a list with the splitted dataset
#' @importFrom rgeos gBuffer
#' @examples
#' #This is an internal function, no example provided
split_by_grid <- function(grid,samples,events,lines,bw){

  ##step1 : creating the spatial trees
  tree_samples <- build_quadtree(samples)
  tree_events <- build_quadtree(events)
  tree_lines <- build_quadtree(lines)

  ##step2 : split the datasets

  selections <- lapply(1:length(grid),function(i){
    square <- grid[i,]
    #selecting the samples in the grid
    sel_samples <- spatial_request(square,tree_samples,samples)
    ##si on a aucun point d'echantillonnage dans la zone, c'est quelle est vide
    if(nrow(sel_samples)==0){
      return(NULL)
    }
    #selecting the events in a buffer
    buff <- gBuffer(square,width=bw)
    sel_events <- spatial_request(buff,tree_events,events)
    #selecting the lines in a buffer
    sel_lines <- spatial_request(buff,tree_lines,lines)
    return(list("samples" = sel_samples,
                "events" = sel_events,
                "lines" = sel_lines))
  })
  #enlevons les quadra vides
  selections <- selections[lengths(selections) != 0]

  for(i in 1:length(selections)){
    selections[[i]]$index <- i
  }

  return(selections)

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
#' @param kernel_func A kernel function, got from select_kernel
#' @param bw The kernel bandwidth (in meters)
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
#' @param max_depth when using the continuous method, a recursive function
#' is applied. When links in network are very short, this might lead to very
#' deep recursion and crash. To avoid it, it is possible to set here a maximym
#' recusion depth. Okabe recommand a value of 3.
#' approximation.
#' @param verbose A boolean, indicating if the function should print messages
#' about process.
#' @importFrom igraph adjacent_vertices get.edge.ids
#' @return A numerci vector with the nkde values
#' @examples
#' #This is an internal function, no example provided
nkde_worker <- function(lines, events, samples, kernel_func, bw, method, div, digits, tol, max_depth, verbose = FALSE){

  ##si on a pas d'evenement, on renvoit que des 0
  if(nrow(events)==0){
    values <- rep(0,nrow(samples))
    return(values)
  }


  ##step1 : creating the graph
  if(verbose){
    "        build graph ..."
  }
  graph_result <- build_graph(lines,digits = digits,line_weight = "length")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  edges <- graph_result$spedges

  ##step2 : for each sample, find its belonging line
  snapped_samples <- maptools::snapPointsToLines(samples,edges,idField = "edge_id")
  samples$edge_id <- snapped_samples$nearest_line_id


  ##step3 : finding for each event, its node
  events$vertex_id <- closest_points(events, nodes)

  ##step4 : adding the spatial coordinates to samples and nodes
  XY_nodes <- sp::coordinates(nodes)
  nodes$X_coords <- XY_nodes[,1]
  nodes$Y_coords <- XY_nodes[,2]

  XY_samples <- sp::coordinates(samples)
  samples$X_coords <- XY_samples[,1]
  samples$Y_coords <- XY_samples[,2]

  ##step5 : adding a local oid for samples and events
  events$oid <- 1:nrow(events)
  samples$oid <- 1:nrow(samples)

  ##step6 : starting the calculations !

  if(verbose){
    "        calculating NKDE values ..."
  }

  if(method == "simple"){
    if(verbose){
      values <- simple_nkde(graph, events, samples, bw, kernel_func, nodes, edges)
    }else{
      invisible(capture.output(values <- simple_nkde(graph, events, samples, bw, kernel_func, nodes, edges)))
    }

  }
  if(method == "discontinuous"){
    if(verbose){
      values <- discontinuous_nkde(graph, events, samples, bw, kernel_func, nodes, edges)
    }else{
      invisible(capture.output(values <- discontinuous_nkde(graph, events, samples, bw, kernel_func, nodes, edges)))
    }

  }
  if(method=="continuous"){
    ## we have to call the package with the Rcpp functions
    ## this is necessary because this function can be used in a multicore context
    #library(spNetworkCpp)
    ### preparing the complementary object for the rcpp function
    # the neighbours list
    neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")
    neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})
    # the edge list
    lists <- lapply(1:length(neighbour_list),function(i){
      n1 <- i
      neighbours <- neighbour_list[[n1]]
      eids <- cbind(rep(n1,length(neighbours)),neighbours)
      eids <- c(t(eids))
      edges_id <- as.numeric(get.edge.ids(graph,eids))
      node_names1 <- paste(n1,neighbours,sep="_")
      l1 <- as.list(edges_id)
      names(l1) <- node_names1
      return(l1)
    })
    edge_list <- unlist(lists,recursive = F)

    ##and finally calculating the values
    values <- spNetworkCpp::continuous_nkde_cpp(edge_list,neighbour_list, events$vertex_id, events$weight,
                                                  samples@data, bw, kernel_func, nodes@data, graph_result$linelist, max_depth, verbose)

  }

  ##step7 : adjusting the kernel values !


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
#' In a geographical mind, it migh be seen as sort of distance decay
#' function like used in Geographically Weighted Regression.\cr
#' \cr
#' The grid_shape parameter allows to split the calculus of the NKDE according
#' to a grid dividing the study area. It might be necessary for big dataset
#' to reduce the memory used. If the grid_shape is c(1,1), then a full network
#' is build for the area. If the grid_shape is c(2,2), then the area is
#' split in 4 rectangles. For each rectangle, the sample points falling in the
#' rectangle are used, the events in a radius of the bandwidth length are used
#' and the lines in a radius of the bandwidth length are used. The results are
#' combined at the end and ordered to match the original order of the samples.
#' \cr\cr
#' The geographical coordinates of the start and end nodes are used to build
#' the network. To avoid troubles with digits, we truncate the coordinates
#' according to the digit parameter. A minimal loss of precision is expected
#' but results in a fast construction of the network.\cr\cr
#' To calculate the distances on the network, all the sampling points and the
#' events are added as vertices. To reduce the size of the network, it is
#' possible to reduce the number of vertex in the network by adding the events
#' and the sampling points at the extremity of the lines if they are close to
#' them. This is controled by the parameter tol.
#'
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
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#' events within the radius around each sampling point), "bw" (the bandwith)
#' "none" (the simple sum).
#' @param max_depth when using the continuous method, a recursive function
#' is applied. When links in network are very short, this might lead to very
#' deep recursion and crash. To avoid it, it is possible to set here a maximym
#' recusion depth. Okabe recommand a value of 3, but it highly depends on the
#' length of the lines of the network. 15 should be a good approximation in
#' most cases, but might be time consuming.
#' @param digits The number of digits to keep in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines extremities. When
#' points are closer, they are added at the extermity of the lines.
#' @param grid_shape A vector of two values indicating how the study area
#' must be splitted when performing the calculus (see details). Defaut is c(1,1)
#' @param verbose A boolean, indicating if the function should print messages
#' about process.
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
#'                   method = "discontinuous", digits = 1, tol = 1,
#'                   grid_shape = c(1,1),
#'                   verbose=FALSE)
nkde <- function(lines, events, w, samples, kernel_name, bw, method, div="bw", max_depth = 15, digits=5, tol=0.1, grid_shape=c(1,1), verbose=TRUE){

  ##step0
  kernel_func <- select_kernel(kernel_name)

  ##step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
    data <- prepare_data(samples, lines, events,w,digits,tol)
  }else{
    invisible(capture.output(data <- prepare_data(samples, lines, events,w,digits,tol)))
  }

  lines <- data$lines
  samples <- data$samples
  events <- data$events

  ##step2 : creating the grid
  grid <- build_grid(grid_shape,lines)

  ##step3 : splitting the dataset with each rectangle
  selections <- split_by_grid(grid,samples,events,lines,bw)

  ##step4 : calculating the values
  n_quadra <- length(selections)
  dfs <- lapply(1:n_quadra,function(i){

    sel <- selections[[i]]

    if(verbose){
      print(paste("    quadra ",i,"/",n_quadra,sep=""))
    }

    values <- nkde_worker(sel$lines, sel$events,
                          sel$samples, kernel_func,
                          bw, method, div, digits,
                          tol, max_depth, verbose)

    df <- data.frame("goid"=sel$samples$goid,
                     "k" = values)
    return(df)
  })

  if(verbose){
    print("combining the results ...")
  }

  ##step5 : combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]
  return(tot_df$k)

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
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#' events within the radius around each sampling point), "bw" (the bandwith)
#' "none" (the simple sum).
#' @param max_depth when using the continuous method, a recursive function
#' is applied. When links in network are very short, this might lead to very
#' deep recursion and crash. To avoid it, it is possible to set here a maximym
#' recusion depth. Okabe recommand a value of 3, but it highly depends on the
#' length of the lines of the network. 15 should be a good approximation in
#' most cases, but might be time consuming.
#' @param digits The number of digits to keep in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines extremities. When
#' points are closer, they are added at the extermity of the lines.
#' @param grid_shape A vector of two values indicating how the study area
#' must be splitted when performing the calculus (see details). Defaut is c(1,1)
#' @param verbose A boolean, indicating if the function should print messages
#' about process.
#' @return A vector of values, they are the density estimates at samplings
#' points
#' @export
#' @examples
#' data(mtl_network)
#' data(bike_accidents)
#' future::plan(future::multiprocess(workers=4))
#' lixels <- lixelize_lines(mtl_network,200,mindist = 50)
#' samples <- lines_center(lixels)
#' densities <- nkde.mc(mtl_network,
#'                   events = bike_accidents,
#'                   w = rep(1,nrow(bike_accidents)),
#'                   samples = samples,
#'                   kernel_name = "quartic",
#'                   bw = 300, div= "bw",
#'                   method = "discontinuous", digits = 1, tol = 1,
#'                   grid_shape = c(3,3),
#'                   verbose=FALSE)
#' \dontshow{
#'    ## R CMD check: make sure any open connections are closed afterward
#'    if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
nkde.mc <- function(lines, events, w, samples, kernel_name, bw, method, div="bw", max_depth = 15, digits=5, tol=0.1, grid_shape=c(1,1), verbose=TRUE){

  ##step0
  kernel_func <- select_kernel(kernel_name)

  ##step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
    data <- prepare_data(samples, lines, events,w,digits,tol)
  }else{
    invisible(capture.output(data <- prepare_data(samples, lines, events,w,digits,tol)))
  }
  lines <- data$lines
  samples <- data$samples
  events <- data$events

  ##step2 : creating the grid
  grid <- build_grid(grid_shape,lines)

  ##step3 : splitting the dataset with each rectangle
  selections <- split_by_grid(grid,samples,events,lines,bw)

  ##step4 : calculating the values

  if(verbose){
    print("starting the kernel calculations ...")
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {

        invisible(capture.output(values <- nkde_worker(sel$lines, sel$events,
                              sel$samples, kernel_func,
                              bw, method, div, digits,
                              tol, max_depth, verbose)))

        df <- data.frame("goid"=sel$samples$goid,
                         "k" = values)
        p(sprintf("i=%g", sel$index))
        return(df)
      })
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {

      invisible(capture.output(values <- nkde_worker(sel$lines, sel$events,
                            sel$samples, kernel_func,
                            bw, method, div, digits,
                            tol, max_depth, verbose)))

      df <- data.frame("goid"=sel$samples$goid,
                       "k" = values)
      return(df)
    })
  }


  if(verbose){
    print("combining the results ...")
  }
  ##step5 : combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]
  return(tot_df$k)

}


