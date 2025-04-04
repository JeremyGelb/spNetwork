#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Functions to perform the simple TNKDE ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Simple TNKDE algorithm
#'
#' @description Function to perform the simple tnkde.
#'
#' @param graph a graph object from igraph representing the network
#' @param events a feature collection of points representing the events. It must be
#' snapped on the network, and be nodes of the network. A column vertex_id
#' must indicate for each event its corresponding node
#' @param samples a feature collection of points representing the sampling points.
#' The samples must be snapped on the network. A column edge_id must indicate
#' for each sample on which edge it is snapped.
#' @param samples_time a numeric vector indicating when the densities must be
#' sampled
#' @param bws_net a vector indicating the network kernel bandwidth (in meters) for each
#' event
#' @param bws_time a vector indicating the time kernel bandwidth for each
#' event
#' @param kernel_func a function obtained with the function select_kernel
#' @param nodes a feature collection of points representing the nodes of the network
#' @param edges a feature collection of linestrings representing the edges of the network
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#' events within the radius around each sampling point), "bw" (the bandwidth)
#' "none" (the simple sum).
#' @return a list of two matrices. The first one ins the matrix of the densities,
#' the rows are the samples and the columns the time. The second has the same
#' dimensions and contains the number of events influencing each sample
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
simple_tnkde <- function(graph, events, samples, samples_time, bws_net, bws_time, kernel_func, nodes, edges, div){

  ## step 1 set all values to 0
  base_k <- matrix(0, nrow = nrow(samples), ncol = length(samples_time))
  base_count <- matrix(0, nrow = nrow(samples), ncol = length(samples_time))

  #if the number of event is 0, then return only 0 values
  if(nrow(events)==0){
    return(list("sum_k" = base_k,
                "n" = base_count
                ))
  }

  # sample_tree <- build_quadtree(samples)
  # edges_tree <- build_quadtree(edges)

  # we find the snapped geometries
  event_nodes <- nodes[events$vertex_id,]
  event_nodes$order_id <- as.character(1:nrow(event_nodes))

  # I must precalculate the spatial intersections for the samples for each event
  samples_inter <- st_join(samples, st_buffer(event_nodes['order_id'], bws_net))
  samples_inter <- split(samples_inter, samples_inter$order_id)

  # and for the edges
  edges_inter <- st_join(edges, st_buffer(event_nodes['order_id'], (bws_net + +0.1*bws_net)))
  edges_inter <- split(edges_inter, edges_inter$order_id)



  ## NOTE: we can iterate several time on the same vertex as a start point
  # It would be efficient to store the values of known vertex, but only if
  # we will find this vertex several times
  count_verts <- table(events$vertex_id)
  dupp_verts <- as.numeric(unique(names(count_verts)[count_verts > 1]))
  stored_dens_values <- list()

  ## step2 iterate over each event
  pb <- txtProgressBar(min = 0, max = nrow(events), style = 3)
  for(i in 1:nrow(events)){
    setTxtProgressBar(pb, i)
    bw_net <- bws_net[[i]]
    bw_time <- bws_time[[i]]
    #extracting the starting values
    e <- events[i,]
    y <- e$vertex_id
    w <- e$weight

    selected_samples <- samples_inter[[as.character(i)]]
    selected_edges <- edges_inter[[as.character(i)]]

    # calculating the network density values ! (reusing when possible)
    if(y %in% dupp_verts){
      if(is.null(stored_dens_values[y][[1]])){
        # samples_k <- ess_kernel(graph,y, bw_net, kernel_func, samples, nodes, edges)
    		samples_k <- ess_kernel(graph,y, bw_net, kernel_func, ok_samples = selected_samples,
                                  nodes = nodes, ok_edges = selected_edges, N = nrow(samples))
        stored_dens_values[y] <- samples_k
      }else{
        samples_k <- stored_dens_values[y][[1]]
      }
    }else{
      samples_k <- ess_kernel(graph,y, bw_net, kernel_func, ok_samples = selected_samples,
                              nodes = nodes, ok_edges = selected_edges)
    }

    # calculating the temporal density values
    et <- e$time
    time_k <- kernel_func(abs(et - samples_time), bw_time)

    # applying the scaling if we are applying a bw scaling
    if(div == "bw"){
      time_k <- time_k/bw_time
      samples_k <- samples_k/bw_net
    }

    # creating the matrix of density values
    mat_k <- sapply(time_k, function(x){
      return(x * samples_k)
    })

    samples_count <- ifelse(mat_k>0,1,0)
    base_k <- mat_k * w + base_k
    base_count <- base_count + (samples_count * w)
  }

  # applying the scaling if we are applying a n scaling
  if(div == "n"){
    base_k <- base_k/base_count
  }

  return(list(
    "sum_k" = base_k,
    "n" = base_count
  ))
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### main worker functions for tnkde ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title TNKDE worker
#'
#' @description The worker function for tnkde and tnkde.mc
#'
#' @param lines A feature collection of linestrings with the sampling points. The
#' geometries must be simple Linestrings (may crash if some geometries
#'  are invalid)
#' @param events_loc A feature collection of points representing the aggergated events on the
#' network. The points will be snapped on the network.
#' @param events A feature collection of points representing the base events on the
#' network
#' @param samples_loc A feature collection of points representing the locations for
#' which the densities will be estimated.
#' @param samples_time A numeric vector representing when each density will be
#' estimated
#' @param kernel_name The name of the kernel to use
#' @param bw_net The global network kernel bandwidth
#' @param bw_time The global time kernel bandwidth
#' @param bws_net The network kernel bandwidth (in meters) for each event
#' @param bws_time The time bandwidth for each event
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#' events within the radius around each sampling point), "bw" (the bandwidth)
#' "none" (the simple sum).
#' @param digits The number of digits to keep in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines extremities. When
#' points are closer, they are added at the extremity of the lines.
#' @param sparse A Boolean indicating if sparse or regular matrices should be
#' used by the Rcpp functions. Regular matrices are faster, but require more
#' memory and could lead to error, in particular with multiprocessing. Sparse
#' matrices are slower, but require much less memory.
#' @param max_depth When using the continuous and discontinuous methods, the
#' calculation time and memory use can go wild  if the network has a lot of
#' small edges (area with a lot of intersections and a lot of events). To
#' avoid it, it is possible to set here a maximum depth. Considering that the
#' kernel is divided at intersections, a value of 8 should yield good
#' estimates. A larger value can be used without problem for the discontinuous
#' method. For the continuous method, a larger value will strongly impact
#' calculation speed.
#' @param verbose A Boolean, indicating if the function should print messages
#' about the process.
#' @importFrom igraph adjacent_vertices get.edge.ids
#' @return A numeric matrix with the nkde values
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
tnkde_worker <- function(lines, events_loc, events, samples_loc, samples_time, kernel_name,bw_net, bw_time, bws_net, bws_time, method, div, digits, tol, sparse, max_depth, verbose = FALSE){

  # if we do not have event in that space, just return 0 values
  if(is.null(events)){
    values <- matrix(0, nrow = nrow(samples_loc), ncol = length(samples_time))
    return(values)
  }
  if(nrow(events)==0){
    values <- matrix(0, nrow = nrow(samples_loc), ncol = length(samples_time))
    return(values)
  }

  ## step1 creating the graph
  if(verbose){
    print("    build graph ...")
  }
  graph_result <- build_graph(lines, digits = digits,line_weight = "length")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  edges <- graph_result$spedges

  snapped_samples <- snapPointsToLines2(samples_loc, edges, snap_dist = bw_net, idField = "edge_id")
  samples_loc$edge_id <- snapped_samples$nearest_line_id

  ## step3 finding for each event, its corresponding node
  events_loc$vertex_id <- closest_points(events_loc, nodes)
  events$vertex_id <- events_loc$vertex_id[match(events$goid,events_loc$goid)]

  ## step4 adding the spatial coordinates to samples and nodes
  XY_nodes <- st_coordinates(nodes)
  nodes$X_coords <- XY_nodes[,1]
  nodes$Y_coords <- XY_nodes[,2]

  XY_samples <- st_coordinates(samples_loc)
  samples_loc$X_coords <- XY_samples[,1]
  samples_loc$Y_coords <- XY_samples[,2]

  ## step5 adding a local oid for samples and events
  events$oid <- 1:nrow(events)
  samples_loc$oid <- 1:nrow(samples_loc)

  ## step6 starting the calculations !

  if(verbose){
    print("        calculating NKDE values ...")
  }

  if(method == "simple"){
    # the cas of the simple method, no c++ here
    kernel_func <- select_kernel(kernel_name)
    if(verbose){
      values <- simple_tnkde(graph, events, samples_loc, samples_time, bws_net, bws_time, kernel_func, nodes, edges, div)
    }else{
      invisible(capture.output(values <- simple_tnkde(graph, events, samples_loc, samples_time, bws_net, bws_time, kernel_func, nodes, edges, div)))
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
        values <- spNetwork::continuous_tnkde_cpp_arma_sparse(neighbour_list = neighbour_list,
                                                              events = events$vertex_id,
                                                              weights = events$weight,
                                                              events_time = events$time,
                                                              samples = st_drop_geometry(samples_loc),
                                                              samples_time = samples_time,
                                                              bws_net = bws_net,
                                                              bws_time = bws_time,
                                                              kernel_name = kernel_name,
                                                              nodes = st_drop_geometry(nodes), line_list = graph_result$linelist,
                                                              max_depth = max_depth, verbose = verbose,
                                                              div = div)
      }else{
        values <- spNetwork::continuous_tnkde_cpp_arma(neighbour_list = neighbour_list,
                                                       events = events$vertex_id,
                                                       weights = events$weight,
                                                       events_time = events$time,
                                                       samples = st_drop_geometry(samples_loc),
                                                       samples_time = samples_time,
                                                       bws_net = bws_net,
                                                       bws_time = bws_time,
                                                       kernel_name = kernel_name,
                                                       nodes = st_drop_geometry(nodes),
                                                       line_list = graph_result$linelist,
                                                       max_depth = max_depth, verbose = verbose,
                                                       div = div)
      }

    }

    if(method == "discontinuous"){
      if(sparse){
        values <- spNetwork::discontinuous_tnkde_cpp_arma_sparse(neighbour_list = neighbour_list,
                                                                 events = events$vertex_id,
                                                                 weights = events$weight,
                                                                 events_time = events$time,
                                                                 samples = st_drop_geometry(samples_loc),
                                                                 samples_time = samples_time,
                                                                 bws_net = bws_net,
                                                                 bws_time = bws_time,
                                                                 kernel_name = kernel_name,
                                                                 nodes = st_drop_geometry(nodes),
                                                                 line_list = graph_result$linelist,
                                                                 max_depth = max_depth, verbose = verbose,
                                                                 div = div)
      }else{
        values <- spNetwork::discontinuous_tnkde_cpp_arma(neighbour_list = neighbour_list,
                                                          events = events$vertex_id,
                                                          weights = events$weight,
                                                          events_time = events$time,
                                                          samples = st_drop_geometry(samples_loc),
                                                          samples_time = samples_time,
                                                          bws_net = bws_net,
                                                          bws_time = bws_time,
                                                          kernel_name = kernel_name,
                                                          nodes = st_drop_geometry(nodes), line_list = graph_result$linelist,
                                                          max_depth = max_depth, verbose = verbose,
                                                          div = div)
      }


    }

  }

  ## step7 adjusting the kernel values !
  ## NOTE : the values are already adjusted !
  return(values$sum_k)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### main function for tnkde ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Temporal Network Kernel density estimate
#'
#' @description Calculate the Temporal Network Kernel Density Estimate based on a network of lines,
#' sampling points in space and times, and events in space and time.
#'
#' @details
#' **Temporal Network Kernel Density Estimate**\cr
#' The TNKDE is an extension of the NKDE considering both the location of events on the network and
#' in time. Thus, density estimation (density sampling) can be done along lines of the network and
#' at different time. It can be used with the three NKDE (simple, discontinuous and continuous).
#' \cr\cr
#' **density in time and space**\cr
#' Two bandwidths must be provided, one for the network distance and one for the
#' time distance. They are both used to calculate the contribution of each event
#' to each sampling point. Let us consider one event E and a sample S. dnet(E,S)
#' is the contribution to network density of E at S location and dtime(E,S) is
#' the contribution to time density of E at S time. The total contribution is
#' thus dnet(E,S) * dtime(E,S). If one of the two densities is 0, then the total
#' density is 0 because the sampling point is out of the covered area by the
#' event in time or in the network space.
#' \cr\cr
#'
#' **adaptive bandwidth**\cr
#' It is possible to use an adaptive bandwidth both on the network and in time.
#' Adaptive bandwidths are calculated using the Abramson’s smoothing regimen
#' \insertCite{abramson1982bandwidth}{spNetwork}. To do so, the original fixed
#' bandwidths must be specified (bw_net and bw_time parameters).
#' The maximum size of the two local bandwidths can be limited with
#' the parameters trim_bw_net and trim_bw_time.
#' \cr\cr
#'
#' **Diggle correction factor**\cr
#' A set of events can be limited in both space (limits of the study
#' area) and time ( beginning and ending of the data collection period). These
#' limits induce lower densities at the border of the set of events, because
#' they are not sampled outside the limits. It is possible to apply the Diggle
#' correction factor \insertCite{diggle1985kernel}{spNetwork} in both the
#' network and time spaces to minimize this effect.
#' \cr\cr
#'
#' **Separated or simultaneous adaptive bandwidth**\cr
#' When the parameter adaptive is TRUE, one can choose between using separated
#' calculation of network and temporal bandwidths, and calculating them
#' simultaneously. In the first case (default), the network bandwidths are
#' determined for each event by considering only their locations and the time
#' bandwidths are determined by considering only there time stamps. In the second
#' case, for each event, the spatio-temporal density at its location on the
#' network and in time is estimated and used to determine both the network and
#' temporal bandwidths. This second approach must be preferred if the events are
#' characterized by a high level of spatio-temporal autocorrelation.
#'
#' @template nkde_params-arg
#' @template diggle_corr-arg
#' @template tnkde-args
#' @param samples_loc A feature collection of points representing the locations for
#' which the densities will be estimated.
#' @param samples_time A numeric vector indicating when the densities will be sampled
#' @template nkde_geoms-args
#' @template grid_shape-arg
#' @template sparse-arg
#' @template check-arg
#' @template verbose-arg
#' @param adaptive_separate A boolean indicating if the adaptive bandwidths
#'   for the time and the network dimensions must be calculated separately (TRUE) or in
#'   interaction (FALSE)
#' @return A matrix with the estimated density for each sample point (rows) at
#'   each timestamp (columns). If adaptive = TRUE, the function returns a list
#'   with two slots: k (the matrix with the density values) and events (a
#'   feature collection of points with the local bandwidths).
#' @export
#' @examples
#' \donttest{
#' # loading the data
#' data(mtl_network)
#' data(bike_accidents)
#'
#' # converting the Date field to a numeric field (counting days)
#' bike_accidents$Time <- as.POSIXct(bike_accidents$Date, format = "%Y/%m/%d")
#' start <- as.POSIXct("2016/01/01", format = "%Y/%m/%d")
#' bike_accidents$Time <- difftime(bike_accidents$Time, start, units = "days")
#' bike_accidents$Time <- as.numeric(bike_accidents$Time)
#'
#' # creating sample points
#' lixels <- lixelize_lines(mtl_network, 50)
#' sample_points <- lines_center(lixels)
#'
#' # choosing sample in times (every 10 days)
#' sample_time <- seq(0, max(bike_accidents$Time), 10)
#'
#' # calculating the densities
#' tnkde_densities <- tnkde(lines = mtl_network,
#'     events = bike_accidents, time_field = "Time",
#'     w = rep(1, nrow(bike_accidents)),
#'     samples_loc = sample_points,
#'     samples_time = sample_time,
#'     kernel_name = "quartic",
#'     bw_net = 700, bw_time = 60, adaptive = TRUE,
#'     trim_bw_net = 900, trim_bw_time = 80,
#'     method = "discontinuous", div = "bw",
#'     max_depth = 10, digits = 2, tol = 0.01,
#'     agg = 15, grid_shape = c(1,1),
#'     verbose  = FALSE)
#'}
tnkde <- function(lines, events, time_field, w, samples_loc, samples_time, kernel_name, bw_net, bw_time, adaptive=FALSE, adaptive_separate = TRUE, trim_bw_net=NULL, trim_bw_time=NULL, method, div="bw", diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE, grid_shape=c(1,1), verbose=TRUE, check=TRUE){

  ## step0 basic checks
  if(verbose){
    print("checking inputs ...")
  }
  if((kernel_name %in% c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov','uniform'))==FALSE){
    stop('the name of the kernel function must be one of c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov" ,"uniform")')
  }

  if(method %in% c("simple","continuous","discontinuous") == FALSE){
    stop('the method must be one of c("simple","continuous","discontinuous"')
  }
  if(method == "continuous" & kernel_name == "gaussian"){
    stop("using the continuous NKDE and the gaussian kernel function can yield negative values for densities because the gaussian kernel does not integrate to 1 within the bandiwdth, please consider using the quartic kernel instead")
  }

  if(min(bw_net)<=0 | min(bw_time) <= 0){
    stop("the network and time bandwidths for the kernel must be superior to 0")
  }

  if(adaptive & (length(bw_net) > 1 | length(bw_time) > 1)){
    stop("When adaptive is TRUE, only a global bandwidth must be given for bw_net and bw_time")
  }

  if(adaptive & (is.null(trim_bw_net) | is.null(trim_bw_time))){
    stop("if adaptive is TRUE, a value for trim_bw must be supplied")
  }
  if(diggle_correction & is.null(study_area)){
    stop("the study_area must be defined if the Diggle correction factor is used")
  }
  if(check){
    check_geometries(lines, samples_loc, events, study_area)
  }
  if(time_field %in% names(events) == FALSE){
    stop(paste0("The column ",time_field," is not a column of events"))
  }

  events$time <- events[[time_field]]
  events$wid <- 1:nrow(events)


  ## step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
  }
  # NOTE : after data preparation, the aggregated events must be kept
  # for creating nodes in the network, but must remain separated when
  # we will iterate over them.

  data <- prepare_data(samples_loc, lines, events, w , digits, tol, agg)
  lines <- data$lines
  samples_loc <- data$samples
  events_loc <- data$events
  events$weight <- w
  events$bws_net <- bw_net
  events$bws_time <- bw_time

  # I must find for each original event, which new location it matches (after aggregation)
  #idx <- FNN::knnx.index(st_coordinates(events_loc),st_coordinates(events), k = 1)
  idx <- closest_points(events, events_loc)
  events$goid <- events_loc$goid[idx]

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples_loc,events))

  ## adaptive bandwidth !
  if(adaptive==FALSE){
    bws_net <- events$bws_net
    bws_time <- events$bws_time
  }else{
    if(adaptive_separate == TRUE){
      ## we want to use an adaptive bw in the network space
      bws_net_all <- adaptive_bw.mc(grid = grid,
                             events = events_loc,
                             lines = lines,
                             bw = bw_net,
                             trim_bw = trim_bw_net,
                             method,
                             kernel_name, max_depth, tol, digits, sparse, verbose)
      ## and in the time space
      bws_time <- adaptive_bw_1d(events$time, w, bw_time, kernel_name)
      bws_net <- bws_net_all[events$goid]
    }else{
      interaction_bws <- adaptive_bw_tnkde(grid, events_loc, events, lines,
                                           bw_net, bw_time ,trim_bw_net, trim_bw_time,
                                           method, kernel_name, max_depth, div,
                                           tol, digits, sparse, verbose)
      bws_time <- interaction_bws$bws_time
      bws_net <- interaction_bws$bws_net
    }
  }

  ## calculating the correction factor
  if(diggle_correction){
    if(verbose){
      print("Calculating the correction factor")
    }
    # corr_factor_net <- correction_factor(study_area,events_loc,lines,method, bws_net, kernel_name, tol, digits, max_depth, sparse)
    # corr_factor_time <- correction_factor_time(events[[time_field]], samples_time, bws_time, kernel_name)
    # # the final density is (dnet*corr_net * dtime*corr_time)*w
    # # because we only have multiplications, I can do the product here
    # corr_factor <- corr_factor_net * corr_factor_time

    outside_mass_time <- correction_factor_time(events[[time_field]], samples_time, bws_time, kernel_name)
    corr_factor_net <- correction_factor(study_area,events_loc,lines,method, bws_net, kernel_name, tol, digits, max_depth, sparse)
    corr_factor_time <- correction_factor_time(events[[time_field]], samples_time, bws_time, kernel_name)
    # the final density is (dnet*corr_net * dtime*corr_time)*w
    # because we only have multiplications, I can do the product here
    # note : mutliplying the correction factor is equivalent to calculate the total inside mass
    # and taking its inverse
    corr_factor <- corr_factor_net * corr_factor_time

    # OLD MASTER VERSION
    # outside_mass_net <- 1-(1/(corr_factor_net))
    # corr_factor <- 1/(1-(outside_mass_net * outside_mass_time))

  }else{
    corr_factor <- rep(1,nrow(events))
  }
  ## duplicating the bws and correction factor if needed for aggregated events
  corr_factor <- corr_factor[events$goid]

  events$weight <- events$weight * corr_factor

  events$bw_net <- bws_net
  events$bw_time <- bws_time

  max_bw <- max(bws_net)
  ## step3 splitting the dataset with each rectangle
  ## NOTE: WE ARE SPLITTING THE LINES HERE
  selections <- split_by_grid(grid, samples_loc, events_loc, lines, max_bw, tol, digits)

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

    # in selection we have the events location, but I also need the base events
    sel_events <- subset(events, events$goid %in% sel$events$goid)

    # nice tnkde worker function !
    values <- tnkde_worker(lines = sel$lines,
                          events_loc = sel$events,
                          events = sel_events,
                          samples_loc = sel$samples,
                          samples_time = samples_time,
                          kernel_name = kernel_name,
                          bw_net = max(bw_net),
                          bw_time = max(bw_time),
                          bws_net = sel_events$bw_net,
                          bws_time = sel_events$bw_time,
                          method = method,
                          div = div,
                          digits = digits,
                          tol = tol,
                          sparse = sparse,
                          max_depth = max_depth,
                          verbose = verbose)

    # NOTE : the first column store the goid !
    mat_result <- cbind(sel$samples$goid, values)

    return(mat_result)
  })

  if(verbose){
    print("combining the results ...")
  }

  ## step5  combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df[,1]),]
  ## oups, case were we have only one sample point...
  if(is.null(dim(tot_df))){
    dim(tot_df) <- c(1,length(tot_df))
  }

  if(adaptive){
    return(list("events" = events,
                "k" = tot_df[,2:ncol(tot_df)]))
  }else{
    return(tot_df[,2:ncol(tot_df)])
  }
}


#' @title Temporal Network Kernel density estimate (multicore)
#'
#' @description Calculate the Temporal Network Kernel Density Estimate based on a network of lines,
#' sampling points in space and times, and events in space and time with multicore support.
#'
#' @details
#' For details, see help(tnkde) and help(nkde)
#'
#' @template nkde_params-arg
#' @template diggle_corr-arg
#' @template tnkde-args
#' @param samples_loc A feature collection of points representing the locations for
#' which the densities will be estimated.
#' @param samples_time A numeric vector indicating when the densities will be sampled
#' @template nkde_geoms-args
#' @template grid_shape-arg
#' @template sparse-arg
#' @template check-arg
#' @template verbose-arg
#' @param adaptive_separate A boolean indicating if the adaptive bandwidths for
#'   the time and the network dimensions must be calculated separately (TRUE) or
#'   in interaction (FALSE)
#' @return A matrix with the estimated density for each sample point (rows) at
#'   each timestamp (columns). If adaptive = TRUE, the function returns a list
#'   with two slots: k (the matrix with the density values) and events (a
#'   feature collection of points with the local bandwidths).
#' @export
#' @examples
#' \donttest{
#' # loading the data
#' data(mtl_network)
#' data(bike_accidents)
#'
#' # converting the Date field to a numeric field (counting days)
#' bike_accidents$Time <- as.POSIXct(bike_accidents$Date, format = "%Y/%m/%d")
#' start <- as.POSIXct("2016/01/01", format = "%Y/%m/%d")
#' bike_accidents$Time <- difftime(bike_accidents$Time, start, units = "days")
#' bike_accidents$Time <- as.numeric(bike_accidents$Time)
#'
#' # creating sample points
#' lixels <- lixelize_lines(mtl_network, 50)
#' sample_points <- lines_center(lixels)
#'
#' # choosing sample in times (every 10 days)
#' sample_time <- seq(0, max(bike_accidents$Time), 10)
#'
#' future::plan(future::multisession(workers=1))
#'
#' # calculating the densities
#' tnkde_densities <- tnkde.mc(lines = mtl_network,
#'     events = bike_accidents, time_field = "Time",
#'     w = rep(1, nrow(bike_accidents)),
#'     samples_loc = sample_points,
#'     samples_time = sample_time,
#'     kernel_name = "quartic",
#'     bw_net = 700, bw_time = 60, adaptive = TRUE,
#'     trim_bw_net = 900, trim_bw_time = 80,
#'     method = "discontinuous", div = "bw",
#'     max_depth = 10, digits = 2, tol = 0.01,
#'     agg = 15, grid_shape = c(1,1),
#'     verbose  = FALSE)
#'
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'}
tnkde.mc <- function(lines, events, time_field, w, samples_loc, samples_time, kernel_name, bw_net, bw_time, adaptive=FALSE, adaptive_separate = TRUE, trim_bw_net=NULL, trim_bw_time=NULL, method, div="bw", diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE, grid_shape=c(1,1), verbose=TRUE, check=TRUE){

  ## step0 basic checks
  if(verbose){
    print("checking inputs ...")
  }
  if((kernel_name %in% c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov','uniform'))==FALSE){
    stop('the name of the kernel function must be one of c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov" ,"uniform")')
  }

  if(method %in% c("simple","continuous","discontinuous") == FALSE){
    stop('the method must be one of c("simple","continuous","discontinuous"')
  }
  if(method == "continuous" & kernel_name == "gaussian"){
    stop("using the continuous NKDE and the gaussian kernel function can yield negative values for densities because the gaussian kernel does not integrate to 1 within the bandiwdth, please consider using the quartic kernel instead")
  }

  if(min(bw_net)<=0 | min(bw_time) <= 0){
    stop("the network and time bandwidths for the kernel must be superior to 0")
  }

  if(adaptive & (length(bw_net) > 1 | length(bw_time) > 1)){
    stop("When adaptive is TRUE, only a global bandwidth must be given for bw_net and bw_time")
  }

  if(adaptive & (is.null(trim_bw_net) | is.null(trim_bw_time))){
    stop("if adaptive is TRUE, a value for trim_bw must be supplied")
  }
  if(diggle_correction & is.null(study_area)){
    stop("the study_area must be defined if the Diggle correction factor is used")
  }
  if(check){
    check_geometries(lines, samples_loc, events, study_area)
  }
  if(time_field %in% names(events) == FALSE){
    stop(paste0("The column ",time_field," is not a column of events"))
  }

  events$time <- events[[time_field]]
  events$wid <- 1:nrow(events)

  ## step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
  }
  # NOTE : after data preparation, the aggregated events must be kept
  # for creating nodes in the network, but must remain separated when
  # we will iterate over them.

  data <- prepare_data(samples_loc, lines, events, w , digits, tol, agg)
  lines <- data$lines
  samples_loc <- data$samples
  events_loc <- data$events
  events$weight <- w
  events$bws_net <- bw_net
  events$bws_time <- bw_time

  # I must find for each original event, which new location it matches (after aggregation)
  #idx <- FNN::knnx.index(st_coordinates(events_loc),st_coordinates(events), k = 1)
  idx <- closest_points(events, events_loc)
  events$goid <- events_loc$goid[idx]

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples_loc,events))

    ## adaptive bandwidth !
  if(adaptive==FALSE){
    bws_net <- events$bws_net
    bws_time <- events$bws_time
  }else{
    if(adaptive_separate  == TRUE){
      ## we want to use an adaptive bw in the network space
      bws_net <- adaptive_bw.mc(grid, events_loc, lines, bw_net, trim_bw_net, method,
                                kernel_name, max_depth, tol, digits, sparse, verbose)
      ## and in the time space
      bws_time <- adaptive_bw_1d(events$time, w, bw_time, kernel_name)
      bws_net <- bws_net[events$goid]
    }else{
      interaction_bws <- adaptive_bw_tnkde.mc(grid, events_loc, events, lines,
                                           bw_net, bw_time ,trim_bw_net, trim_bw_time,
                                           method, kernel_name, max_depth, div,
                                           tol, digits, sparse, verbose)
      bws_time <- interaction_bws$bws_time
      bws_net <- interaction_bws$bws_net
    }
  }

  ## calculating the correction factor
  if(diggle_correction){
    if(verbose){
      print("Calculating the correction factor")
    }
    # corr_factor_net <- correction_factor(study_area,events_loc,lines,method, bws_net, kernel_name, tol, digits, max_depth, sparse)
    # corr_factor_time <- correction_factor_time(events[[time_field]], samples_time, bws_time, kernel_name)
    # the final correction factor is 1 / (1 - outside density)
    # and ontisde density is outside_time * outside_network
    # because we only have multiplications, I can do the product here

    outside_mass_time <- correction_factor_time(events[[time_field]], samples_time, bws_time, kernel_name)
    corr_factor_net <- correction_factor(study_area,events_loc,lines,method, bws_net, kernel_name, tol, digits, max_depth, sparse)
    # NOTE : we need to convert back the correction factor to mass outside network
    outside_mass_net <- 1-(1/(corr_factor_net))
    corr_factor <- 1/(1-(outside_mass_net * outside_mass_time))
  }else{
    corr_factor <- rep(1,nrow(events))
  }
  ## duplicating the bws and correction factor if needed for aggregated events
  corr_factor <- corr_factor[events$goid]

  events$weight <- events$weight * corr_factor

  events$bw_net <- bws_net
  events$bw_time <- bws_time

  max_bw <- max(bws_net)
  ## step3 splitting the dataset with each rectangle
  selections <- split_by_grid.mc(grid, samples_loc, events_loc, lines, max_bw, tol, digits)

  ## step 4 calculating the values
  # save.image(".Rproj.user/error_occured_here.rda")
  if(verbose){
    print("starting the kernel calculations ...")
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {

        # in selection we have the events location, but I also need the base events
        sel_events <- subset(events, events$goid %in% sel$events$goid)

        # TODO : a nice tnkde worker function !
        values <- tnkde_worker(lines = sel$lines,
                               events_loc = sel$events,
                               events = sel_events,
                               samples_loc = sel$samples,
                               samples_time = samples_time,
                               kernel_name = kernel_name,
                               bw_net = bw_net,
                               bw_time = bw_time,
                               bws_net = sel_events$bw_net,
                               bws_time = sel_events$bw_time,
                               method = method,
                               div = div,
                               digits = digits,
                               tol = tol,
                               sparse = sparse,
                               max_depth = max_depth,
                               verbose = verbose)

        # NOTE : the first column store the goid !
        mat_result <- cbind(sel$samples$goid, values)
        p(sprintf("i=%g", sel$index))

        return(mat_result)
      })
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {

      # in selection we have the events location, but I also need the base events
      sel_events <- subset(events, events$goid %in% sel$events$goid)

      # TODO : a nice tnkde worker function !
      values <- tnkde_worker(lines = sel$lines,
                             events_loc = sel$events,
                             events = sel_events,
                             samples_loc = sel$samples,
                             samples_time = samples_time,
                             kernel_name = kernel_name,
                             bw_net = bw_net,
                             bw_time = bw_time,
                             bws_net = sel_events$bw_net,
                             bws_time = sel_events$bw_time,
                             method = method,
                             div = div,
                             digits = digits,
                             tol = tol,
                             sparse = sparse,
                             max_depth = max_depth,
                             verbose = verbose)

      # NOTE : the first column store the goid !
      mat_result <- cbind(sel$samples$goid, values)
      return(mat_result)
    })
  }

  if(verbose){
    print("combining the results ...")
  }

  ## step5  combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df[,1]),]
  ## oups, case were we have only one sample point...
  if(is.null(dim(tot_df))){
    dim(tot_df) <- c(1,length(tot_df))
  }

  if(adaptive){
    return(list("events" = events,
                "k" = tot_df[,2:ncol(tot_df)]))
  }else{
    return(tot_df[,2:ncol(tot_df)])
  }
}


