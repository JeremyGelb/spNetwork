#' @title Worker function for adaptive bandwidth for TNDE
#'
#' @description The worker function to calculate Adaptive bandwidths according
#'   to Abramson’s smoothing regimen for TNKDE with a space-time interaction (INTERNAL).
#'
#' @param lines A feature collection of linestrings representing the underlying
#'   network
#' @param quad_events a feature collection of points indicating for which events
#'   the densities must be calculated
#' @param events_loc A feature collection of points representing the location of
#'   the events
#' @param events A feature collection of points representing the events.
#'   Multiple events can share the same location. They are linked by the goid
#'   column
#' @param w A numeric vector with the weight of the events
#' @param kernel_name The name of the kernel to use (string)
#' @param bw_net A float indicating the fixed network bandwidth
#' @param bw_time A float indicating the fixed time bandwidth
#' @param method The type of NKDE to use (string)
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#'   events within the radius around each sampling point), "bw" (the bandwidth)
#'   "none" (the simple sum).
#' @param digits The number of digits to retain from the spatial coordinates. It
#'   ensures that topology is good when building the network. Default is 3. Too
#'   high a precision (high number of digits) might break some connections
#' @param tol A float indicating the minimum distance between the events and the
#'   lines' extremities when adding the point to the network. When points are
#'   closer, they are added at the extremity of the lines.
#' @template verbose-arg
#' @template sparse-arg
#' @param max_depth An integer, the maximum depth to reach for continuous and
#'   discontinuous NKDE
#' @return A vector with the local bandwidths
#' @export
#' @examples
#' #This is an internal function, no example provided
worker_adaptive_bw_tnkde <- function(lines,
                                     quad_events, events_loc, events, w,
                                     kernel_name, bw_net, bw_time, method, div,
                                     digits, tol, sparse, max_depth, verbose = FALSE){

  # if we do not have event in that space, just return NULL
  if(nrow(events)==0){
    return(NULL)
  }

  ## step1 creating the graph
  graph_result <- build_graph(lines,digits = digits,line_weight = "length")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  edges <- graph_result$spedges

  ## step2 finding for each event, its corresponding node
  ## NOTE : there will be less samples than events most of the time
  events_loc$vertex_id <- closest_points(events_loc, nodes)

  events_loc2 <- st_drop_geometry(events_loc)
  events2 <- st_drop_geometry(events)
  quad_events2 <- st_drop_geometry(quad_events)

  #first a join for all the events in the bw
  vertex_id <- NULL # avoid a NOTE
  i.vertex_id <- NULL # avoid a NOTE
  setDT(events2)[events_loc2, on = "goid", vertex_id := i.vertex_id]

  #and a second join for the quad_events
  setDT(quad_events2)[events_loc2, on = "goid", vertex_id := i.vertex_id]

  ## step3 starting the calculations !
  neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")
  neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})


  kernel_values <- adaptive_bw_tnkde_cpp(method = method,
                                         neighbour_list = neighbour_list,
                                         sel_events = quad_events2$vertex_id,
                                         sel_events_wid = quad_events2$wid,
                                         sel_events_time = quad_events2$time,
                                         events = events2$vertex_id,
                                         events_wid = events2$wid,
                                         events_time = events2$time,
                                         weights = w,
                                         bws_net = bw_net,
                                         bws_time = bw_time,
                                         kernel_name = kernel_name,
                                         line_list = graph_result$linelist,
                                         max_depth = max_depth,
                                         min_tol =.Machine$double.xmin)
  ## and adding the self weight
  if (method != "continuous"){
    kfun <- select_kernel(kernel_name)
    kernel_values <- kernel_values + ((kfun(0, bw_net)*kfun(0, bw_time))/(bw_net * bw_time))
  }


  ## at that point, we have a list of numeric vectors or a list of dataframes, one for each bw
  return(kernel_values)
}




#' @title Adaptive bandwidth for TNDE
#'
#' @description Function to calculate Adaptive bandwidths according to Abramson’s smoothing regimen
#' for TNKDE with a space-time interaction.
#'
#' @param grid A spatial grid to split the data within
#' @param events A feature collection of points representing the events points
#' @param lines A feature collection of linestrings representing the network
#' @param bw_time The fixed kernel bandwidth for the time dimension
#' @param bw_net The fixed kernel bandwidth for the network dimension
#' @param trim_bw_time The maximum size of local bandiwidths for time dimension
#' @param trim_bw_net The maximum size of local bandiwidths for network dimension
#' @param method The method to use when calculating the NKDE
#' @param kernel_name The name of the kernel to use
#' @param max_depth The maximum recursion depth
#' @param div The divisor to use for kernels
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @param sparse A Boolean indicating if sparse matrix should be used
#' @param verbose A Boolean indicating if update messages should be printed
#' @return A vector with the local bandwidths
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
adaptive_bw_tnkde <- function(grid, events_loc, events, lines,
                              bw_net, bw_time ,trim_bw_net, trim_bw_time,
                              method, kernel_name, max_depth, div,
                              tol, digits, sparse, verbose){

  ##step 1 split the datas !
  selections <- split_by_grid_abw(grid, events_loc, lines, trim_bw_net, tol, digits)

  ## step 2 calculating the temp NKDE values
  if(verbose){
    print("start calculating the kernel values for the adaptive bandwidth...")
  }
  n_quadra <- length(selections)

  dfs <- lapply(1:n_quadra,function(i){
    sel <- selections[[i]]

    #sel_events <- subset(events, events$goid %in% sel$events$goid)
    sel_events <- subset(events, events$goid %in% sel$samples$goid)

    # nice tnkde worker function !
    values <- worker_adaptive_bw_tnkde(lines = sel$lines,
                                  #quad_events = sel$samples,
                                  quad_events = sel_events,
                                  events_loc = sel$events,
                                  events = sel_events,
                                  w = sel_events$weight,
                                  kernel_name = kernel_name,
                                  bw_net = bw_net, bw_time = bw_time,
                                  method = method, div = div,
                                  digits = digits, tol = tol, sparse = sparse,
                                  max_depth = max_depth, verbose = verbose)
    # values will be a simple numeric vector

    # NOTE : the first column store the goid !
    mat_result <- cbind(sel_events$goid, as.vector(values))
    return(mat_result)
  })

  ## step 3  combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- data.frame(tot_df[order(tot_df[,1]),])
  names(tot_df) <- c("goid", "k")

  ## step 4 calculating the new bandwidth !
  delta <- calc_gamma(tot_df$k)
  new_net_bw <- bw_net * (tot_df$k**(-1/2) * delta**(-1))
  new_net_bw <- ifelse(new_net_bw<trim_bw_net, new_net_bw, trim_bw_net)

  new_time_bw <- bw_time * (tot_df$k**(-1/2) * delta**(-1))
  new_time_bw <- ifelse(new_time_bw<trim_bw_time, new_time_bw, trim_bw_time)
  return(list("bws_net" = new_net_bw,
              "bws_time" = new_time_bw))
}



#' @title Adaptive bandwidth for TNDE (multicore)
#'
#' @description Function to calculate Adaptive bandwidths according to Abramson’s smoothing regimen
#' for TNKDE with a space-time interaction with multicore support.
#'
#' @param grid A spatial grid to split the data within
#' @param events A feature collection of points representing the events points
#' @param lines A feature collection of linestrings representing the network
#' @param bw_time The fixed kernel bandwidth for the time dimension
#' @param bw_net The fixed kernel bandwidth for the network dimension
#' @param trim_bw_time The maximum size of local bandiwidths for time dimension
#' @param trim_bw_net The maximum size of local bandiwidths for network dimension
#' @param method The method to use when calculating the NKDE
#' @param kernel_name The name of the kernel to use
#' @param max_depth The maximum recursion depth
#' @param div The divisor to use for kernels
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @param sparse A Boolean indicating if sparse matrix should be used
#' @param verbose A Boolean indicating if update messages should be printed
#' @return A vector with the local bandwidths
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
adaptive_bw_tnkde.mc <- function(grid, events_loc, events, lines,
                              bw_net, bw_time ,trim_bw_net, trim_bw_time,
                              method, kernel_name, max_depth, div,
                              tol, digits, sparse, verbose){

  ##step 1 split the datas !
  selections <- split_by_grid_abw(grid, events_loc, lines, trim_bw_net, tol, digits)

  ## step 2 calculating the temp NKDE values
  if(verbose){
    print("start calculating the kernel values for the adaptive bandwidth...")
  }

  if(verbose == FALSE){
    dfs <- future.apply::future_lapply(selections, function(sel) {

      sel_events <- subset(events, events$goid %in% sel$samples$goid)

      values <- spNetwork::worker_adaptive_bw_tnkde(lines = sel$lines,
                                         # quad_events = sel_events, #sel_samples
                                         quad_events = sel_events,
                                         events_loc = sel$events,
                                         events = sel_events,
                                         w = sel_events$weight,
                                         kernel_name = kernel_name,
                                         bw_net = bw_net, bw_time = bw_time,
                                         method = method, div = div,
                                         digits = digits, tol = tol, sparse = sparse,
                                         max_depth = max_depth, verbose = FALSE)
      # values will be a simple numeric vector

      # NOTE : the first column store the goid !
      #mat_result <- cbind(sel$samples$goid, as.vector(values))
      mat_result <- cbind(sel_events$goid, as.vector(values))
      return(mat_result)

    }, future.packages = c("spNetwork"))
  }else {
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {

        sel_events <- subset(events, events$goid %in% sel$samples$goid)

        values <- spNetwork::worker_adaptive_bw_tnkde(lines = sel$lines,
                                           quad_events = sel_events,
                                           events_loc = sel$events,
                                           events = sel_events,
                                           w = sel_events$weight,
                                           kernel_name = kernel_name,
                                           bw_net = bw_net, bw_time = bw_time,
                                           method = method, div = div,
                                           digits = digits, tol = tol, sparse = sparse,
                                           max_depth = max_depth, verbose = FALSE)
        # values will be a simple numeric vector

        # NOTE : the first column store the goid !
        mat_result <- cbind(sel_events$goid, as.vector(values))
        p(sprintf("i=%g", sel$index))
        return(mat_result)

      }, future.packages = c("spNetwork"))
    })
  }


  ## step 3  combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- data.frame(tot_df[order(tot_df[,1]),])
  names(tot_df) <- c("goid", "k")

  ## step 4 calculating the new bandwidth !
  delta <- calc_gamma(tot_df$k)
  new_net_bw <- bw_net * (tot_df$k**(-1/2) * delta**(-1))
  new_net_bw <- ifelse(new_net_bw<trim_bw_net, new_net_bw, trim_bw_net)

  new_time_bw <- bw_time * (tot_df$k**(-1/2) * delta**(-1))
  new_time_bw <- ifelse(new_time_bw<trim_bw_time, new_time_bw, trim_bw_time)
  return(list("bws_net" = new_net_bw,
              "bws_time" = new_time_bw))
}

