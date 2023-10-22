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
#' @param bw_time The fixed kernel bandwidth for the time dimension. Can also be a vector
#' if several bandwidth must be used.
#' @param bw_net The fixed kernel bandwidth for the network dimension. Can also be a vector
#' if several bandwidth must be used.
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
#' @return A vector with the local bandwidths or an array if bw_net and bw_time are vectors
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
  graph_result <- build_graph(lines, digits = digits, line_weight = "length")
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
  print(graph)
  neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")
  neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})
  print("here is the neighbour list : ")
  print(neighbour_list)

  kernel_values <- adaptive_bw_tnkde_cpp2(method = method,
                                         neighbour_list = neighbour_list,
                                         sel_events = quad_events2$vertex_id,
                                         sel_events_wid = quad_events2$wid,
                                         sel_events_time = quad_events2$time,
                                         events = events_loc2$vertex_id,
                                         events_wid = events_loc2$wid,
                                         events_time = events_loc2$time,
                                         weights = events_loc2$weight,
                                         bws_net = bw_net,
                                         bws_time = bw_time,
                                         kernel_name = kernel_name,
                                         line_list = graph_result$linelist,
                                         max_depth = max_depth,
                                         min_tol =.Machine$double.xmin)
  # the self weight is already added from
  # the function above if the method is continuous
  if(method != "continuous"){
    ## and adding the self weight
    kfun <- select_kernel(kernel_name)
    bws_prod_mat <- outer(bw_net, bw_time)
    kern_prod_mat <- outer(kfun(0, bw_net), kfun(0, bw_time))
    # NOTE :
    # WHEN two events are merged and get a weight of 2 we must
    # adjust their weight accordingly here
    w2 <- events_loc2$weight[match(quad_events2$goid,events_loc2$goid)]
    kern_prod_arr <- sapply(w2, function(x){
      (x * kern_prod_mat) / bws_prod_mat
    }, simplify = "array")
    #kernel_values <- kernel_values + replicate(dim(kernel_values)[[3]], new_mat)
    kernel_values <- kernel_values + kern_prod_arr
    # HERE : I have validated that I get the same result from the new method with adaptive_bw_tnkde_cpp2
    # the result is an array rather than a vector and can accomodate several bws for network and time
    # VS the old version
    #kernel_values <- kernel_values + ((kfun(0, bw_net)*kfun(0, bw_time))/(bw_net * bw_time))


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
#' @param bw_time The fixed kernel bandwidth for the time dimension. Can also be a vector
#' if several bandwidth must be used.
#' @param bw_net The fixed kernel bandwidth for the network dimension. Can also be a vector
#' if several bandwidth must be used.
#' @param trim_bw_time The maximum size of local bandwidths for time dimension. Must be a vector if bw_net is a vector
#' @param trim_bw_net The maximum size of local bandwidths for network dimension. Must be a vector if bw_net is a vector
#' @param method The method to use when calculating the NKDE
#' @param kernel_name The name of the kernel to use
#' @param max_depth The maximum recursion depth
#' @param div The divisor to use for kernels
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @param sparse A Boolean indicating if sparse matrix should be used
#' @param verbose A Boolean indicating if update messages should be printed
#' @return A vector with the local bandwidths, or an array if bw_time and bw_net are vectors.
#' In that case, the array has the following dimensions : length(bw_net) X length(bw_time) X nrow(events)
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
adaptive_bw_tnkde <- function(grid, events_loc, events, lines,
                              bw_net, bw_time, trim_bw_net, trim_bw_time,
                              method, kernel_name, max_depth, div,
                              tol, digits, sparse, verbose){

  ##step 1 split the datas !
  selections <- split_by_grid_abw(grid, events_loc, lines, max(trim_bw_net), tol, digits)

  ## step 2 calculating the temp NKDE values
  if(verbose){
    print("start calculating the kernel values for the adaptive bandwidth...")
  }
  n_quadra <- length(selections)

  dfs <- lapply(1:n_quadra,function(i){
    print(i)
    sel <- selections[[i]]

    # sel_events gives us the events inside the quadra considered
    sel_events <- subset(events, events$goid %in% sel$samples$goid)

    # nice tnkde worker function !
    values <- worker_adaptive_bw_tnkde(lines = sel$lines,
                                  quad_events = sel_events, # the events for wich we need to calculate the density
                                  events_loc = sel$events, # the locations of all the events in the quadra
                                  events = sel$samples,
                                  w = sel_events$weight,
                                  kernel_name = kernel_name,
                                  bw_net = bw_net,
                                  bw_time = bw_time,
                                  method = method,
                                  div = div,
                                  digits = digits,
                                  tol = tol,
                                  sparse = sparse,
                                  max_depth = max_depth,
                                  verbose = verbose)
    # values will always be an array now
    # if length(bw_net) == 1 and length(bw_time) == 1, then it could simplified to a vector
    # with c(values)

    # # NOTE : the first column store the goid !
    # mat_result <- cbind(sel_events$goid, as.vector(values))
    # return(mat_result)
    # we return a list, the firtst element is the array with the k values
    # and the second is the goid of each event (third dimension of the array)
    return(list(values, sel_events$wid))

  })

  ## step 3  combining the results
  tot_arr <- do.call(abind::abind, lapply(dfs, function(x){x[[1]]}))
  all_wids <- do.call(c, lapply(dfs, function(x){x[[2]]}))

  # tot_df <- do.call(rbind,dfs)
  # tot_df <- data.frame(tot_df[order(tot_df[,1]),])
  # names(tot_df) <- c("goid", "k")

  ## step 4 calculating the new bandwidth !
  final_bws_net <- array(0,dim = dim(tot_arr))
  final_bws_time <- array(0,dim = dim(tot_arr))

  print("here are the estimated densities before calculating local bws")
  print(tot_arr)
  for(i in 1:length(bw_net)){
    for(j in 1:length(bw_time)){
      k <- tot_arr[i,j,]
      k <- k[order(all_wids)]
      delta <- calc_gamma(k)
      new_net_bw <- bw_net[[i]] * (k**(-1/2) * delta**(-1))
      new_net_bw <- ifelse(new_net_bw<trim_bw_net[[i]], new_net_bw, trim_bw_net[[i]])
      final_bws_net[i,j,] <- new_net_bw

      new_time_bw <- bw_time[[j]] * (k**(-1/2) * delta**(-1))
      new_time_bw <- ifelse(new_time_bw<trim_bw_time[[j]], new_time_bw, trim_bw_time[[j]])
      final_bws_time[i,j,] <- new_time_bw
    }
  }

  # delta <- calc_gamma(tot_df$k)
  # new_net_bw <- bw_net * (tot_df$k**(-1/2) * delta**(-1))
  # new_net_bw <- ifelse(new_net_bw<trim_bw_net, new_net_bw, trim_bw_net)
  #
  # new_time_bw <- bw_time * (tot_df$k**(-1/2) * delta**(-1))
  # new_time_bw <- ifelse(new_time_bw<trim_bw_time, new_time_bw, trim_bw_time)

  # let me add a small conversion here to avoid modifying all my code everywhere
  if((length(bw_net) + length(bw_time)) == 2){
    final_bws_net <- c(final_bws_net)
    final_bws_time <- c(final_bws_time)
  }

  return(list("bws_net" = final_bws_net,
              "bws_time" = final_bws_time))
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

  ## step 1 split the data !
  selections <- split_by_grid_abw(grid, events_loc, lines, max(trim_bw_net), tol, digits)

  ## step 2 calculating the temp NKDE values
  if(verbose){
    print("start calculating the kernel values for the adaptive bandwidth...")
  }

  if(verbose == FALSE){
    dfs <- future.apply::future_lapply(selections, function(sel) {

      sel_events <- subset(events, events$goid %in% sel$samples$goid)

      values <- spNetwork::worker_adaptive_bw_tnkde(lines = sel$lines,
                                         # quad_events = sel_events, #sel_samples
                                         quad_events = sel_events, # the events for wich we want to calculate the densities
                                         events_loc = sel$events, # all the events contributing to the densities (possibly outside the quadra)
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
      #return(mat_result)
      return(list(values, sel_events$wid))

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
        #mat_result <- cbind(sel_events$goid, as.vector(values))
        p(sprintf("i=%g", sel$index))
        #return(mat_result)
        return(list(values, sel_events$wid))

      }, future.packages = c("spNetwork"))
    })
  }


  # ## step 3  combining the results
  # tot_df <- do.call(rbind,dfs)
  # tot_df <- data.frame(tot_df[order(tot_df[,1]),])
  # names(tot_df) <- c("goid", "k")
  #
  # ## step 4 calculating the new bandwidth !
  # delta <- calc_gamma(tot_df$k)
  # new_net_bw <- bw_net * (tot_df$k**(-1/2) * delta**(-1))
  # new_net_bw <- ifelse(new_net_bw<trim_bw_net, new_net_bw, trim_bw_net)
  #
  # new_time_bw <- bw_time * (tot_df$k**(-1/2) * delta**(-1))
  # new_time_bw <- ifelse(new_time_bw<trim_bw_time, new_time_bw, trim_bw_time)
  # return(list("bws_net" = new_net_bw,
  #             "bws_time" = new_time_bw))

  ## step 3  combining the results
  tot_arr <- do.call(abind::abind, lapply(dfs, function(x){x[[1]]}))
  all_wids <- do.call(c, lapply(dfs, function(x){x[[2]]}))

  # tot_df <- do.call(rbind,dfs)
  # tot_df <- data.frame(tot_df[order(tot_df[,1]),])
  # names(tot_df) <- c("goid", "k")

  ## step 4 calculating the new bandwidth !
  final_bws_net <- array(0,dim = dim(tot_arr))
  final_bws_time <- array(0,dim = dim(tot_arr))

  print("here are the estimated densities before calculating local bws")
  print(tot_arr)
  for(i in 1:length(bw_net)){
    for(j in 1:length(bw_time)){
      k <- tot_arr[i,j,]
      k <- k[order(all_wids)]
      delta <- calc_gamma(k)
      new_net_bw <- bw_net[[i]] * (k**(-1/2) * delta**(-1))
      new_net_bw <- ifelse(new_net_bw<trim_bw_net[[i]], new_net_bw, trim_bw_net[[i]])
      final_bws_net[i,j,] <- new_net_bw

      new_time_bw <- bw_time[[j]] * (k**(-1/2) * delta**(-1))
      new_time_bw <- ifelse(new_time_bw<trim_bw_time[[j]], new_time_bw, trim_bw_time[[j]])
      final_bws_time[i,j,] <- new_time_bw
    }
  }

  # delta <- calc_gamma(tot_df$k)
  # new_net_bw <- bw_net * (tot_df$k**(-1/2) * delta**(-1))
  # new_net_bw <- ifelse(new_net_bw<trim_bw_net, new_net_bw, trim_bw_net)
  #
  # new_time_bw <- bw_time * (tot_df$k**(-1/2) * delta**(-1))
  # new_time_bw <- ifelse(new_time_bw<trim_bw_time, new_time_bw, trim_bw_time)

  # let me add a small conversion here to avoid modifying all my code everywhere
  if((length(bw_net) + length(bw_time)) == 2){
    final_bws_net <- c(final_bws_net)
    final_bws_time <- c(final_bws_time)
  }

  return(list("bws_net" = final_bws_net,
              "bws_time" = final_bws_time))

}

