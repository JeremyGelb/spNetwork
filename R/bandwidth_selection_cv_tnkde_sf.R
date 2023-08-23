#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### shared functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Check function for parameters in bandwidth selection methods
#'
#' @description A check function for bandwidth selection methods raising an error if a parameter is not valid
#'
#' @param check A boolean indicating if the geometries must be checked
#' @param lines A feature collection of linestrings representing the underlying network
#' @param samples A feature collection of points representing the sample location
#' @param events a feature collection of points representing the events
#' @param kernel_name The name of the kernel to use
#' @param method The name of the NKDE to use
#' @template bw_tnkde_selection-args
#' @template diggle_corr-arg
#' @keywords internal
#' @examples
#' # no example provided, this is an internal function
bw_checks <- function(check,lines,samples,events,
                      kernel_name, method, bw_net_range = NULL, bw_time_range = NULL,
                      bw_net_step = NULL, bw_time_step = NULL,
                      diggle_correction = FALSE, study_area = NULL){

  if((kernel_name %in% c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov','uniform'))==FALSE){
    stop('the name of the kernel function must be one of c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov" ,"uniform")')
  }

  if(method %in% c("simple","continuous","discontinuous") == FALSE){
    stop('the method must be one of c("simple","continuous","discontinuous"')
  }
  if(method == "continuous" & kernel_name == "gaussian"){
    stop("using the continuous NKDE and the gaussian kernel function can yield negative values for densities because the gaussian kernel does not integrate to 1 within the bandwidth, please consider using the quartic kernel instead")
  }

  if(is.null(bw_net_range) == FALSE){
    if(min(bw_net_range)<=0){
      stop("the bandwidths for the kernel must be superior to 0")
    }
  }
  if(is.null(bw_time_range) == FALSE){
    if(min(bw_time_range)<=0){
      stop("the bandwidths for the kernel must be superior to 0")
    }
  }

  if(is.null(bw_net_step) == FALSE){
    if(bw_net_step<=0){
      stop("the step between two bandwidths must be greater than 0")
    }
  }
  if(is.null(bw_time_step) == FALSE){
    if(bw_time_step<=0){
      stop("the step between two bandwidths must be greater than 0")
    }
  }

  if(diggle_correction & is.null(study_area)){
    stop("the study_area must be defined if the Diggle correction factor is used")
  }
  if(check){
    check_geometries(lines,samples,events,study_area)
  }
}


#' @title Time and Network bandwidth correction calculation
#'
#' @description Caclulating the border correction factor for both time and network bandwidths
#'
#' @param net_bws A vector of network bandwidths
#' @param time_bws A vector of time bandwidths
#' @template diggle_corr-arg
#' @param events A feature collection of points representing the events
#' @param events_loc A feature collection of points representing the unique location of events
#' @param lines A feature collection of linestrings representing the underlying lines of the network
#' @param method The name of a NKDE method
#' @param kernel_name The name of the kernel to use
#' @param events A feature collection of points representing the events
#' @param kernel_name The name of the kernel to use
#' @param method The name of the NKDE to use
#' @param tol  float indicating the minimum distance between the events and the
#'   lines' extremities when adding the point to the network. When points are
#'   closer, they are added at the extremity of the lines.
#' @param digits An integer, the number of digits to keep for the spatial coordinates
#' @param max_depth The maximal depth for continuous or discontinuous NKDE
#' @template sparse-arg
#' @keywords internal
#' @examples
#' # no example provided, this is an internal function
bw_tnkde_corr_factor <- function(net_bws, time_bws, diggle_correction, study_area, events, events_loc, lines,
                                 method, kernel_name, tol, digits, max_depth, sparse){
  net_bws_corr <- lapply(net_bws, function(bw){
    if(diggle_correction){
      bws <- rep(bw,nrow(events_loc))
      # network corr_factor
      corr_factor <- correction_factor(study_area, events_loc, lines, method, bws,
                                       kernel_name, tol, digits, max_depth, sparse)

      corr_factor <- corr_factor[events$goid]
    }else{
      corr_factor<- rep(1,nrow(events))
    }
    return(corr_factor)
  })

  ## calculating time corr_factors

  time_bws_corr <- lapply(time_bws, function(bw){
    if(diggle_correction){
      bws <- rep(bw,nrow(events))
      # network corr_factor
      corr_factor <- correction_factor_time(events$time, events$time, bws, kernel_name)
    }else{
      corr_factor<- rep(1,nrow(events))
    }
    return(corr_factor)
  })
  return(list(net_bws_corr, time_bws_corr))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### The single core version ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Bandwidth selection by likelihood cross validation for temporal NKDE
#'
#' @description Calculate for multiple network and time bandwidths the cross validation likelihood to
#' select an appropriate bandwidth in a data-driven approach
#'
#' @details  The function calculates the likelihood cross validation score for several time and network
#' bandwidths in order to find the most appropriate one. The general idea is to find the pair of
#' bandwidths that would produce the most similar results if one event is removed from
#' the dataset (leave one out cross validation). We use here the shortcut formula as
#' described by the package spatstat \insertCite{spatstatpkg}{spNetwork}.
#'
#' LCV(h) = sum[i] log(lambda[-i](x[i]))
#'
#' Where the sum is taken for all events x[i] and where lambda[-i](x[i]) is the leave-one-out kernel
#' estimate at x[i] for a bandwidth h. A lower value indicates a better bandwidth.
#'
#' @references{
#'     \insertAllCited{}
#' }
#'
#' @template bw_tnkde_selection-args
#' @template nkde_params-arg
#' @template diggle_corr-arg
#' @template nkde_geoms-args
#' @param time_field The name of the field in events indicating when the events occurred. It must be a numeric field
#' @template sparse-arg
#' @template grid_shape-arg
#' @param sub_sample A float between 0 and 1 indicating the percentage of quadra
#' to keep in the calculus. For large datasets, it may be useful to limit the
#' bandwidth evaluation and thus reduce calculation time.
#' @template verbose-arg
#' @template check-arg
#' @param zero_strat A string indicating what to do when density is 0 when calculating LOO density estimate for an isolated event.
#' "min_double" (default) replace the 0 value by the minimum double possible on the machine. "remove" will remove them from the final
#' score. The first approach penalizes more strongly the small bandwidths.
#' @return A matrix with the cross validation score.  Each row corresponds to a network
#' bandwidth and each column to a time bandwidth (the higher the better).
#' @export
#' @examples
#' \donttest{
#' # loading the data
#' data(mtl_network)
#' data(bike_accidents)
#'
#' # converting the Date field to a numeric field (counting days)
#' bike_accidents$Time <- as.POSIXct(bike_accidents$Date, format = "%Y/%m/%d")
#' bike_accidents$Time <- difftime(bike_accidents$Time, min(bike_accidents$Time), units = "days")
#' bike_accidents$Time <- as.numeric(bike_accidents$Time)
#' bike_accidents <- subset(bike_accidents, bike_accidents$Time>=89)
#'
#' # calculating the cross validation values
#' cv_scores <- bws_tnkde_cv_likelihood_calc(
#'   bw_net_range = c(100,1000),
#'   bw_net_step = 100,
#'   bw_time_range = c(10,60),
#'   bw_time_step = 5,
#'   lines = mtl_network,
#'   events = bike_accidents,
#'   time_field = "Time",
#'   w = rep(1, nrow(bike_accidents)),
#'   kernel_name = "quartic",
#'   method = "discontinuous",
#'   diggle_correction = FALSE,
#'   study_area = NULL,
#'   max_depth = 10,
#'   digits = 2,
#'   tol = 0.1,
#'   agg = 15,
#'   sparse=TRUE,
#'   grid_shape=c(1,1),
#'   sub_sample=1,
#'   verbose = FALSE,
#'   check = TRUE)
#'}
bws_tnkde_cv_likelihood_calc <- function(bw_net_range, bw_net_step,
                                         bw_time_range, bw_time_step,
                                         lines, events, time_field,
                                         w, kernel_name, method,
                                         diggle_correction = FALSE, study_area = NULL,
                                         max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE,
                                         zero_strat = "min_double",
                                         grid_shape=c(1,1), sub_sample=1, verbose=TRUE, check=TRUE){

  ## step0 basic checks
  samples <- events
  if(time_field %in% names(events) == FALSE){
    stop("time_field must be the name of a numeric column in events")
  }
  events$time <- events[[time_field]]
  events$weight <- w
  div <- "bw"
  events$wid <- 1:nrow(events)

  ## checking inputs
  bw_checks(check, lines, samples, events,
                  kernel_name, method, bw_net_range = bw_net_range, bw_time_range = bw_time_range,
                  bw_net_step = bw_net_step, bw_time_step = bw_time_step,
                  diggle_correction = diggle_correction, study_area = study_area)

  if(zero_strat %in% c("min_double", "remove") == FALSE){
    stop("zero_strat argument must be one of c('min_double', 'remove')")
  }

  ## step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
  }

  data <- prepare_data(samples, lines, events, w , digits,tol,agg)
  lines <- data$lines
  events_loc <- data$events

  #idx <- FNN::knnx.index(st_coordinates(events_loc),st_coordinates(events), k = 1)
  idx <- closest_points(events, events_loc)
  events$goid <- events_loc$goid[idx]

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))

  ## calculating the correction factor for each bw combinaisons
  net_bws <- seq(from = bw_net_range[[1]], to = bw_net_range[[2]], by = bw_net_step)
  time_bws <- seq(from = bw_time_range[[1]], to = bw_time_range[[2]], by = bw_time_step)

  if(verbose){
    print("Calculating the correction factor if required")
  }

  ## calculating network corr_factors
  corr_factors <- bw_tnkde_corr_factor(net_bws, time_bws, diggle_correction, study_area, events, events_loc, lines,
                                   method, kernel_name, tol, digits, max_depth, sparse)
  net_bws_corr <- corr_factors[[1]]
  time_bws_corr <- corr_factors[[2]]


  ## NB : the weights can change because of the different BW in time and space
  ## The weights will be passed to c++, to is must be in an easy format, like an array
  ## event_weights(rows = bws_net, cols = bws_time, slices = events)
  events_weight <- array(0, dim = c(length(net_bws), length(time_bws), nrow(events)))

  for (i in 1:length(time_bws_corr)){
    bw_time <- time_bws[[i]]
    corr_time <- time_bws_corr[[i]]
    for (j in 1:length(net_bws_corr)){
      bw_net <- net_bws[[j]]
      corr_net <- net_bws_corr[[j]]
      events_weight[j,i,] <- events$weight * corr_time * corr_net
    }
  }

  max_bw_net <- max(net_bws)

  ## step3 splitting the dataset with each rectangle
  # NB : here we select the events in the gris (samples) and the events locations in the buffer (events_loc)
  selections <- split_by_grid(grid, events, events_loc, lines,max_bw_net, tol, digits, split_all = FALSE)

  ## sub sampling the quadra if required
  if (sub_sample < 1){
    nb <- ceiling(length(selections) * sub_sample)
    selections <- selections[sample(1:length(selections),size = nb,replace = FALSE)]
  }

  ## step 4 calculating the CV values
  if(verbose){
    print("start calculating the CV values ...")
  }

  n_quadra <- length(selections)

  if (verbose){
    pb <- txtProgressBar(min = 0, max = n_quadra, style = 3)
  }
  dfs <- lapply(1:n_quadra,function(i){

    sel <- selections[[i]]

    # the events_loc must cover the quadra and the bw
    sel_events_loc <- sel$events

    # idem for all the events
    sel_events <- subset(events, events$goid %in% sel_events_loc$goid)

    # but I also need to know on which events I must calculate the densities (in the quadra)
    quad_events <- sel$samples
    sel_weights <- events_weight[,,sel_events$wid]

    values <- tnkde_worker_bw_sel(sel$lines,
                                  quad_events,
                                  sel_events_loc,
                                  sel_events,
                                  sel_weights,
                                  kernel_name, net_bws,
                                  time_bws,
                                  method, div, digits,
                                  tol,sparse, max_depth, verbose)


    if(verbose){
      setTxtProgressBar(pb, i)
    }
    return(values)
  })

  # removing NULL elements in list of cubes
  dfs[sapply(dfs, is.null)] <- NULL
  dfs$along <- 3
  final_array <- do.call(abind::abind, dfs)
  if(zero_strat == "min_double"){
    bin_arr <- final_array == 0
    final_array[bin_arr] <- .Machine$double.xmin
    final_mat <- rowSums(log(final_array), dims = 2) / dim(final_array)[[3]]
  }else{
    bin_arr <- final_array == 0
    final_array[bin_arr] <- 1
    final_mat <- rowSums(log(final_array), dims = 2) / rowSums(!bin_arr, dims = 2)
  }

  colnames(final_mat) <- time_bws
  rownames(final_mat) <- net_bws

  return(final_mat)
}




#' @title Bandwidth selection by likelihood cross validation for temporal NKDE (multicore)
#'
#' @description Calculate for multiple network and time bandwidths the cross validation likelihood to
#' select an appropriate bandwidth in a data-driven approach with multicore support
#'
#' @details See the function bws_tnkde_cv_likelihood_calc for more details. Note that the calculation is split
#' according to the grid_shape argument. If the grid_shape is c(1,1) then only one process can be used.
#'
#' @template bw_tnkde_selection-args
#' @template nkde_params-arg
#' @template diggle_corr-arg
#' @template nkde_geoms-args
#' @param time_field The name of the field in events indicating when the events occurred. It must be a numeric field
#' @template sparse-arg
#' @template grid_shape-arg
#' @param sub_sample A float between 0 and 1 indicating the percentage of quadra
#' to keep in the calculus. For large datasets, it may be useful to limit the
#' bandwidth evaluation and thus reduce calculation time.
#' @template verbose-arg
#' @template check-arg
#' @param zero_strat A string indicating what to do when density is 0 when calculating LOO density estimate for an isolated event.
#' "min_double" (default) replace the 0 value by the minimum double possible on the machine. "remove" will remove them from the final
#' score. The first approach penalizes more strongly the small bandwidths.
#' @return A matrix with the cross validation score.  Each row corresponds to a network
#' bandwidth and each column to a time bandwidth (the higher the better).
#' @export
#' @examples
#' \donttest{
#' # loading the data
#' data(mtl_network)
#' data(bike_accidents)
#'
#' # converting the Date field to a numeric field (counting days)
#' bike_accidents$Time <- as.POSIXct(bike_accidents$Date, format = "%Y/%m/%d")
#' bike_accidents$Time <- difftime(bike_accidents$Time, min(bike_accidents$Time), units = "days")
#' bike_accidents$Time <- as.numeric(bike_accidents$Time)
#' bike_accidents <- subset(bike_accidents, bike_accidents$Time>=89)
#'
#' future::plan(future::multisession(workers=1))
#'
#' # calculating the cross validation values
#' cv_scores <- bws_tnkde_cv_likelihood_calc.mc(
#'   bw_net_range = c(100,1000),
#'   bw_net_step = 100,
#'   bw_time_range = c(10,60),
#'   bw_time_step = 5,
#'   lines = mtl_network,
#'   events = bike_accidents,
#'   time_field = "Time",
#'   w = rep(1, nrow(bike_accidents)),
#'   kernel_name = "quartic",
#'   method = "discontinuous",
#'   diggle_correction = FALSE,
#'   study_area = NULL,
#'   max_depth = 10,
#'   digits = 2,
#'   tol = 0.1,
#'   agg = 15,
#'   sparse=TRUE,
#'   grid_shape=c(1,1),
#'   sub_sample=1,
#'   verbose = FALSE,
#'   check = TRUE)
#'
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'}
bws_tnkde_cv_likelihood_calc.mc <- function(bw_net_range, bw_net_step,
                                         bw_time_range, bw_time_step,
                                         lines, events, time_field,
                                         w, kernel_name, method,
                                         diggle_correction = FALSE, study_area = NULL,
                                         max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE,
                                         zero_strat = "min_double",
                                         grid_shape=c(1,1), sub_sample=1, verbose=TRUE, check=TRUE){

  ## step0 basic checks
  samples <- events
  if(time_field %in% names(events) == FALSE){
    stop("time_field must be the name of a numeric column in events")
  }
  events$time <- events[[time_field]]
  events$weight <- w
  div <- "bw"
  events$wid <- 1:nrow(events)

  ## checking inputs
  bw_checks(check, lines, samples, events,
            kernel_name, method, bw_net_range = bw_net_range, bw_time_range = bw_time_range,
            bw_net_step = bw_net_step, bw_time_step = bw_time_step,
            diggle_correction = diggle_correction, study_area = study_area)

  if(zero_strat %in% c("min_double", "remove") == FALSE){
    stop("zero_strat argument must be one of c('min_double', 'remove')")
  }

  ## step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
  }

  data <- prepare_data(samples, lines, events, w , digits,tol,agg)
  lines <- data$lines
  events_loc <- data$events

  #idx <- FNN::knnx.index(st_coordinates(events_loc),st_coordinates(events), k = 1)
  idx <- closest_points(events, events_loc)
  events$goid <- events_loc$goid[idx]

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))

  ## calculating the correction factor for each bw combinaisons
  net_bws <- seq(from = bw_net_range[[1]], to = bw_net_range[[2]], by = bw_net_step)
  time_bws <- seq(from = bw_time_range[[1]], to = bw_time_range[[2]], by = bw_time_step)

  if(verbose){
    print("Calculating the correction factor if required")
  }

  ## calculating network corr_factors
  corr_factors <- bw_tnkde_corr_factor(net_bws, time_bws, diggle_correction, study_area, events, events_loc, lines,
                                       method, kernel_name, tol, digits, max_depth, sparse)
  net_bws_corr <- corr_factors[[1]]
  time_bws_corr <- corr_factors[[2]]


  ## NB : the weights can change because of the different BW in time and space
  ## The weights will be passed to c++, to is must be in an easy format, like an array
  ## event_weights(rows = bws_net, cols = bws_time, slices = events)
  events_weight <- array(0, dim = c(length(net_bws), length(time_bws), nrow(events)))

  for (i in 1:length(time_bws_corr)){
    bw_time <- time_bws[[i]]
    corr_time <- time_bws_corr[[i]]
    for (j in 1:length(net_bws_corr)){
      bw_net <- net_bws[[j]]
      corr_net <- net_bws_corr[[j]]
      events_weight[j,i,] <- events$weight * corr_time * corr_net
    }
  }

  max_bw_net <- max(net_bws)

  ## step3 splitting the dataset with each rectangle
  # NB : here we select the events in the gris (samples) and the events locations in the buffer (events_loc)
  if(verbose){
    print("splitting the data by the grid...")
  }
  selections <- split_by_grid.mc(grid, events, events_loc, lines,max_bw_net, tol, digits, split_all = FALSE)

  ## sub sampling the quadra if required
  if (sub_sample < 1){
    nb <- ceiling(length(selections) * sub_sample)
    selections <- selections[sample(1:length(selections),size = nb,replace = FALSE)]
  }

  ## step 4 calculating the CV values
  if(verbose){
    print("start calculating the CV values ...")
  }

  n_quadra <- length(selections)


  if(verbose){
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {

        # the events_loc must cover the quadra and the bw
        sel_events_loc <- sel$events

        # idem for all the events
        sel_events <- subset(events, events$goid %in% sel_events_loc$goid)

        # but I also need to know on which events I must calculate the densities (in the quadra)
        quad_events <- sel$samples
        sel_weights <- events_weight[,,sel_events$wid]

        values <- spNetwork::tnkde_worker_bw_sel(
                                      sel$lines,
                                      quad_events,
                                      sel_events_loc,
                                      sel_events,
                                      sel_weights,
                                      kernel_name, net_bws, time_bws,
                                      method, div, digits,
                                      tol,sparse, max_depth, verbose)

        p(sprintf("i=%g", sel$index))

        return(values)

      })
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {

      # the events_loc must cover the quadra and the bw
      sel_events_loc <- sel$events

      # idem for all the events
      sel_events <- subset(events, events$goid %in% sel_events_loc$goid)

      # but I also need to know on which events I must calculate the densities (in the quadra)
      quad_events <- sel$samples
      sel_weights <- events_weight[,,sel_events$wid]

      values <- spNetwork::tnkde_worker_bw_sel(sel$lines, quad_events, sel_events_loc, sel_events, sel_weights,
                                    kernel_name, net_bws, time_bws,
                                    method, div, digits,
                                    tol,sparse, max_depth, verbose)
      return(values)
    })
  }

  # removing NULL elements in list of cubes
  dfs[sapply(dfs, is.null)] <- NULL
  dfs$along <- 3

  final_array <- do.call(abind::abind, dfs)
  if(zero_strat == "min_double"){
    bin_arr <- final_array == 0
    final_array[bin_arr] <- .Machine$double.xmin
    final_mat <- rowSums(log(final_array), dims = 2) / dim(final_array)[[3]]
  }else{
    bin_arr <- final_array == 0
    final_array[bin_arr] <- 1
    final_mat <- rowSums(log(final_array), dims = 2) / rowSums(!bin_arr, dims = 2)
  }

  # add <- function(x) Reduce("+", x)
  #
  # final_mat <- add(dfs)
  colnames(final_mat) <- time_bws
  rownames(final_mat) <- net_bws

  return(final_mat)
}





#' @title Worker function fo Bandwidth selection by likelihood cross validation for temporal NKDE
#'
#' @description Calculate for multiple network and time bandwidths the cross validation likelihood to
#' select an appropriate bandwidth in a data-driven approach (INTERNAL)
#' @param lines A feature collection of linestrings representing the underlying network
#' @param quad_events a feature collection of points indicating for which events the densities must be calculated
#' @param events_loc A feature collection of points representing the location of the events
#' @param events A feature collection of points representing the events. Multiple events can share
#' the same location. They are linked by the goid column
#' @param w A numeric array with the weight of the events for each pair of bandwidth
#' @param kernel_name The name of the kernel to use (string)
#' @param bws_net A numeric vector with the network bandwidths
#' @param bws_time A numeric vector with the time bandwidths
#' @param div The type of divisor (not used currently)
#' @param max_depth The maximum depth of recursion
#' @param method The type of NKDE to use (string)
#' @param digits The number of digits to retain from the spatial coordinates. It
#'   ensures that topology is good when building the network. Default is 3. Too high a
#'   precision (high number of digits) might break some connections
#' @param tol A float indicating the minimum distance between the events and the
#'   lines' extremities when adding the point to the network. When points are
#'   closer, they are added at the extremity of the lines.
#' @template sparse-arg
#' @param verbose A boolean
#' @param cvl A boolean indicating if the cvl method (TRUE) or the loo (FALSE) method must be used
#' @return An array with the CV score for each pair of bandiwdths (rows and lines) for each event (slices)
#' @export
#' @examples
#' # no example provided, this is an internal function
tnkde_worker_bw_sel <- function(lines, quad_events, events_loc, events, w, kernel_name, bws_net, bws_time, method, div, digits, tol, sparse, max_depth, verbose = FALSE, cvl = FALSE){

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

  kernel_values <- tnkde_get_loo_values(method,
                                        neighbour_list,
                                        quad_events2$vertex_id,
                                        quad_events2$wid,
                                        quad_events2$time,
                                        events2$vertex_id,
                                        events2$wid, events2$time, w,
                                        bws_net, bws_time,
                                        kernel_name, graph_result$linelist, max_depth,
                                        .Machine$double.xmin
  )

  # kernel_values is supposed to be an array (bw_net, bw_time, events)
  return(kernel_values)
}
