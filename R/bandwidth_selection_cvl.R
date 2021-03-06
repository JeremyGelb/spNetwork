#' @title Bandwidth selection by Cronie and Van Lieshout's Criterion
#'
#' @description Calculate for multiple bandiwdths the Cronie and Van Lieshout's Criterion to
#' select an appropriate bandwidth in a data driven approach.
#'
#' @details The Cronie and Van Lieshout's Criterion \insertCite{cronie2018non}{spNetwork}
#' find the optimal bandwidth by minimizing the difference between the size of the observation
#' window and the sum of the reciprocal of the estimated kernel density at the
#' events locations. In the network case, the size of the study area is the sum
#' of the length of each line in the network. Thus, it is important to only
#' use the necessary parts of the network.
#'
#' @references{
#'     \insertAllCited{}
#' }
#'
#' @param bw_range The range of the bandwidths to consider, given as a numeric
#' vector of two values: c(bandwidth_min, bandwidth_max)
#' @param bw_step The step between each bandwidth to calculate given as a float
#' @param lines A SpatialLinesDataFrame representing the underlying network. The
#' geometries must be a SpatialLinesDataFrame (may crash if some geometries
#'  are invalid)
#' @param events events A SpatialPointsDataFrame representing the events on the
#' network. The points will be snapped on the network.
#' @param w A vector representing the weight of each event
#' @param kernel_name The name of the kernel to use. Must be one of triangle,
#' gaussian, tricube, cosine ,triweight, quartic, epanechnikov or uniform.
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param diggle_correction A Boolean indicating if the correction factor
#' for edge effect must be used.
#' @param study_area A SpatialPolygonsDataFrame or a SpatialPolygon
#' representing the limits of the study area.
#' @param max_depth when using the continuous and discontinuous methods, the
#' calculation time and memory use can go wild  if the network has many
#' small edges (area with many of intersections and many events). To
#' avoid it, it is possible to set here a maximum depth. Considering that the
#' kernel is divided at intersections, a value of 10 should yield good
#' estimates in most cases. A larger value can be used without problem for the
#' discontinuous method. For the continuous method, a larger value will
#' strongly impact calculation speed.
#' @param digits The number of digits to retain in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines' extremities. When
#' points are closer, they are added at the extremity of the lines.
#' @param agg A double indicating if the events must be aggregated within a distance.
#' If NULL, the events are aggregated by rounding the coordinates.
#' @param sparse A Boolean indicating if sparse or regular matrix should be
#' used by the Rcpp functions. Regular matrix are faster, but require more
#' memory and could lead to error, in particular with multiprocessing. Sparse
#' matrix are slower, but require much less memory (not used for the moment).
#' @param grid_shape A vector of two values indicating how the study area
#' must be split when performing the calculus (see details). Default is c(1,1).
#' @param sub_sample A float between 0 and 1 indicating the percentage of quadra
#' to keep in the calculus. For large datasets, it may be useful to limit the
#' bandwidth evaluation and thus reduce calculation time.
#' @param verbose A Boolean, indicating if the function should print messages
#' about process.
#' @param check A Boolean indicating if the geometry checks must be run before
#' calculating the density.
#' @return A dataframe with two columns, one for the bandwidths and the second for
#' the Cronie and Van Lieshout's Criterion.
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
#' bike_accidents <- rgdal::readOGR(eventsgpkg,layer="bike_accidents", verbose=FALSE)
#' cv_scores <- bw_cvl_calc(c(200,400),50,
#'                                mtl_network, bike_accidents,
#'                                rep(1,nrow(bike_accidents)),
#'                                "quartic", "discontinuous",
#'                                diggle_correction = FALSE, study_area = NULL,
#'                                max_depth = 8,
#'                                digits=2, tol=0.1, agg=5,
#'                                sparse=TRUE, grid_shape=c(1,1),
#'                                sub_sample = 1, verbose=TRUE, check=TRUE)
#' }
bw_cvl_calc <- function(bw_range,bw_step,lines, events, w, kernel_name, method,  diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE, grid_shape=c(1,1), sub_sample=1, verbose=TRUE, check=TRUE){

  ## step0 basic checks
  samples <- events
  div <- "bw"

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

  if(min(bw_range)<=0){
    stop("the bandwidth for the kernel must be superior to 0")
  }

  if(bw_step<=0){
    stop("the step between two bandwidths must be greater than 0")
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
  Wl <- rgeos::gLength(lines) / 1000
  samples <- data$events
  events <- data$events

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))

  ## calculating the correction factor for each bw
  all_bws <- seq(min(bw_range),max(bw_range),bw_step)
  if(verbose){
    print("Calculating the correction factor if required")
  }
  for (bw in all_bws){
    if(diggle_correction){
      bws <- rep(bw,nrow(events))
      corr_factor <- correction_factor(study_area,events,lines,method,bws, kernel_name, tol, digits, max_depth, sparse)
    }else{
      corr_factor<- rep(1,nrow(events))
    }

    events[[paste("weight_",bw,sep="")]] <- events$weight * corr_factor
    samples[[paste("weight_",bw,sep="")]] <- samples$weight * corr_factor

  }
  max_bw <- max(bw_range)

  ## step3 splitting the dataset with each rectangle
  selections <- split_by_grid(grid,samples,events,lines,max_bw, tol, digits, split_all = FALSE)

  ## sub sampling the quadra if required
  if (sub_sample < 1){
    nb <- ceiling(length(selections) * sub_sample)
    selections <- selections[sample(1:length(selections),size = nb,replace = F)]
  }

  ## step 4 calculating the CVl values
  if(verbose){
    print("start calculating the CVl values ...")
  }

  n_quadra <- length(selections)

  if (verbose){
    pb <- txtProgressBar(min = 0, max = n_quadra, style = 3)
  }
  dfs <- lapply(1:n_quadra,function(i){
    sel <- selections[[i]]

    values <- nkde_worker_bw_sel_cvl(sel$lines, sel$events,
                                 sel$samples, kernel_name, all_bws,
                                 method, div, digits,
                                 tol,sparse, max_depth, verbose)

    if(verbose){
      setTxtProgressBar(pb, i)
    }

    ## on combine les resultats de chaque bw dans une matrice
    df <- data.frame(do.call(cbind,values))
    names(df) <- paste("k",1:ncol(df),sep="")
    df$goid <- sel$samples$goid
    return(df)
  })

  ## step5  combining the results for each quadra
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]


  ## calculer les valeurs de CVl
  cv_scores <- sapply(1:length(all_bws), function(i){
    bw <- all_bws[[i]]
    kvalues <- tot_df[,i]
    score <- (Wl - sum(1/kvalues))**2
    return(score)
  })

  finaldf <- data.frame(
    "bw" = all_bws,
    "cvl_scores" = cv_scores
  )

  return(finaldf)
}


#' @title Bandwidth selection by Cronie and Van Lieshout's Criterion (multicore version)
#'
#' @description Calculate for multiple bandiwdths the Cronie and Van Lieshout's Criterion to
#' select an appropriate bandwidth in a data driven approach. A future plan can be used
#' to split the work across cores. The different quadras generated in accordance with the
#' argument grid_shape are used for the parallelisation. So if only one quadra is
#' generate (grid_shape = c(1,1)), the function will still use only one core. This also
#' affect the progress bar.
#'
#'
#' @details For more details, see help(bw_cvl_calc)
#'
#' @param bw_range The range of the bandwidths to consider, given as a numeric
#' vector of two values: c(bandwidth_min, bandwidth_max)
#' @param bw_step The step between each bandwidth to calculate given as a float
#' @param lines A SpatialLinesDataFrame representing the underlying network. The
#' geometries must be a SpatialLinesDataFrame (may crash if some geometries
#'  are invalid)
#' @param events events A SpatialPointsDataFrame representing the events on the
#' network. The points will be snapped on the network.
#' @param w A vector representing the weight of each event
#' @param kernel_name The name of the kernel to use. Must be one of triangle,
#' gaussian, tricube, cosine ,triweight, quartic, epanechnikov or uniform.
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param diggle_correction A Boolean indicating if the correction factor
#' for edge effect must be used.
#' @param study_area A SpatialPolygonsDataFrame or a SpatialPolygon
#' representing the limits of the study area.
#' @param max_depth when using the continuous and discontinuous methods, the
#' calculation time and memory use can go wild  if the network has many
#' small edges (area with many of intersections and many events). To
#' avoid it, it is possible to set here a maximum depth. Considering that the
#' kernel is divided at intersections, a value of 10 should yield good
#' estimates in most cases. A larger value can be used without problem for the
#' discontinuous method. For the continuous method, a larger value will
#' strongly impact calculation speed.
#' @param digits The number of digits to retain in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines' extremities. When
#' points are closer, they are added at the extremity of the lines.
#' @param agg A double indicating if the events must be aggregated within a distance.
#' If NULL, the events are aggregated by rounding the coordinates.
#' @param sparse A Boolean indicating if sparse or regular matrix should be
#' used by the Rcpp functions. Regular matrix are faster, but require more
#' memory and could lead to error, in particular with multiprocessing. Sparse
#' matrix are slower, but require much less memory (not used for the moment).
#' @param grid_shape A vector of two values indicating how the study area
#' must be split when performing the calculus (see details). Default is c(1,1).
#' @param sub_sample A float between 0 and 1 indicating the percentage of quadra
#' to keep in the calculus. For large datasets, it may be useful to limit the
#' bandwidth evaluation and thus reduce calculation time.
#' @param verbose A Boolean, indicating if the function should print messages
#' about process.
#' @param check A Boolean indicating if the geometry checks must be run before
#' calculating the density.
#' @return A dataframe with two columns, one for the bandwidths and the second for
#' the Cronie and Van Lieshout's Criterion.
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
#' bike_accidents <- rgdal::readOGR(eventsgpkg,layer="bike_accidents", verbose=FALSE)
#' future::plan(future::multisession(workers=2))
#' cv_scores <- bw_cvl_calc(c(200,400),50,
#'                                mtl_network, bike_accidents,
#'                                rep(1,nrow(bike_accidents)),
#'                                "quartic", "discontinuous",
#'                                diggle_correction = FALSE, study_area = NULL,
#'                                max_depth = 8,
#'                                digits=2, tol=0.1, agg=5,
#'                                sparse=TRUE, grid_shape=c(1,1),
#'                                sub_sample = 1, verbose=TRUE, check=TRUE)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
bw_cvl_calc.mc <- function(bw_range,bw_step,lines, events, w, kernel_name, method,  diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE, grid_shape=c(1,1), sub_sample=1, verbose=TRUE, check=TRUE){

  ## step0 basic checks
  samples <- events
  div <- "bw"

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

  if(min(bw_range)<=0){
    stop("the bandwidth for the kernel must be superior to 0")
  }

  if(bw_step<=0){
    stop("the step between two bandwidths must be greater than 0")
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
  Wl <- rgeos::gLength(lines) / 1000
  samples <- data$events
  events <- data$events

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))
  all_bws <- seq(min(bw_range),max(bw_range),bw_step)

  ## calculating the correction factor for each bw
  if(verbose){
    print("Calculating the correction factor if required")
  }
  for (bw in all_bws){
    if(diggle_correction){
      bws <- rep(bw,nrow(events))
      corr_factor <- correction_factor(study_area,events,lines,method,bws, kernel_name, tol, digits, max_depth, sparse)
    }else{
      corr_factor<- rep(1,nrow(events))
    }

    events[[paste("weight_",bw,sep="")]] <- events$weight * corr_factor
    samples[[paste("weight_",bw,sep="")]] <- samples$weight * corr_factor

  }

  max_bw <- max(bw_range)

  ## step3 splitting the dataset with each rectangle
  selections <- split_by_grid.mc(grid,samples,events,lines,max_bw, tol, digits, split_all = FALSE)

  ## sub sampling the quadra if required
  if (sub_sample < 1){
    nb <- ceiling(length(selections) * sub_sample)
    selections <- selections[sample(1:length(selections),size = nb,replace = F)]
  }

  ## step 4 calculating the CVl values
  if(verbose){
    print("start calculating the CVl values ...")
  }

  n_quadra <- length(selections)

  if(verbose){
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {

        values <- nkde_worker_bw_sel(sel$lines, sel$events,
                                     sel$samples, kernel_name, all_bws,
                                     method, div, digits,
                                     tol,sparse, max_depth, verbose)
        p(sprintf("i=%g", sel$index))
        ## on combine les resultats de chaque bw dans une matrice
        df <- data.frame(do.call(cbind,values))
        names(df) <- paste("k",1:ncol(df),sep="")
        df$goid <- sel$samples$goid
        return(df)

      })
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {

      values <- nkde_worker_bw_sel(sel$lines, sel$events,
                                   sel$samples, kernel_name, all_bws,
                                   method, div, digits,
                                   tol,sparse, max_depth, verbose)

      df <- data.frame(do.call(cbind,values))
      names(df) <- paste("k",1:ncol(df),sep="")
      df$goid <- sel$samples$goid
      return(df)

    })
  }


  ## step5  combining the results for each quadra
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]

  ## calculer les valeurs de CVl
  cv_scores <- sapply(1:length(all_bws), function(i){
    bw <- all_bws[[i]]
    kvalues <- tot_df[,i]
    score <- (Wl - sum(1/kvalues))**2
    return(score)
  })

  finaldf <- data.frame(
    "bw" = all_bws,
    "cvl_scores" = cv_scores
  )

  return(finaldf)
}





#' @title Worker function for bandwidth selection by Cronie and Van Lieshout's Criterion
#'
#' @description The worker function for bandwidth selection by Cronie and Van Lieshout's Criterion
#'
#' @param lines A SpatialLinesDataFrame representing the underlying network. The
#' geometries must be a SpatialLinesDataFrame (may crash if some geometries
#'  are invalid)
#' @param events A SpatialPointsDataFrame representing the events on the
#' network.
#' @param samples A SpatialPointsDataFrame representing the samples on the
#' network.
#' @param kernel_name The name of the kernel to use. Must be one of triangle,
#' gaussian, tricube, cosine ,triweight, quartic, epanechnikov or uniform.
#' @param bws A vector with all the bandiwdths to test.
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param div The divisor to use (should always be dist here).
#' @param max_depth when using the continuous and discontinuous methods, the
#' calculation time and memory use can go wild  if the network has many
#' small edges (area with many of intersections and many events). To
#' avoid it, it is possible to set here a maximum depth. Considering that the
#' kernel is divided at intersections, a value of 10 should yield good
#' estimates in most cases. A larger value can be used without problem for the
#' discontinuous method. For the continuous method, a larger value will
#' strongly impact calculation speed.
#' @param digits The number of digits to retain in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines' extremities. When
#' points are closer, they are added at the extremity of the lines.
#' @param agg A double indicating if the events must be aggregated within a distance.
#' If NULL, the events are aggregated by rounding the coordinates.
#' @param sparse A Boolean indicating if sparse or regular matrix should be
#' used by the Rcpp functions. Regular matrix are faster, but require more
#' memory and could lead to error, in particular with multiprocessing. Sparse
#' matrix are slower, but require much less memory (not used for the moment).
#' @param verbose A Boolean, indicating if the function should print messages
#' about process.
#' @return A list of dataframes (continuous kernel) or a list of numeric vectors (other kernels).
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
nkde_worker_bw_sel_cvl <- function(lines, events, samples, kernel_name, bws, method, div, digits, tol, sparse, max_depth, verbose = FALSE){

  # if we do not have event in that space, just return 0 values
  if(nrow(events)==0){
    values <- lapply(bws,function(i){rep(0,nrow(samples))})
    return(values)
  }

  ## step1 creating the graph
  graph_result <- build_graph(lines,digits = digits,line_weight = "length")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  edges <- graph_result$spedges

  ## step2 finding for each event, its corresponding node
  ## NOTE : there will be less samples than events most of the time
  ## because of the avoidance of island effects.
  events$vertex_id <- closest_points(events, nodes)
  samples$vertex_id <- closest_points(samples, nodes)
  events$oid <- 1:nrow(events)

  ## step3 starting the calculations !
  neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")
  neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})

  ## we calculate the nkde values for each bw provided
  kernel_values <- lapply(bws, function(bw){

    repbws <- rep(bw,nrow(events))

    if(method == "simple"){
      values <- spNetwork::get_loo_values_simple(neighbour_list, samples$vertex_id, samples[[paste("weight_",bw,sep="")]],
                                                 events$vertex_id, events[[paste("weight_",bw,sep="")]],
                                                 repbws, kernel_name, graph_result$linelist, max_depth)
    }else if(method=="continuous"){
      ##and finally calculating the values
      # NOTE : Values is a dataframe with two columns :
      # the kvalue calculated for each event when all the events are considered
      # the kvalue specific at each event (when all other events are not considered)
      values <- spNetwork::get_loo_values_continuous(neighbour_list,  samples$vertex_id, samples[[paste("weight_",bw,sep="")]],
                                                     events$vertex_id, events[[paste("weight_",bw,sep="")]],
                                                     repbws, kernel_name, graph_result$linelist, max_depth)
      values <- values$sum_k

    }else if(method == "discontinuous"){
      values <- spNetwork::get_loo_values_discontinuous(neighbour_list, samples$vertex_id, samples[[paste("weight_",bw,sep="")]],
                                                        events$vertex_id, events[[paste("weight_",bw,sep="")]], repbws,
                                                        kernel_name, graph_result$linelist, max_depth)
    }

    ## step7 adjusting the kernel values !
    # dividing by bw is crucial, other wise, larger BW are always better !
    return(values * (1/bw))
  })

  ## at that point, we have a list of numeric vectors for each bw
  return(kernel_values)
}


