# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### helpers for k-functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Preparing results for K functions
#'
#' @description Prepare the final results at the end of the execution of the main
#' functions calculating K or G functions.
#'
#' @param k_vals a numeric vector with the real K values
#' @param g_vals a numeric vector with the real g values
#' @param all_values a list with the simulated K and G values that must be arranged.
#' @param conf_int the confidence interval parameter.
#' @param calc_g_func a boolean indicating if the G function has been calculated.
#' @param cross a boolean indicating if we have calculated a simple (FALSE) or a cross function.
#' @param dist_seq a numeric vector representing the distance used for calculation
#' @param return_sims a boolean, indicating if the simulations must be returned
#'
#' @return A list with the following values :
#'
#' * plotk: A ggplot2 object representing the values of the k-function
#' * plotg: A ggplot2 object representing the values of the g-function
#' * values: A DataFrame with the values used to build the plots
#'
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot geom_ribbon geom_path aes_string labs
#' @md
#' @keywords internal
#' @examples
#' # no example, this is an internal function
prep_kfuncs_results <- function(k_vals, g_vals, all_values, conf_int, calc_g_func, cross, dist_seq, return_sims){

  ## step8 : extract the k_vals and g_vals matrices
  if(calc_g_func){
    k_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,1])}))
    g_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,2])}))
  }else{
    k_mat <- do.call(cbind,all_values)
  }


  ## step9 : calculating the summary stats
  upper <- 1-conf_int / 2
  lower <- conf_int / 2
  k_stats <- apply(k_mat,MARGIN = 1, function(i){
    return(quantile(i,probs = c(lower,upper)))
  })

  plot_df <- data.frame(
    "obs_k" = k_vals,
    "lower_k" = k_stats[1,],
    "upper_k" = k_stats[2,],
    "distances" = dist_seq
  )

  plotk <- ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_k",ymax = "upper_k"),
                fill = grDevices::rgb(0.1,0.1,0.1), alpha=0.4)+
    geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")

  if(cross){
    plotk <- plotk + labs(x = "distances",
                          y = "empirical cross K-function")
  }else{
    plotk <- plotk + labs(x = "distances",
                          y = "empirical K-function")
  }




  obj <- list(
    "plotk" = plotk,
    "values" = plot_df
  )


  if(calc_g_func){
    g_stats <- apply(g_mat,MARGIN = 1, function(i){
      return(quantile(i,probs = c(lower,upper)))
    })

    plot_df$obs_g <- g_vals
    plot_df$lower_g <- g_stats[1,]
    plot_df$upper_g <- g_stats[2,]

    plotg <- ggplot(plot_df)+
      geom_ribbon(aes_string(x = "distances", ymin="lower_g", ymax = "upper_g"),
                  fill = grDevices::rgb(0.1,0.1,0.1),alpha=0.4, )+
      geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
                linetype="dashed")+
      geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
                linetype="dashed")+
      geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")

    if(cross){
      plotg <- plotg + labs(x = "distances",
                            y = "empirical cross G-function")
    }else{
      plotg <- plotg + labs(x = "distances",
                            y = "empirical G-function")
    }



    obj$plotg <- plotg
    obj$values <- plot_df

  }


  if(return_sims){
    obj$sim_k_values <- k_mat
    if(calc_g_func){
      obj$sim_g_values <- g_mat
    }
  }
  return(obj)
}




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### execution k functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Network k and g functions (maturing)
#'
#' @description Calculate the k and g functions for a set of points on a
#'   network (maturing).
#'
#' @details The k-function is a method to characterize the dispersion of a set
#'   of points. For each point, the numbers of other points in subsequent radii
#'   are calculated. This empirical k-function can be more or less clustered
#'   than a k-function obtained if the points were randomly located in space. In
#'   a network, the network distance is used instead of the Euclidean distance.
#'   This function uses Monte Carlo simulations to assess if the points are
#'   clustered or dispersed, and gives the results as a line plot. If the line
#'   of the observed k-function is higher than the shaded area representing the
#'   values of the simulations, then the points are more clustered than what we
#'   can expect from randomness and vice-versa. The function also calculates the
#'   g-function, a modified version of the k-function using rings instead of
#'   disks. The width of the ring must be chosen. The main interest is to avoid
#'   the cumulative effect of the classical k-function. This function is maturing,
#'   it works as expected (unit tests) but will probably be modified in the
#'   future releases (gain speed, advanced features, etc.).
#'
#' @template kfunctions-arg
#' @template common_kfunctions-arg
#' @param return_sims a boolean indicating if the simulated k and g values must also
#' be returned as matrices
#'
#' @return A list with the following values :
#'
#' * plotk: A ggplot2 object representing the values of the k-function
#' * plotg: A ggplot2 object representing the values of the g-function
#' * values: A DataFrame with the values used to build the plots
#'
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot geom_ribbon geom_path aes_string labs
#' @export
#' @md
#' @examples
#' \donttest{
#' data(main_network_mtl)
#' data(mtl_libraries)
#' result <- kfunctions(main_network_mtl, mtl_libraries,
#'      start = 0, end = 2500, step = 10,
#'      width = 200, nsim = 50,
#'      conf_int = 0.05, tol = 0.1, agg = NULL,
#'      verbose = FALSE)
#' }
new_kfunctions <- function(lines, points,
                       start,
                       end,
                       step,
                       width,
                       nsim,
                       conf_int = 0.05,
                       digits = 2,
                       tol = 0.1,
                       resolution = NULL,
                       agg = NULL,
                       verbose = TRUE,
                       return_sims = FALSE,
                       calc_g_func = TRUE
                       ){

  ## step0 : clean the points
  if (verbose){
    print("Preparing data ...")
  }
  n <- nrow(points)
  points$goid <- seq_len(nrow(points))
  points$weight <- rep(1,nrow(points))
  points <- clean_events(points,digits,agg)

  probs <- NULL

  ## step1 : clean the lines
  if(is.null(probs)){
    lines$probs <- 1
  }else{
    lines$probs <- probs
  }

  lines$length <- as.numeric(st_length(lines))
  lines <- subset(lines, lines$length>0)
  lines$oid <- seq_len(nrow(lines))


  ## step2 : adding the points to the lines
  if (verbose){
    print("Snapping points on lines ...")
  }
  snapped_events <- snapPointsToLines2(points,lines,idField = "oid")
  new_lines <- split_lines_at_vertex(lines, snapped_events,
                                     snapped_events$nearest_line_id, tol)

  # new_lines <- add_vertices_lines(lines,snapped_events,
  #                                 snapped_events$nearest_line_id, tol)

  ## step3 : splitting the lines
  if (verbose){
    print("Building graph ...")
  }
  #new_lines <- simple_lines(new_lines)
  new_lines$length <- as.numeric(st_length(new_lines))
  new_lines <- subset(new_lines,new_lines$length>0)

  new_lines <- remove_loop_lines(new_lines,digits)

  new_lines$oid <- seq_len(nrow(new_lines))
  new_lines <- new_lines[c("length","oid","probs")]
  Lt <- sum(as.numeric(st_length(new_lines)))

  new_lines$weight <- as.numeric(st_length(new_lines))

  ## step4 : building the graph for the real case
  graph_result <- build_graph_cppr(new_lines,digits = digits,
                              line_weight = "weight",
                              attrs = TRUE, direction = NULL)

  graph <- graph_result$graph
  nodes <- graph_result$spvertices


  node_id <- closest_points(snapped_events, nodes)
  snapped_events$vertex_id <- nodes$ref[node_id]

  ## step5 : calculating the distance matrix
  dist_mat <- cppRouting::get_distance_matrix(graph, from = snapped_events$vertex_id, to = snapped_events$vertex_id)

  ## step6 : calcualte the kfunction and the g function
  if (verbose){
    print("Calculating k and g functions ...")
  }

  #diag(dist_mat) <- end+width+1

  if(calc_g_func){
    # if required, we also calculate the G function
    k_g_vals <- kgfunc_cpp2(dist_mat = dist_mat,
                            start = start,
                            end = end,
                            step = step,
                            widt = width,
                            Lt = Lt,
                            n = n,
                            wc = snapped_events$weight,
                            wr = snapped_events$weight)
    k_vals <- k_g_vals[,1]
    g_vals <- k_g_vals[,2]
  }else{
    # otherwise, only the K function
    k_vals <- kfunc_cpp2(dist_mat,start,end,step,Lt,n,wc = snapped_events$weight, wr = snapped_events$weight)
  }


  ## step7 : generate the permutations
  if (verbose){
    print("Calculating the simulations ...")
  }

  if(verbose){
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }

  # the case where we can must edit the network and create small chunks
  if (is.null(resolution)==FALSE){

    # we start by creating the lixels
    new_lines2 <- lixelize_lines(new_lines, resolution, mindist = resolution / 2.0)
    new_lines2$weight <- as.numeric(st_length(new_lines2))


    # we can now rebuild the network with these new lines
    graph_result <- build_graph_cppr(new_lines2,digits = digits,
                                     line_weight = "weight",
                                     attrs = TRUE, direction = NULL)

    graph <- graph_result$graph
    nodes <- graph_result$spvertices


    node_id <- closest_points(snapped_events, nodes)
    snapped_events$vertex_id <- nodes$id[node_id]

  }

  # we can now do the simulation by permutating the location of the events on
  # the network.

  all_values <- lapply(1:nsim, function(i){

    vidx <- sample(graph$dict$ref, size = sum(snapped_events$weight))
    dist_mat <- cppRouting::get_distance_matrix(graph, from = vidx, to = vidx)

    #diag(dist_mat) <- end+width+1

    if(calc_g_func){
      k_g_vals <- kgfunc_cpp2(dist_mat,start,end,step,width,Lt,n,wc = snapped_events$weight, wr = snapped_events$weight)
      return(k_g_vals)
    }else{
      k_vals <- kfunc_cpp2(dist_mat,start,end,step,Lt,n,wc = snapped_events$weight, wr = snapped_events$weight)
      return(k_vals)
    }

  })


  ## Then we can prepare the results
  obj <- prep_kfuncs_results(k_vals, g_vals, all_values, conf_int, calc_g_func,
                             cross = FALSE, dist_seq = seq(start, end, step),
                             return_sims = return_sims)


  return(obj)
}




#' @title Network k and g functions (multicore)
#'
#' @description Calculate the k and g functions for a set of points on a network
#'   with multicore support. For details, please see the function kfunctions.
#'   (maturing)
#'
#' @details For details, please look at the function kfunctions.
#'
#' @template kfunctions-arg
#' @template common_kfunctions-arg
#' @param return_sims a boolean indicating if the simulated k and g values must also
#' be returned as matrices
#' @template grid_shape-arg
#'
#' @return A list with the following values :
#'
#' * plotk: A ggplot2 object representing the values of the k-function
#' * plotg: A ggplot2 object representing the values of the g-function
#' * values: A DataFrame with the values used to build the plots
#'
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot geom_ribbon geom_path aes_string labs
#' @export
#' @md
#' @examples
#' \donttest{
#' data(main_network_mtl)
#' data(mtl_libraries)
#' result <- kfunctions(main_network_mtl, mtl_libraries,
#'      start = 0, end = 2500, step = 10,
#'      width = 200, nsim = 50,
#'      conf_int = 0.05, tol = 0.1, agg = NULL,
#'      verbose = FALSE)
#' }
new_kfunctions.mc <- function(lines, points,
                           start,
                           end,
                           step,
                           width,
                           nsim,
                           conf_int = 0.05,
                           digits = 2,
                           tol = 0.1,
                           resolution = NULL,
                           agg = NULL,
                           verbose = TRUE,
                           return_sims = FALSE,
                           calc_g_func = TRUE,
                           grid_shape = c(1,1)){

  ## step0 : clean the points
  if (verbose){
    print("Preparing data ...")
  }
  n <- nrow(points)
  points$goid <- seq_len(nrow(points))
  points$weight <- rep(1,nrow(points))
  points <- clean_events(points,digits,agg)


  ## step1 : clean the lines
  lines <- remove_loop_lines(lines,digits)
  lines$length <- as.numeric(st_length(lines))
  lines <- subset(lines, lines$length>0)
  lines$oid <- seq_len(nrow(lines))


  ##______________________
  # Gridding process

  ## we will now define a grid
  grid <- build_grid(grid_shape = grid_shape, spatial = list(lines, points))

  ## we create two spatial indices to make the requests faster
  tree_points <- build_quadtree(points)
  tree_lines <- build_quadtree(lines)


  ## and split my data according to this grid
  selections <- lapply(1:nrow(grid),function(i){
    square <- grid[i,]
    # selecting the points in the grid
    sel_points <- spatial_request(square,tree_points,points)
    # if there is no sampling points in the rectangle, then return NULL
    if(nrow(sel_points)==0){
      return(NULL)
    }
    # selecting the events in a buffer
    bw <- end + width
    buff <- st_buffer(square,dist = bw)

    sel_points2 <- spatial_request(buff,tree_points,points)

    # selecting the lines in a buffer
    sel_lines <- spatial_request(buff,tree_lines,lines)
    sel_lines$oid <- 1:nrow(sel_lines)

    return(list(
      'lines' = sel_lines,
      'origins' = sel_points,
      'destinations' = sel_points2
    ))

  })

  ##______________________
  # preparing values for randomisations
  Lt <- sum(as.numeric(st_length(lines)))
  N <- sum(points$weight)

  selections <- selections[lengths(selections) != 0]


  ##______________________
  # applying the function on the grid

  # in this step, we will perform the network k and g functions on each element
  # of the grid. The idea is to obtain the counting matrices for each square and then
  # to calculate the means at each distance. The randomization is also done locally
  # by sampling the vertices. We must take into account the local length of the network
  # to do so.

  results <- lapply(1:length(selections), function(i){

    # we extract the selections first
    sel <- selections[[i]]
    sub_lines <- sel$lines
    origins <- sel$origins
    destinations <- sel$destinations

    # then, we snapp the points on the lines
    snapped_events <- snapPointsToLines2(destinations, sub_lines, idField = "oid")
    new_lines <- split_lines_at_vertex(sub_lines, snapped_events,
                                       snapped_events$nearest_line_id, tol)

    # we create a graph
    new_lines$length <- as.numeric(st_length(new_lines))
    new_lines <- subset(new_lines,new_lines$length>0)



    new_lines$oid <- seq_len(nrow(new_lines))
    new_lines <- new_lines[c("length","oid")]
    local_Lt <- sum(as.numeric(st_length(new_lines)))

    new_lines$weight <- as.numeric(st_length(new_lines))

    graph_result <- build_graph_cppr(new_lines,digits = digits,
                                     line_weight = "weight",
                                     attrs = TRUE, direction = NULL)

    graph <- graph_result$graph
    nodes <- graph_result$spvertices

    # we must find the the locations of the events on the graph
    node_id <- closest_points(snapped_events, nodes)
    snapped_events$vertex_id <- nodes$ref[node_id]

    start_nodes <- subset(snapped_events, snapped_events$goid %in% origins$goid)

    # and calculate the distance matrix
    dist_mat <- cppRouting::get_distance_matrix(graph, from = start_nodes$vertex_id, to = snapped_events$vertex_id)
    values <- list()
    values$sumw <- sum(start_nodes$weight)

    # with this distance matrix, we can calculate the counting matrix for this square
    if(calc_g_func){
      matrices <- kgfunc_counting(dist_mat,
                                  wc = snapped_events$weight,
                                  wr = start_nodes$weight,
                                  breaks = rev(seq_num3(start, end, step)),
                                  width = width / 2)
      values$k_counts <- matrices[[1]]
      values$g_counts <- matrices[[2]]

    }else{
      values$k_counts <- kfunc_counting(dist_mat,
                                  wc = snapped_events$weight,
                                  wr = start_nodes$weight,
                                  breaks = rev(seq_num3(start, end, step))
                                 )

    }


    ## Excellent ! The next step is to do the permutation locally. To do so, we must first apply the
    ## sampling strategy provided by the user.
    if (is.null(resolution)==FALSE){

      # we start by creating the lixels
      new_lines2 <- lixelize_lines(new_lines, resolution, mindist = resolution / 2.0)
      new_lines2$weight <- as.numeric(st_length(new_lines2))


      # we can now rebuild the network with these new lines
      graph_result <- build_graph_cppr(new_lines2,digits = digits,
                                       line_weight = "weight",
                                       attrs = TRUE, direction = NULL)

      graph <- graph_result$graph
      nodes <- graph_result$spvertices


      node_id <- closest_points(snapped_events, nodes)
      snapped_events$vertex_id <- nodes$id[node_id]

    }

    sample_size <- round((local_Lt / Lt) * N)

    my_breaks <- rev(seq_num3(start, end, step))

    simulations <- lapply(1:nsim, function(j){
      vidx <- sample(graph$dict$ref, size = sample_size)
      dist_mat <- cppRouting::get_distance_matrix(graph, from = vidx, to = vidx)

      values <- list()
      values$sumw <- sample_size

      # with this distance matrix, we can calculate the counting matrix for this square
      if(calc_g_func){
        matrices <- kgfunc_counting(dist_mat,
                                    wc = rep(1, ncol(dist_mat)),
                                    wr = rep(1, nrow(dist_mat)),
                                    breaks = my_breaks,
                                    width = width / 2)
        values$k_counts <- matrices[[1]]
        values$g_counts <- matrices[[2]]

      }else{
        values$k_counts <- kfunc_counting(dist_mat,
                                          wc = rep(1, ncol(dist_mat)),
                                          wr = rep(1, nrow(dist_mat)),
                                          breaks = my_breaks)
      }
      return(values)


    })

    values$simulations <- simulations

    return(values)

  })

  ## its is time to combine the values obtained by the different tiles

  # let me start with the counts obtained for the k function
  k_counts <- do.call(rbind, lapply(results, function(x){x$k_counts}))

  # the countings must be transformed in kvalues. For each distance, it is the mean of
  # the reached points multiplied by t1
  t1 <- 1.0/((N-1)/Lt)
  k_vals <- colSums(k_counts) / N * t1
  k_vals <- rev(k_vals)

  if(calc_g_func){
    g_counts <- do.call(rbind, lapply(results, function(x){x$g_counts}))
    g_vals <- colSums(g_counts) / N * t1
    g_vals <- rev(g_vals)
  }



  # and then the simulations
  # all values must be a list with an element per simulation. Each element will
  # be a matrix of two columns. One with the k values, and the second with the g values

  all_values <- lapply(1:nsim, function(i){
    k_sim_vals <- do.call(rbind, lapply(results, function(x){
      x$simulations[[i]]$k_counts
    }))
    sims_k <- colSums(k_sim_vals) / N * t1
    sims_k <- rev(sims_k)
    if(calc_g_func){
      g_sim_vals <- do.call(rbind, lapply(results, function(x){
        x$simulations[[i]]$g_counts
      }))
      sims_g <- colSums(g_sim_vals) / N * t1
      sims_g <- rev(sims_g)
      return(cbind(sims_k,sims_g))
    }else{
      return(sims_k)
    }
  })


  ## Then we can prepare the results
  obj <- prep_kfuncs_results(k_vals, g_vals, all_values, conf_int, calc_g_func,
                             cross = FALSE, dist_seq = seq(start, end, step),
                             return_sims = return_sims)


  return(obj)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### execution cross k functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Network cross k and g functions (maturing)
#'
#' @description Calculate the cross k and g functions for a set of points on a
#'   network. (maturing)
#'
#' @details The cross k-function is a method to characterize the dispersion of a
#'   set of points (A) around a second set of points (B). For each point in B,
#'   the numbers of other points in A in subsequent radii are calculated. This
#'   empirical cross k-function can be more or less clustered than a cross
#'   k-function obtained if the points in A were randomly located around points
#'   in B. In a network, the network distance is used instead of the Euclidean
#'   distance. This function uses Monte Carlo simulations to assess if the
#'   points are clustered or dispersed and gives the results as a line plot. If
#'   the line of the observed cross k-function is higher than the shaded area
#'   representing the values of the simulations, then the points in A are more
#'   clustered around points in B than what we can expect from randomness and
#'   vice-versa. The function also calculates the cross g-function, a modified
#'   version of the cross k-function using rings instead of disks. The width of
#'   the ring must be chosen. The main interest is to avoid the cumulative
#'   effect of the classical k-function. Note that the cross k-function of
#'   points A around B is not necessarily the same as the cross k-function of
#'   points B around A. This function is maturing, it works as expected (unit
#'   tests) but will probably be modified in the future releases (gain speed,
#'   advanced features, etc.).
#'
#' @template kross_kfunctions-arg
#' @template common_kfunctions-arg
#' @param return_sims a boolean indicating if the simulated k and g values must also
#' be returned as matrices
#'
#'
#' @return A list with the following values : \cr \item{plotk}{ A
#'   ggplot2 object representing the values of the cross k-function}
#'   \item{plotg}{ A ggplot2 object representing the values of the cross
#'   g-function} \item{values}{ A DataFrame with the values used to build the
#'   plots}
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot geom_ribbon geom_path aes_string labs
#' @importFrom grDevices rgb
#' @export
#' @examples
#' \donttest{
#' data(main_network_mtl)
#' data(mtl_libraries)
#' data(mtl_theatres)
#' result <- cross_kfunctions(main_network_mtl, mtl_theatres, mtl_libraries,
#'                            start = 0, end = 2500, step = 10, width = 250,
#'                            nsim = 50, conf_int = 0.05, digits = 2,
#'                            tol = 0.1, agg = NULL, verbose = FALSE)
#' }
new_cross_kfunctions <- function(lines, pointsA, pointsB,
                             start, end, step, width,
                             nsim, conf_int = 0.05,
                             digits = 2, tol = 0.1,
                             resolution = NULL, agg = NULL,
                             verbose = TRUE, return_sims = FALSE, calc_g_func = TRUE){

  ## step0 : clean the points
  if(verbose){
    print("Preparing data ...")
  }

  probs <- NULL

  st_geometry(pointsA) <- 'geometry'
  st_geometry(pointsB) <- 'geometry'

  pointsA$weight <- rep(1,nrow(pointsA))
  pointsA <- clean_events(pointsA,digits,agg)
  pointsA$goid <- seq_len(nrow(pointsA))

  pointsB$weight <- rep(1,nrow(pointsB))
  pointsB <- clean_events(pointsB,digits,agg)
  pointsB$goid <- seq_len(nrow(pointsB))

  na <- sum(pointsA$weight)
  nb <- sum(pointsB$weight)

  ## step1 : clean the lines
  lines$length <- as.numeric(st_length(lines))
  lines <- subset(lines, lines$length>0)
  lines$oid <- seq_len(nrow(lines))

  ## step2 : adding the points to the lines
  if(verbose){
    print("Snapping points on lines ...")
  }
  pointsA$type <- "A"
  pointsB$type <- "B"
  all_events <- rbind(pointsA[c("type","goid","weight")],
                      pointsB[c("type","goid","weight")])

  snapped_events <- snapPointsToLines2(all_events, lines, idField = "oid")


  new_lines <- split_lines_at_vertex(lines, snapped_events,
                                     snapped_events$nearest_line_id, tol)


  ## step3 : splitting the lines
  if(verbose){
    print("Building graph ...")
  }
  #new_lines <- simple_lines(new_lines)
  new_lines$length <- as.numeric(st_length(new_lines))
  new_lines <- subset(new_lines,new_lines$length>0)

  new_lines <- remove_loop_lines(new_lines,digits)

  new_lines$oid <- seq_len(nrow(new_lines))
  new_lines <- new_lines[c("length","oid")]
  Lt <- sum(as.numeric(st_length(new_lines)))

  new_lines$weight <- as.numeric(st_length(new_lines))

  ## step4 : building the graph for the real case

  graph_result <- build_graph_cppr(new_lines,digits = digits,
                                   line_weight = "weight",
                                   attrs = TRUE, direction = NULL)


  graph <- graph_result$graph
  nodes <- graph_result$spvertices


  node_id <- closest_points(snapped_events, nodes)
  snapped_events$vertex_id <- nodes$ref[node_id]


  ## step5 : calculating the distance matrix
  snappedA <- subset(snapped_events, snapped_events$type == 'A')
  snappedB <- subset(snapped_events, snapped_events$type == 'B')
  dist_mat <- cppRouting::get_distance_matrix(graph, from = snappedB$vertex_id, to = snappedA$vertex_id)


  ## step6 : calculating the k and g function
  if(calc_g_func){
    k_g_vals <- kgfunc_cpp2(dist_mat,start,end,step,width,Lt,na,wc = snappedA$weight, wr = snappedB$weight, cross = TRUE)
    k_vals <- k_g_vals[,1]
    g_vals <- k_g_vals[,2]
  }else{
    k_vals <- kfunc_cpp2(dist_mat,start,end,step,Lt,na,wc = snappedA$weight, wr = snappedB$weight, cross = TRUE)
  }


  ## step7 : generate the permutations
  if (verbose){
    print("Calculating the simulations ...")
  }

  if(verbose){
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }

  # the case where we can must edit the network and create small chunks
  if (is.null(resolution)==FALSE){

    # we start by creating the lixels
    new_lines2 <- lixelize_lines(new_lines, resolution, mindist = resolution / 2.0)
    new_lines2$weight <- as.numeric(st_length(new_lines2))


    # we can now rebuild the network with these new lines
    graph_result <- build_graph_cppr(new_lines2,digits = digits,
                                     line_weight = "weight",
                                     attrs = TRUE, direction = NULL)

    graph <- graph_result$graph
    nodes <- graph_result$spvertices


    node_id <- closest_points(snapped_events, nodes)
    snapped_events$vertex_id <- nodes$id[node_id]
    snappedA <- subset(snapped_events, snapped_events$type == 'A')
    snappedB <- subset(snapped_events, snapped_events$type == 'B')

  }

  # we can now do the simulation by permutating the location of the events on
  # the network.

  all_values <- lapply(1:nsim, function(i){

    vidxA <- sample(graph$dict$ref, size = sum(snappedA$weight))
    vidxB <- sample(graph$dict$ref, size = sum(snappedB$weight))
    dist_mat <- cppRouting::get_distance_matrix(graph, from = vidxB, to = vidxA)

    if(calc_g_func){
      k_g_vals <- kgfunc_cpp2(dist_mat,start,end,step,width,Lt,na,wc = snappedA$weight, wr = snappedB$weight, cross = TRUE)
      return(k_g_vals)
    }else{
      k_vals <- kfunc_cpp2(dist_mat,start,end,step,Lt,na,wc = snappedA$weight, wr = snappedB$weight, cross = TRUE)
      return(k_vals)
    }

  })

  ## Then we can prepare the results
  obj <- prep_kfuncs_results(k_vals, g_vals, all_values, conf_int, calc_g_func,
                             cross = TRUE, dist_seq = seq(start, end, step),
                             return_sims = return_sims
                             )
  return(obj)
}




#' @title Network cross k and g functions (maturing, multicore)
#'
#' @description Calculate the cross k and g functions for a set of points on a
#'   network. For more details, see the document of the function cross_kfunctions.
#'
#'
#' @template kross_kfunctions-arg
#' @template common_kfunctions-arg
#' @param return_sims a boolean indicating if the simulated k and g values must also
#' be returned as matrices
#' @template grid_shape-arg
#'
#' @return A list with the following values : \cr \item{plotk}{ A
#'   ggplot2 object representing the values of the cross k-function}
#'   \item{plotg}{ A ggplot2 object representing the values of the cross
#'   g-function} \item{values}{ A DataFrame with the values used to build the
#'   plots}
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot geom_ribbon geom_path aes_string labs
#' @importFrom grDevices rgb
#' @export
#' @examples
#' \donttest{
#' data(main_network_mtl)
#' data(mtl_libraries)
#' data(mtl_theatres)
#' result <- new_cross_kfunctions.mc(main_network_mtl, mtl_theatres, mtl_libraries,
#'                            start = 0, end = 2500, step = 10, width = 250,
#'                            nsim = 50, conf_int = 0.05, digits = 2,
#'                            tol = 0.1, agg = NULL, verbose = FALSE)
#' }
new_cross_kfunctions.mc <- function(lines,
                                    pointsA,
                                    pointsB,
                                    start,
                                    end,
                                    step,
                                    width,
                                    nsim,
                                    conf_int = 0.05,
                                    digits = 2,
                                    tol = 0.1,
                                    resolution = NULL,
                                    agg = NULL,
                                    verbose = TRUE,
                                    return_sims = FALSE,
                                    calc_g_func = TRUE,
                                    grid_shape = c(1,1)){

  ## step0 : clean the points
  if(verbose){
    print("Preparing data ...")
  }

  probs <- NULL

  st_geometry(pointsA) <- 'geometry'
  st_geometry(pointsB) <- 'geometry'

  pointsA$weight <- rep(1,nrow(pointsA))
  pointsA <- clean_events(pointsA,digits,agg)
  pointsA$goid <- seq_len(nrow(pointsA))

  pointsB$weight <- rep(1,nrow(pointsB))
  pointsB <- clean_events(pointsB,digits,agg)
  pointsB$goid <- seq_len(nrow(pointsB))

  na <- sum(pointsA$weight)
  nb <- sum(pointsB$weight)

  ## step1 : clean the lines
  lines$length <- as.numeric(st_length(lines))
  lines <- subset(lines, lines$length>0)
  lines$oid <- seq_len(nrow(lines))


  ##______________________
  # Gridding process

  ## we will now define a grid
  grid <- build_grid(grid_shape = grid_shape, spatial = list(lines, pointsA, pointsB))

  ## we create two spatial indices to make the requests faster
  tree_pointsA <- build_quadtree(pointsA)
  tree_pointsB <- build_quadtree(pointsB)
  tree_lines <- build_quadtree(lines)


  ## and split my data according to this grid
  selections <- lapply(1:nrow(grid),function(i){
    square <- grid[i,]

    # selecting the points in the grid
    sel_points <- spatial_request(square,tree_pointsB,pointsB)
    # if there is no sampling points in the rectangle, then return NULL
    if(nrow(sel_points)==0){
      return(NULL)
    }
    # selecting the events in a buffer
    bw <- end + width
    buff <- st_buffer(square,dist = bw)

    sel_points2 <- spatial_request(buff,tree_pointsA,pointsA)

    # selecting the lines in a buffer
    sel_lines <- spatial_request(buff,tree_lines,lines)
    sel_lines$oid <- 1:nrow(sel_lines)

    return(list(
      'lines' = sel_lines,
      'origins' = sel_points,
      'destinations' = sel_points2
    ))

  })

  ##______________________
  # preparing values for randomisations
  Lt <- sum(as.numeric(st_length(lines)))

  selections <- selections[lengths(selections) != 0]


  ##______________________
  # applying the function on the grid

  # in this step, we will perform the network k and g functions on each element
  # of the grid. The idea is to obtain the counting matrices for each square and then
  # to calculate the means at each distance. The randomization is also done locally
  # by sampling the vertices. We must take into account the local length of the network
  # to do so.

  results <- lapply(1:length(selections), function(i){

    # we extract the selections first
    sel <- selections[[i]]
    sub_lines <- sel$lines
    origins <- sel$origins
    destinations <- sel$destinations

    # then, we snapp the points on the lines
    origins$type <- "B"
    destinations$type <- "A"

    all_events <- rbind(origins[c("type","goid","weight")],
                        destinations[c("type","goid","weight")])

    snapped_events <- snapPointsToLines2(all_events, lines, idField = "oid")

    new_lines <- split_lines_at_vertex(sub_lines, snapped_events,
                                       snapped_events$nearest_line_id, tol)

    # we create a graph
    new_lines$length <- as.numeric(st_length(new_lines))
    new_lines <- subset(new_lines,new_lines$length>0)



    new_lines$oid <- seq_len(nrow(new_lines))
    new_lines <- new_lines[c("length","oid")]
    local_Lt <- sum(as.numeric(st_length(new_lines)))

    new_lines$weight <- as.numeric(st_length(new_lines))

    graph_result <- build_graph_cppr(new_lines,digits = digits,
                                     line_weight = "weight",
                                     attrs = TRUE, direction = NULL)

    graph <- graph_result$graph
    nodes <- graph_result$spvertices

    # we must find the the locations of the events on the graph
    node_id <- closest_points(snapped_events, nodes)
    snapped_events$vertex_id <- nodes$ref[node_id]

    start_nodes <- subset(snapped_events, snapped_events$type == 'B')
    end_nodes <- subset(snapped_events, snapped_events$type == 'A')

    # and calculate the distance matrix
    dist_mat <- cppRouting::get_distance_matrix(graph, from = start_nodes$vertex_id, to = end_nodes$vertex_id)
    values <- list()
    values$sumw <- sum(start_nodes$weight)

    # with this distance matrix, we can calculate the counting matrix for this square
    if(calc_g_func){
      matrices <- kgfunc_counting(dist_mat,
                                  wc = snapped_events$weight,
                                  wr = start_nodes$weight,
                                  breaks = rev(seq_num3(start, end, step)),
                                  width = width / 2, cross = TRUE)
      values$k_counts <- matrices[[1]]
      values$g_counts <- matrices[[2]]

    }else{
      values$k_counts <- kfunc_counting(dist_mat,
                                        wc = snapped_events$weight,
                                        wr = start_nodes$weight,
                                        breaks = rev(seq_num3(start, end, step)),
                                        cross = TRUE
      )

    }


    ## Excellent ! The next step is to do the permutation locally. To do so, we must first apply the
    ## sampling strategy provided by the user.
    if (is.null(resolution)==FALSE){

      # we start by creating the lixels
      new_lines2 <- lixelize_lines(new_lines, resolution, mindist = resolution / 2.0)
      new_lines2$weight <- as.numeric(st_length(new_lines2))


      # we can now rebuild the network with these new lines
      graph_result <- build_graph_cppr(new_lines2,digits = digits,
                                       line_weight = "weight",
                                       attrs = TRUE, direction = NULL)

      graph <- graph_result$graph
      nodes <- graph_result$spvertices


      node_id <- closest_points(snapped_events, nodes)
      snapped_events$vertex_id <- nodes$id[node_id]

    }


    #sample_size <- round((local_Lt / Lt) * N)
    frac <- local_Lt / Lt

    my_breaks <- rev(seq_num3(start, end, step))

    simulations <- lapply(1:nsim, function(j){

      vidxA <- sample(graph$dict$ref, size = round(na * frac) )
      vidxB <- sample(graph$dict$ref, size = round(nb * frac))

      dist_mat <- cppRouting::get_distance_matrix(graph, from = vidxB, to = vidxA)

      values <- list()
      #values$sumw <- sample_size

      # with this distance matrix, we can calculate the counting matrix for this square
      if(calc_g_func){
        matrices <- kgfunc_counting(dist_mat,
                                    wc = rep(1, ncol(dist_mat)),
                                    wr = rep(1, nrow(dist_mat)),
                                    breaks = my_breaks,
                                    width = width / 2, cross = TRUE)
        values$k_counts <- matrices[[1]]
        values$g_counts <- matrices[[2]]

      }else{
        values$k_counts <- kfunc_counting(dist_mat,
                                          wc = rep(1, ncol(dist_mat)),
                                          wr = rep(1, nrow(dist_mat)),
                                          breaks = my_breaks, cross = TRUE)
      }
      return(values)


    })

    values$simulations <- simulations

    return(values)

  })

  ## its is time to combine the values obtained by the different tiles

  # let me start with the counts obtained for the k function
  k_counts <- do.call(rbind, lapply(results, function(x){x$k_counts}))

  # the countings must be transformed in kvalues. For each distance, it is the mean of
  # the reached points multiplied by t1
  t1 <- 1.0/((na)/Lt)
  k_vals <- colSums(k_counts) / nb * t1
  k_vals <- rev(k_vals)

  if(calc_g_func){
    g_counts <- do.call(rbind, lapply(results, function(x){x$g_counts}))
    g_vals <- colSums(g_counts) / nb * t1
    g_vals <- rev(g_vals)
  }



  # and then the simulations
  # all values must be a list with an element per simulation. Each element will
  # be a matrix of two columns. One with the k values, and the second with the g values

  all_values <- lapply(1:nsim, function(i){
    k_sim_vals <- do.call(rbind, lapply(results, function(x){
      x$simulations[[i]]$k_counts
    }))
    sims_k <- colSums(k_sim_vals) / nb * t1
    sims_k <- rev(sims_k)
    if(calc_g_func){
      g_sim_vals <- do.call(rbind, lapply(results, function(x){
        x$simulations[[i]]$g_counts
      }))
      sims_g <- colSums(g_sim_vals) / nb * t1
      sims_g <- rev(sims_g)
      return(cbind(sims_k,sims_g))
    }else{
      return(sims_k)
    }
  })


  ## Then we can prepare the results
  obj <- prep_kfuncs_results(k_vals, g_vals, all_values, conf_int, calc_g_func,
                             cross = FALSE, dist_seq = seq(start, end, step),
                             return_sims = return_sims)


  return(obj)
}

