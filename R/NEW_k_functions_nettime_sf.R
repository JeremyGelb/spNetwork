# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### execution k functions in space-time ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Network k and g functions for spatio-temporal data (experimental, NOT READY FOR USE)
#'
#' @description Calculate the k and g functions for a set of points on a
#'   network and in time (experimental, NOT READY FOR USE).
#'
#' @details The k-function is a method to characterize the dispersion of a set
#'   of points. For each point, the numbers of other points in subsequent radii
#'   are calculated in both space and time. This empirical k-function can be more or less clustered
#'   than a k-function obtained if the points were randomly located . In
#'   a network, the network distance is used instead of the Euclidean distance.
#'   This function uses Monte Carlo simulations to assess if the points are
#'   clustered or dispersed. The function also calculates the
#'   g-function, a modified version of the k-function using rings instead of
#'   disks. The width of the ring must be chosen. The main interest is to avoid
#'   the cumulative effect of the classical k-function. This function is maturing,
#'   it works as expected (unit tests) but will probably be modified in the
#'   future releases (gain speed, advanced features, etc.).
#'
#' @template kfunctions-arg
#' @template common_kfunctions_nt-arg
#' @param points_time A numeric vector indicating when the point occured
#' @param calc_g_func A boolean indicating if the G function must also be calculated
#'
#' @return A list with the following values :
#'
#' * obs_k: A matrix with the observed k-values
#' * lower_k: A matrix with the lower bounds of the simulated k-values
#' * upper_k: A matrix with the upper bounds of the simulated k-values
#' * obs_g: A matrix with the observed g-values
#' * lower_g: A matrix with the lower bounds of the simulated g-values
#' * upper_g: A matrix with the upper bounds of the simulated g-values
#' * distances_net: A vector with the used network distances
#' * distances_time: A vector with the used time distances
#'
#' @importFrom stats quantile
#' @export
#' @md
#' @examples
#' \donttest{
#' data(mtl_network)
#' data(bike_accidents)
#'
#' # converting the Date field to a numeric field (counting days)
#' bike_accidents$Time <- as.POSIXct(bike_accidents$Date, format = "%Y/%m/%d")
#' start <- as.POSIXct("2016/01/01", format = "%Y/%m/%d")
#' bike_accidents$Time <- difftime(bike_accidents$Time, start, units = "days")
#' bike_accidents$Time <- as.numeric(bike_accidents$Time)
#'
#' values <- k_nt_functions(
#'       lines =  mtl_network,
#'       points = bike_accidents,
#'       points_time = bike_accidents$Time,
#'       start_net = 0 ,
#'       end_net = 2000,
#'       step_net = 10,
#'       width_net = 200,
#'       start_time = 0,
#'       end_time = 360,
#'       step_time = 7,
#'       width_time = 14,
#'       nsim = 50,
#'       conf_int = 0.05,
#'       digits = 2,
#'       tol = 0.1,
#'       resolution = NULL,
#'       agg = 15,
#'       verbose = TRUE)
#'}
k_nt_functions <- function(lines,
                           points,
                           points_time,
                           start_net,
                           end_net,
                           step_net,
                           width_net,
                           start_time,
                           end_time,
                           step_time,
                           width_time,
                           nsim,
                           conf_int = 0.05,
                           digits = 2,
                           tol = 0.1,
                           resolution = NULL,
                           agg = NULL,
                           verbose = TRUE,
                           calc_g_func = TRUE
                           ){

  ## step0 : clean the points
  probs <- NULL
  if (verbose){
    print("Preparing data ...")
  }
  n <- nrow(points)
  points$goid <- seq_len(nrow(points))
  points$weight <- rep(1,nrow(points))

  # we aggregate the points to have a fewer number of locations
  # but because of the time dimension, we have to keep them separated in the end
  agg_points <- clean_events(points,digits,agg)
  agg_points$locid <- 1:nrow(agg_points)

  # and we need now to find the new location of the original points
  xy_origin <- st_coordinates(points)
  xy_agg <- st_coordinates(agg_points)
  ids <- dbscan::kNN(xy_agg, k = 1, query = xy_origin)
  points$locid <- ids$id[,1]

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
  snapped_events <- snapPointsToLines2(agg_points, lines, idField = "oid")
  new_lines <- split_lines_at_vertex(lines, snapped_events,
                                     snapped_events$nearest_line_id, tol)

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

  # calculating the extent of the study area
  Lt <- sum(as.numeric(st_length(new_lines)))
  Tt <- max(points_time) - min(points_time)

  new_lines$weight <- as.numeric(st_length(new_lines))

  ## step4 : building the graph for the real case
  graph_result <- build_graph_cppr(new_lines,digits = digits,
                                   line_weight = "weight",
                                   attrs = TRUE, direction = NULL)

  graph <- graph_result$graph
  nodes <- graph_result$spvertices


  snapped_events$vertex_id <- closest_points(snapped_events, nodes)
  points$vertex_id <- snapped_events$vertex_id[match(points$locid, snapped_events$locid)]

  ## step5 : calculating the distance matrix
  dist_mat_net <- cppRouting::get_distance_matrix(graph, from = points$vertex_id, to = points$vertex_id)

  ## and generate a matrix with the time distances !
  dist_mat_time <- as.matrix(stats::dist(points_time))

  ## step6 : calcualte the kfunction and the g function
  if (verbose){
    print("Calculating k and g functions ...")
  }

  if(calc_g_func){

    k_g_vals <- k_g_nt_func_cpp2(dist_mat_net = dist_mat_net,
                                 dist_mat_time = dist_mat_time,
                                 start_net = start_net,
                                 end_net = end_net,
                                 step_net = step_net,
                                 start_time = start_time,
                                 end_time = end_time,
                                 step_time = step_time,
                                 width_net = width_net,
                                 width_time = width_time,
                                 cross = FALSE,
                                 Lt = Lt,
                                 Tt = Tt,
                                 n = n,
                                 wr = rep(1, nrow(dist_mat_net)),
                                 wc = rep(1, ncol(dist_mat_net))
    )
    kvals <- k_g_vals[[1]]
    gvals <- k_g_vals[[2]]

  }else{

    kvals <- k_nt_func_cpp2(dist_mat_net = dist_mat_net,
                                 dist_mat_time = dist_mat_time,
                                 start_net = start_net,
                                 end_net = end_net,
                                 step_net = step_net,
                                 start_time = start_time,
                                 end_time = end_time,
                                 step_time = step_time,
                                 cross = FALSE,
                                 Lt = Lt,
                                 Tt = Tt,
                                 n = n,
                                 wr = rep(1, nrow(dist_mat_net)),
                                 wc = rep(1, ncol(dist_mat_net))
    )
  }




  ## step7 : generate the permutations
  if (verbose){
    print("Calculating the simulations ...")
  }


  if(verbose){
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }


  # if the user provided a resolution on the network, then we must create chunks
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


    snapped_events$vertex_id <- closest_points(snapped_events, nodes)
    points$vertex_id <- snapped_events$vertex_id[match(points$locid, snapped_events$locid)]

  }

  # based on the network, we will randomize the location of the events and there datetime

  all_values <- lapply(1:nsim, function(i){

    if(verbose){
      setTxtProgressBar(pb,i)
    }

    vidx <- sample(graph$dict$ref, size = sum(snapped_events$weight))
    dist_mat <- cppRouting::get_distance_matrix(graph, from = vidx, to = vidx)
    time_mat <- as.matrix(stats::dist(sample(points_time, size = length(points_time), replace = T)))

    if(calc_g_func){

      k_g_vals <- k_g_nt_func_cpp2(dist_mat_net = dist_mat,
                                   dist_mat_time = time_mat,
                                   start_net = start_net,
                                   end_net = end_net,
                                   step_net = step_net,
                                   start_time = start_time,
                                   end_time = end_time,
                                   step_time = step_time,
                                   width_net = width_net,
                                   width_time = width_time,
                                   cross = FALSE,
                                   Lt = Lt,
                                   Tt = Tt,
                                   n = n,
                                   wr = rep(1, nrow(dist_mat)),
                                   wc = rep(1, ncol(dist_mat))
      )
      return(k_g_vals)

    }else{

      kvals <- k_nt_func_cpp2(dist_mat_net = dist_mat,
                              dist_mat_time = time_mat,
                              start_net = start_net,
                              end_net = end_net,
                              step_net = step_net,
                              start_time = start_time,
                              end_time = end_time,
                              step_time = step_time,
                              cross = FALSE,
                              Lt = Lt,
                              Tt = Tt,
                              n = n,
                              wr = rep(1, nrow(dist_mat)),
                              wc = rep(1, ncol(dist_mat))
      )
      return(list(kvals))

    }


  })



  ## step8 : extract the k_vals and g_vals matrices
  # to do so, we structure the values in a matrix

  upper <- 1 - conf_int / 2
  lower <- conf_int / 2

  L1 <- sapply(all_values,function(i){return( c(i[[1]]) )})
  v1 <- apply(L1, MARGIN = 1, FUN = function(x){
    return(quantile(x, probs = c(lower, upper)))
  })
  k_stats_lower <- v1[1,]
  k_stats_upper <- v1[2,]
  dim(k_stats_lower) <- dim(all_values[[1]][[1]])
  dim(k_stats_upper) <- dim(all_values[[1]][[1]])

  if(calc_g_func){

    L2 <- sapply(all_values,function(i){return( c(i[[2]]) )})
    v2 <- apply(L2, MARGIN = 1, FUN = function(x){
      return(quantile(x, probs = c(lower, upper)))
    })
    g_stats_lower <- v2[1,]
    g_stats_upper <- v2[2,]
    dim(g_stats_lower) <- dim(all_values[[1]][[1]])
    dim(g_stats_upper) <- dim(all_values[[1]][[1]])
  }


  seq_net <- seq_num2(start_net, end_net, step_net)
  seq_time <- seq_num2(start_time, end_time, step_time)

  row.names(kvals) <- seq_net
  row.names(k_stats_lower) <- seq_net
  row.names(k_stats_upper) <- seq_net
  colnames(kvals) <- seq_time
  colnames(k_stats_lower) <- seq_time
  colnames(k_stats_upper) <- seq_time

  results <- list(
    "obs_k" = kvals,
    "lower_k" = k_stats_lower,
    "upper_k" = k_stats_upper
  )

  if(calc_g_func){
    row.names(gvals) <- seq_net
    row.names(g_stats_lower) <- seq_net
    row.names(g_stats_upper) <- seq_net
    colnames(gvals) <- seq_time
    colnames(g_stats_lower) <- seq_time
    colnames(g_stats_upper) <- seq_time
    results$obs_g <- gvals
    results$lower_g <- g_stats_lower
    results$upper_g <- g_stats_upper
  }


  results$distances_net <- seq_num2(start_net, end_net, step_net)
  results$distances_time <- seq_num2(start_time, end_time, step_time)


  # see this plot https://www.researchgate.net/publication/320760197_Spatiotemporal_Point_Pattern_Analysis_Using_Ripley%27s_K_Function
  return(results)
}



#' @title Rervese the elements in a matrix
#'
#' @description reverse the order of the elements in a matrix both
#' column and row wise
#'
#' @param mat The matrix to reverse
#' @return A matrix
#' @keywords internal
#' @export
rev_matrix <- function(mat){
  mat2 <- rev(mat)
  dim(mat2) <- dim(mat)
  return(mat2)
}



#' @title Network k and g functions for spatio-temporal data (multicore, experimental, NOT READY FOR USE)
#'
#' @description Calculate the k and g functions for a set of points on a
#'   network and in time (multicore, experimental, NOT READY FOR USE).
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
#' @template common_kfunctions_nt-arg
#' @param points_time A numeric vector indicating when the point occured
#' @param calc_g_func A boolean indicating if the G function must also be calculated
#' @template grid_shape-arg
#'
#' @return A list with the following values :
#'
#' * obs_k: A matrix with the observed k-values
#' * lower_k: A matrix with the lower bounds of the simulated k-values
#' * upper_k: A matrix with the upper bounds of the simulated k-values
#' * obs_g: A matrix with the observed g-values
#' * lower_g: A matrix with the lower bounds of the simulated g-values
#' * upper_g: A matrix with the upper bounds of the simulated g-values
#' * distances_net: A vector with the used network distances
#' * distances_time: A vector with the used time distances
#'
#'
#' @importFrom stats quantile runif
#' @md
#' @export
k_nt_functions.mc <- function(lines,
                              points,
                              points_time,
                              start_net,
                              end_net,
                              step_net,
                              width_net,
                              start_time,
                              end_time,
                              step_time,
                              width_time,
                              nsim,
                              conf_int = 0.05,
                              digits = 2,
                              tol = 0.1,
                              resolution = NULL,
                              agg = NULL,
                              verbose = TRUE,
                              calc_g_func = TRUE,
                              grid_shape = c(1,1)
){

  #___________________________
  ## step0 : clean the points
  probs <- NULL
  
  if (verbose){
    print("Preparing data ...")
  }
  n <- nrow(points)
  points$goid <- seq_len(nrow(points))
  points$weight <- rep(1,nrow(points))

  # we aggregate the points to have a fewer number of locations
  # but because of the time dimension, we have to keep them separated in the end
  agg_points <- clean_events(points,digits,agg)
  agg_points$locid <- 1:nrow(agg_points)

  # and we need now to find the new location of the original points
  xy_origin <- st_coordinates(points)
  xy_agg <- st_coordinates(agg_points)
  ids <- dbscan::kNN(xy_agg, k = 1, query = xy_origin)
  points$locid <- ids$id[,1]


  #___________________________
  ## step1 : clean the lines

  if(is.null(probs)){
    lines$probs <- 1
  }else{
    lines$probs <- probs
  }

  lines$length <- as.numeric(st_length(lines))
  lines <- subset(lines, lines$length>0)
  lines$oid <- seq_len(nrow(lines))

  ##______________________________________
  # Gridding process

  ## we will now define a grid
  grid <- build_grid(grid_shape = grid_shape, spatial = list(lines, points))
  grid$grid <- as.character(grid$oid)
  grid$oid <- NULL

  ## we precalculate the spatial intersections

  # for the origins, we must retain exactly the points in each quadra
  # we will be looking in the locid and find back the corresponding points

  # first, we find for each location point (agg_points) in which quadra it fells in
  inter_points <- st_join(agg_points, grid)
  df <- inter_points[c('locid', 'grid')]
  df <- merge(
    df, sf::st_drop_geometry(points)[c('time', 'weight', 'locid', 'goid')], by = 'locid'
  )

  # NB : the original events can belong only to one quadra
  # if the spatial intersection is perfect, we attribute the point randomly to a quadra
  df <- data.table(df)
  df <- sf::st_as_sf(data.frame(df[, .SD[1:1], c('goid')]), crs = st_crs(points))

  # and then we create a named list splitting the original points
  inter_points <- split(df, df$grid)

  # for the lines and the destinations, we must retain the points
  # that are within the grid plus a buffer
  bw <- end_net + width_net
  grid_buff <- st_buffer(grid,dist = bw)

  # for the lines, we realize a a simple left spatial joined, and we split the data again
  inter_lines <- st_join(lines, grid_buff)
  inter_lines <- split(inter_lines, inter_lines$grid)

  # We must also find the points that fell into the buffers arround the quadras
  # but some points might be dupplicated !
  points$grid <- NULL
  inter_points2 <- st_join(agg_points, grid_buff)
  df2 <- inter_points2[c('locid', 'grid')]
  df2 <- merge(
    df2, sf::st_drop_geometry(points)[c('time', 'weight', 'locid')], by = 'locid'
  )

  # in inter_points2, we have the locations of the destination points
  inter_points2 <- split(df2, df2$grid)

  #NB : we also must know the total length of the network in each quadra.
  # this value will be used latter to calculate the local sample size during
  # randomization

  quadra_lines <- st_intersection(lines, grid)
  quadra_lines$length <- as.numeric(st_length(quadra_lines))
  quadra_lines <- data.table(sf::st_drop_geometry(quadra_lines))
  quad_lt <- data.frame(quadra_lines[,sum(length),by="grid"])


  ## and split my data according to this grid
  selections <- lapply(grid$grid,function(gid){

    # selecting the points in the grid
    sel_points <- inter_points[[gid]]
    # if there is no sampling points in the rectangle, then return NULL
    if(is.null(sel_points)){
      return(NULL)
    }
    # selecting the events in a buffer
    sel_points2 <- inter_points2[[gid]]

    # selecting the lines in a buffer
    sel_lines <- inter_lines[[gid]]

    this_lt <- quad_lt[quad_lt$grid == gid,2]

    return(list(
      'lines' = sel_lines,
      'origins' = sel_points,
      'destinations' = sel_points2,
      'quad_lt' = this_lt
    ))

  })
  selections <- selections[lengths(selections) != 0]

  # calculating the extent of the study area
  Lt <- sum(as.numeric(st_length(lines)))
  Tt <- max(points_time) - min(points_time)
  N <- sum(points$weight)

  p <- (N-1) / (Lt * Tt)


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
    quad_lt <- sel$quad_lt

    # we can then add the points to the lines
    snapped_events <- snapPointsToLines2(destinations, lines, idField = "oid")

    new_lines <- split_lines_at_vertex(lines, snapped_events,
                                       snapped_events$nearest_line_id, tol)

    new_lines$length <- as.numeric(st_length(new_lines))
    new_lines <- subset(new_lines,new_lines$length>0)

    new_lines <- remove_loop_lines(new_lines,digits)

    new_lines$oid <- seq_len(nrow(new_lines))
    new_lines <- new_lines[c("length","oid","probs")]

    new_lines$weight <- as.numeric(st_length(new_lines))
    local_Lt <- sum(new_lines$weight)

    ## building the graph for the real case
    graph_result <- build_graph_cppr(new_lines,digits = digits,
                                     line_weight = "weight",
                                     attrs = TRUE, direction = NULL)

    graph <- graph_result$graph
    nodes <- graph_result$spvertices

    ## we must now find for each origin and destination its vertex id
    origins$vertex_id <- closest_points(origins, nodes)
    destinations$vertex_id <- closest_points(destinations, nodes)

    ## calculating the distance matrix
    dist_mat_net <- cppRouting::get_distance_matrix(graph,
                                                    from = origins$vertex_id,
                                                    to = destinations$vertex_id)

    ## and generate a matrix with the time distances !
    dist_mat_time <- pair_dists(origins$time, destinations$time)

    #____________________________________________________________
    ## with these two matrices, we can generate the countings for the k and g function
    ## NOTE : the self counting is discarded by the c++ function.

    quad_values <- list()
    quad_values$sumw <- sum(origins$weight)

    # with this distance matrix, we can calculate the counting matrix for this square
    if(calc_g_func){
      arrays <- kgfunc_time_counting(dist_mat_net = dist_mat_net,
                                     dist_mat_time = dist_mat_time,
                                     wc = destinations$weight,
                                     wr = origins$weight,
                                     breaks_net =  rev(seq_num3(start_net, end_net, step_net)),
                                     breaks_time = rev(seq_num3(start_time, end_time, step_time)),
                                     width_net = width_net / 2.0,
                                     width_time = width_time / 2.0,
                                     cross = FALSE)

      quad_values$k_counts <- arrays[[1]]
      quad_values$g_counts <- arrays[[2]]

    }else{
      quad_values$k_counts <- kfunc_time_counting(dist_mat_net = dist_mat_net,
                                              dist_mat_time = dist_mat_time,
                                              wc = destinations$weight,
                                              wr = origins$weight,
                                              breaks_net =  rev(seq_num3(start_net, end_net, step_net)),
                                              breaks_time = rev(seq_num3(start_time, end_time, step_time)),
                                              cross = FALSE)

    }

    #____________________________________________________________
    # now that we have the countings, we will use the multiple cpus to generate
    # randomisations.

    # if the user gaves us a resolution, we must split the lines and rebuild the network

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

    }


    # we have now a network from whic we can sample positions and recalculate
    # the distances matrices. Once again, we must only calculate the countings
    # because the real k and g values will require to have all quadra calculated

    dest_sample_size <- round((local_Lt / Lt) * N)
    ori_sample_size <- round((quad_lt / Lt) * N)

    simulations <- lapply(1:nsim, function(j){

      ori_idx <- sample(graph$dict$ref, size = ori_sample_size, replace = TRUE)
      dest_idx <- sample(graph$dict$ref, size = dest_sample_size, replace = TRUE)

      dist_mat_net <- cppRouting::get_distance_matrix(graph, from = ori_idx, to = dest_idx)

      ori_time <- runif(n = ori_sample_size, min = min(points_time), max = max(points_time))
      dest_time <- runif(n = dest_sample_size, min = min(points_time), max = max(points_time))

      dist_time_mat <- pair_dists(ori_time, dest_time)

      values <- list()
      values$sumw <- ori_sample_size

      # with this distance matrix, we can calculate the counting matrix for this square
      # NB : we set cross to TRUE here because we know that the origin and the destinations
      # are not co-located
      if(calc_g_func){
        arrays <- kgfunc_time_counting(dist_mat_net = dist_mat_net,
                                       dist_mat_time = dist_time_mat,
                                       wc = rep(1, dest_sample_size),
                                       wr = rep(1, ori_sample_size),
                                       breaks_net =  rev(seq_num3(start_net, end_net, step_net)),
                                       breaks_time = rev(seq_num3(start_time, end_time, step_time)),
                                       width_net = width_net / 2.0,
                                       width_time = width_time / 2.0,
                                       cross = TRUE)
        values$k_counts <- arrays[[1]]
        values$g_counts <- arrays[[2]]

      }else{
        values$k_counts <- kfunc_time_counting(dist_mat_net = dist_mat_net,
                                               dist_mat_time = dist_mat_time,
                                               wc = rep(1, dest_sample_size),
                                               wr = rep(1, ori_sample_size),
                                               breaks_net =  rev(seq_num3(start_net, end_net, step_net)),
                                               breaks_time = rev(seq_num3(start_time, end_time, step_time)),
                                               cross = TRUE)
      }
      return(values)


    })


    quad_values$simulations <-simulations

    return(quad_values)
  })


  # That was a massive piece of code. We must then process the results we got from each chunk
  # to calculate the k and g values and the results from the simulations.

  # ______________________________
  # obtaining the real k and g values

  # the counts are given as arrays in each chunk. For each tube of the array, we have the counts
  k_counts <- do.call(abind::abind, lapply(results, function(x){x$k_counts}))
  # this is calculating the mean for each tube and return it in the expected format
  kmeans <- t(colMeans(aperm(k_counts, c(3,2,1))))
  kvals <- (1/p) * kmeans

  if(calc_g_func){
    g_counts <- do.call(abind::abind, lapply(results, function(x){x$g_counts}))
    gmeans <- t(colMeans(aperm(g_counts, c(3,2,1))))
    gvals <- (1/p) * gmeans
  }

  upper <- 1 - conf_int / 2
  lower <- conf_int / 2

  # ______________________________
  # obtaining the randomized k and g values

  k_sims <- lapply(1:nsim, function(i){
    this_sim <- do.call(abind::abind, lapply(results, function(x){x$simulations[[i]]$k_counts}))
    this_sim <- t(colMeans(aperm(this_sim, c(3,2,1))))
    this_sim <- (1/p) * this_sim
    return(this_sim)
  })



  # to get the confidence interval for each cell in the matrix of the k values,
  # I will stack their values in a matrix
  k_sims_mat <- sapply(k_sims, FUN = c)
  V1 <- apply(k_sims_mat, MARGIN = 1, FUN = function(x){
    return(quantile(x, probs = c(lower, upper)))
  })
  k_stats_lower <- V1[1,]
  k_stats_upper <- V1[2,]
  dim(k_stats_lower) <- dim(k_sims[[1]])
  dim(k_stats_upper) <- dim(k_sims[[1]])

  if(calc_g_func){

    g_sims <- lapply(1:nsim, function(i){
      this_sim <- do.call(abind::abind, lapply(results, function(x){x$simulations[[i]]$g_counts}))
      this_sim <- t(colMeans(aperm(this_sim, c(3,2,1))))
      this_sim <- (1/p) * this_sim
      return(this_sim)
    })



    # to get the confidence interval for each cell in the matrix of the k values,
    # I will stack their values in a matrix
    g_sims_mat <- sapply(g_sims, FUN = c)
    V2 <- apply(g_sims_mat, MARGIN = 1, FUN = function(x){
      return(quantile(x, probs = c(lower, upper)))
    })
    g_stats_lower <- V2[1,]
    g_stats_upper <- V2[2,]
    dim(g_stats_lower) <- dim(g_sims[[1]])
    dim(g_stats_upper) <- dim(g_sims[[1]])

  }


  seq_net <- seq_num2(start_net, end_net, step_net)
  seq_time <- seq_num2(start_time, end_time, step_time)

  # I just need to reverse all the elements in the matrix
  kvals <- rev_matrix(kvals)
  k_stats_upper <- rev_matrix(k_stats_upper)
  k_stats_lower <- rev_matrix(k_stats_lower)

  row.names(kvals) <- seq_net
  row.names(k_stats_lower) <- seq_net
  row.names(k_stats_upper) <- seq_net
  colnames(kvals) <- seq_time
  colnames(k_stats_lower) <- seq_time
  colnames(k_stats_upper) <- seq_time

  results <- list(
    "obs_k" = kvals,
    "lower_k" = k_stats_lower,
    "upper_k" = k_stats_upper
  )

  if(calc_g_func){

    gvals <- rev_matrix(gvals)
    g_stats_upper <- rev_matrix(g_stats_upper)
    g_stats_lower <- rev_matrix(g_stats_lower)

    row.names(gvals) <- seq_net
    row.names(g_stats_lower) <- seq_net
    row.names(g_stats_upper) <- seq_net
    colnames(gvals) <- seq_time
    colnames(g_stats_lower) <- seq_time
    colnames(g_stats_upper) <- seq_time
    results$obs_g <- gvals
    results$lower_g <- g_stats_lower
    results$upper_g <- g_stats_upper

  }

  results$distances_net <- seq_num2(start_net, end_net, step_net)
  results$distances_time <- seq_num2(start_time, end_time, step_time)


  # see this plot https://www.researchgate.net/publication/320760197_Spatiotemporal_Point_Pattern_Analysis_Using_Ripley%27s_K_Function
  return(results)
}

