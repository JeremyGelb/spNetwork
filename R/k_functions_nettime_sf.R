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
#'
#' @return A list with the following values : \cr
#'  \itemize{
#'       \item{obs_k}{ A matrix with the observed k-values}
#'       \item{lower_k}{ A matrix with the lower bounds of the simulated k-values}
#'       \item{upper_k}{ A matrix with the upper bounds of the simulated k-values}
#'       \item{obs_g}{ A matrix with the observed g-values}
#'       \item{lower_g}{ A matrix with the lower bounds of the simulated g-values}
#'       \item{upper_g}{ A matrix with the upper bounds of the simulated g-values}
#'       \item{distances_net}{ A vector with the used network distances}
#'       \item{distances_time}{ A vector with the used time distances}
#'       }
#' @importFrom stats quantile
#' @export
#' @examples
#' \donttest{
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' bike_accidents <- sf::st_read(eventsgpkg,layer="bike_accidents")
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- sf::st_read(networkgpkg,layer="mtl_network")
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
k_nt_functions <- function(lines, points, points_time,
                           start_net, end_net, step_net, width_net,
                           start_time, end_time, step_time, width_time,
                           nsim, conf_int = 0.05, digits = 2, tol = 0.1,
                           resolution = NULL, agg = NULL, verbose = TRUE){

  ## step0 : clean the points
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
  graph_result <- build_graph(new_lines,digits = digits,
                              line_weight = "weight",
                              attrs = TRUE)

  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
                                                           name = "probs")
  snapped_events$vertex_id <- closest_points(snapped_events, nodes)

  ## step5 : calculating the distance matrix
  dist_mat <- igraph::distances(graph,v = snapped_events$vertex_id,
                                to = snapped_events$vertex_id)

  ## step 5.5 I must now deal with the dupplicated ids
  ## minus 1 because c++ indexing starts at 0
  dist_mat_net <- extend_matrix_by_ids(dist_mat, points$goid, points$locid-1)

  ## and generate a matrix with the time distances !
  dist_mat_time <- as.matrix(stats::dist(points_time))

  ## step6 : calcualte the kfunction and the g function
  if (verbose){
    print("Calculating k and g functions ...")
  }
  k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time,
                          start_net, end_net, step_net,
                          start_time, end_time, step_time,
                          Lt, Tt, n, w = points$weight)

  g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time,
                          start_net, end_net, step_net, width_net,
                          start_time, end_time, step_time, width_time,
                          Lt, Tt, n,w = points$weight)

  ## step7 : generate the permutations
  if (verbose){
    print("Calculating the simulations ...")
  }
  w <- rep(1,times = n)
  if(verbose){
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }

  dims_time <- dim(dist_mat_time)

  # the case where we can simplified the situation
  if (is.null(resolution)==FALSE){
    dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
                                           edge_df = graph_result$spedges,
                                           n = n,
                                           resolution = resolution,
                                           nsim = nsim)

    all_values <- lapply(1:nsim,function(i){
      dist_mat_net <- dist_matrices[[i]]

      dist_mat_time <- runif(n = prod(dims_time),
                             min = min(points_time), max = max(points_time))
      dim(dist_mat_time) <- dims_time

      k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time,
                              start_net, end_net, step_net,
                              start_time, end_time, step_time,
                              Lt, Tt, n, w = w)

      g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time,
                              start_net, end_net, step_net, width_net,
                              start_time, end_time, step_time, width_time,
                              Lt, Tt, n, w = w)
      if(verbose){
        setTxtProgressBar(pb, i)
      }
      return(list(k_vals,g_vals))
    })

  }else{
    # the case where we can not simplified the situation
    all_values <- lapply(1:nsim,function(i){
      dist_mat_net <- randomize_distmatrix(graph_result$graph,graph_result$spedges,n)
      dist_mat_time <- runif(n = prod(dims_time),
                             min = min(points_time), max = max(points_time))
      dim(dist_mat_time) <- dims_time

      k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time,
                              start_net, end_net, step_net,
                              start_time, end_time, step_time,
                              Lt, Tt, n, w = w)

      g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time,
                              start_net, end_net, step_net, width_net,
                              start_time, end_time, step_time, width_time,
                              Lt, Tt, n, w = w)
      if(verbose){
        setTxtProgressBar(pb, i)
      }
      return(list(k_vals,g_vals))
    })
  }

  ## step8 : extract the k_vals and g_vals matrices
  L1 <- lapply(all_values,function(i){return(i[[1]])})
  L2 <- lapply(all_values,function(i){return(i[[2]])})
  L1[["along"]] <- 3
  L2[["along"]] <- 3
  k_mat <- do.call(abind::abind, L1)
  g_mat <- do.call(abind::abind, L2)

  ## step9 : calculating the summary stats
  upper <- 1 - conf_int / 2
  lower <- conf_int / 2
  k_stats_lower <- apply(k_mat,MARGIN = c(1,2), function(i){
    return(quantile(i,probs = lower))
  })
  k_stats_upper <- apply(k_mat,MARGIN = c(1,2), function(i){
    return(quantile(i,probs = upper))
  })
  g_stats_upper <- apply(g_mat,MARGIN = c(1,2), function(i){
    return(quantile(i,probs = upper))
  })
  g_stats_lower <- apply(g_mat,MARGIN = c(1,2), function(i){
    return(quantile(i,probs = lower))
  })

  seq_net <- seq_num2(start_net, end_net, step_net)
  seq_time <- seq_num2(start_time, end_time, step_time)

  row.names(k_vals) <- seq_net
  row.names(k_stats_lower) <- seq_net
  row.names(k_stats_upper) <- seq_net
  row.names(g_vals) <- seq_net
  row.names(g_stats_lower) <- seq_net
  row.names(g_stats_upper) <- seq_net

  colnames(k_vals) <- seq_time
  colnames(k_stats_lower) <- seq_time
  colnames(k_stats_upper) <- seq_time
  colnames(g_vals) <- seq_time
  colnames(g_stats_lower) <- seq_time
  colnames(g_stats_upper) <- seq_time

  results <- list(
    "obs_k" = k_vals,
    "lower_k" = k_stats_lower,
    "upper_k" = k_stats_upper,
    "obs_g" = g_vals,
    "lower_g" = g_stats_lower,
    "upper_g" = g_stats_upper,
    "distances_net" = seq_num2(start_net, end_net, step_net),
    "distances_time" = seq_num2(start_time, end_time, step_time)
  )
  # see this plot https://www.researchgate.net/publication/320760197_Spatiotemporal_Point_Pattern_Analysis_Using_Ripley%27s_K_Function
  return(results)
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
#'
#' @return A list with the following values : \cr
#'  \itemize{
#'       \item{obs_k}{ A matrix with the observed k-values}
#'       \item{lower_k}{ A matrix with the lower bounds of the simulated k-values}
#'       \item{upper_k}{ A matrix with the upper bounds of the simulated k-values}
#'       \item{obs_g}{ A matrix with the observed g-values}
#'       \item{lower_g}{ A matrix with the lower bounds of the simulated g-values}
#'       \item{upper_g}{ A matrix with the upper bounds of the simulated g-values}
#'       \item{distances_net}{ A vector with the used network distances}
#'       \item{distances_time}{ A vector with the used time distances}
#'       }
#' @importFrom stats quantile
#' @export
k_nt_functions.mc <- function(lines, points, points_time,
                           start_net, end_net, step_net, width_net,
                           start_time, end_time, step_time, width_time,
                           nsim, conf_int = 0.05, digits = 2, tol = 0.1,
                           resolution = NULL, agg = NULL, verbose = TRUE){

  ## step0 : clean the points
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
  Lt <- sum(as.numeric(st_length(new_lines)))
  Tt <- max(points_time) - min(points_time)

  new_lines$weight <- as.numeric(st_length(new_lines))

  ## step4 : building the graph for the real case
  graph_result <- build_graph(new_lines,digits = digits,
                              line_weight = "weight",
                              attrs = TRUE)
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
                                                           name = "probs")
  snapped_events$vertex_id <- closest_points(snapped_events, nodes)

  ## step5 : calculating the distance matrix
  dist_mat <- igraph::distances(graph,v = snapped_events$vertex_id,
                                to = snapped_events$vertex_id)

  ## step 5.5 I must now deal with the dupplicated ids
  ## minus 1 because c++ indexing starts at 0
  dist_mat_net <- extend_matrix_by_ids(dist_mat, points$goid, points$locid-1)
  print(round(dist_mat_net))

  ## and generate a matrix with the time distances !
  dist_mat_time <- as.matrix(stats::dist(points_time))

  ## step6 : calcualte the kfunction and the g function
  if (verbose){
    print("Calculating k and g functions ...")
  }

  k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time,
                          start_net, end_net, step_net,
                          start_time, end_time, step_time,
                          Lt, Tt, n, w = points$weight)

  g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time,
                          start_net, end_net, step_net, width_net,
                          start_time, end_time, step_time, width_time,
                          Lt, Tt, n,w = points$weight)

  ## step7 : generate the permutations
  if (verbose){
    print("Calculating the simulations ...")
  }
  w <- rep(1,times = n)
  if(verbose){
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }

  dims_time <- dim(dist_mat_time)
  sim_seq <- 1:nsim
  graph <- graph_result$graph
  edgesdf <- st_drop_geometry(graph_result$spedges)

  # the classical way
  if (is.null(resolution)){
    if(verbose){
      progressr::with_progress({
        p <- progressr::progressor(along = sim_seq)
        all_values <- future.apply::future_lapply(sim_seq, function(i){

          dist_mat_net <- randomize_distmatrix(graph_result$graph,graph_result$spedges,n)
          dist_mat_time2 <- runif(n = prod(dims_time),
                                 min = min(points_time), max = max(points_time))
          dim(dist_mat_time2) <- dims_time

          k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time2,
                                  start_net, end_net, step_net,
                                  start_time, end_time, step_time,
                                  Lt, Tt, n, w = w)

          g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time2,
                                  start_net, end_net, step_net, width_net,
                                  start_time, end_time, step_time, width_time,
                                  Lt, Tt, n, w = w)
          return(list(k_vals,g_vals))
        },future.packages = c("igraph"))
      })
    }else{
      all_values <- future.apply::future_lapply(sim_seq, function(i){
        dist_mat_net <- randomize_distmatrix(graph_result$graph,graph_result$spedges,n)
        dist_mat_time2 <- runif(n = prod(dims_time),
                                min = min(points_time), max = max(points_time))
        dim(dist_mat_time2) <- dims_time

        k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time2,
                                start_net, end_net, step_net,
                                start_time, end_time, step_time,
                                Lt, Tt, n, w = w)

        g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time2,
                                start_net, end_net, step_net, width_net,
                                start_time, end_time, step_time, width_time,
                                Lt, Tt, n, w = w)
        p(sprintf("i=%g", i))
        return(list(k_vals,g_vals))
      },future.packages = c("igraph"))

    }
  }else{
    # the simplified way
    ## first : generating the matrices
    if(verbose){
      print("generating the randomized distance matrices...")
    }
    dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
                                           edge_df = graph_result$spedges,
                                           n = n,
                                           resolution = resolution,
                                           nsim = nsim)
    if(verbose){
      print("calculating the k and g functions for the randomized matrices..")
    }
    all_values <- future.apply::future_lapply(dist_matrices, function(dist_mat_net){
      dist_mat_time2 <- runif(n = prod(dims_time),
                             min = min(points_time), max = max(points_time))
      dim(dist_mat_time2) <- dims_time

      k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time2,
                              start_net, end_net, step_net,
                              start_time, end_time, step_time,
                              Lt, Tt, n, w = w)

      g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time2,
                              start_net, end_net, step_net, width_net,
                              start_time, end_time, step_time, width_time,
                              Lt, Tt, n, w = w)
      return(list(k_vals,g_vals))
    },future.packages = c("igraph"))

  }

  ## step8 : extract the k_vals and g_vals matrices
  L1 <- lapply(all_values,function(i){return(i[[1]])})
  L2 <- lapply(all_values,function(i){return(i[[2]])})
  L1[["along"]] <- 3
  L2[["along"]] <- 3
  k_mat <- do.call(abind::abind, L1)
  g_mat <- do.call(abind::abind, L2)

  ## step9 : calculating the summary stats
  upper <- 1 - conf_int / 2
  lower <- conf_int / 2
  k_stats_lower <- apply(k_mat,MARGIN = c(1,2), function(i){
    return(quantile(i,probs = lower))
  })
  k_stats_upper <- apply(k_mat,MARGIN = c(1,2), function(i){
    return(quantile(i,probs = upper))
  })
  g_stats_upper <- apply(g_mat,MARGIN = c(1,2), function(i){
    return(quantile(i,probs = upper))
  })
  g_stats_lower <- apply(g_mat,MARGIN = c(1,2), function(i){
    return(quantile(i,probs = lower))
  })

  seq_net <- seq_num2(start_net, end_net, step_net)
  seq_time <- seq_num2(start_time, end_time, step_time)

  row.names(k_vals) <- seq_net
  row.names(k_stats_lower) <- seq_net
  row.names(k_stats_upper) <- seq_net
  row.names(g_vals) <- seq_net
  row.names(g_stats_lower) <- seq_net
  row.names(g_stats_upper) <- seq_net

  colnames(k_vals) <- seq_time
  colnames(k_stats_lower) <- seq_time
  colnames(k_stats_upper) <- seq_time
  colnames(g_vals) <- seq_time
  colnames(g_stats_lower) <- seq_time
  colnames(g_stats_upper) <- seq_time

  results <- list(
    "obs_k" = k_vals,
    "lower_k" = k_stats_lower,
    "upper_k" = k_stats_upper,
    "obs_g" = g_vals,
    "lower_g" = g_stats_lower,
    "upper_g" = g_stats_upper,
    "distances_net" = seq_num2(start_net, end_net, step_net),
    "distances_time" = seq_num2(start_time, end_time, step_time)
  )
  # see this plot https://www.researchgate.net/publication/320760197_Spatiotemporal_Point_Pattern_Analysis_Using_Ripley%27s_K_Function
  return(results)
}

