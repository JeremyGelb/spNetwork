# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### base k functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# OLDER FUNCTION REPLACED BY SOME C++
# kfunc <- function(dist_mat,start,end,step,Lt,n,w){
#   breaks <- seq(start,end,step)
#   #t1 <- Lt/(n*(n-1))
#   t1 <- (n-1)/Lt
#   k_values <- sapply(breaks,function(dist){
#     int_mat <- t(t(dist_mat<=dist) * w)
#     diag(int_mat) <- 0
#     tot <- (rowSums(int_mat))
#     k <- t1 * sum(tot)
#     return(k)
#   })
#   return(k_values)
# }

# OLDER FUNCTION REPLACED BY SOME C++
# gfunc <- function(dist_mat,start,end,step,width,Lt,n,w){
#   breaks <- seq(start,end,step)
#   width <- width/2
#   #t1 <- Lt/(n*(n-1))
#   t1 <- (n-1)/Lt
#   k_values <- sapply(breaks,function(dist){
#     int_mat <- t(t(dist_mat<=dist+width & dist_mat>=dist-width) * w)
#     diag(int_mat) <- 0
#     tot <- (rowSums(int_mat))
#     k <- t1 * sum(tot)
#     return(k)
#   })
#   return(k_values)
# }


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### base cross k functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# OLDER FUNCTION REPLACED BY SOME C++
# cross_kfunc <- function(dist_mat,start,end,step,Lt,na,nb,wa,wb){
#   breaks <- seq(start,end,step)
#   #t1 <- Lt/(na*nb)
#   t1 <- (na/Lt)
#
#   # note : in the matrix, the rows are the b points
#   # and the columns are the a points
#   k_values <- sapply(breaks,function(dist){
#     # applying the a weight (row wise)
#     int_mat <- sweep((dist_mat<=dist), MARGIN = 2, FUN = "*", wa)
#     # applying the b weight (col wise)
#     #int_mat <- t(t(int_mat) * wb)
#     int_mat <- sweep(int_mat,MARGIN = 1, FUN = "*", wb)
#     tot <- (rowSums(int_mat))
#     k <- t1 * sum(tot)
#     return(k)
#   })
#   return(k_values)
# }


# OLDER FUNCTION REPLACED BY SOME C++
# cross_gfun <- function(dist_mat,start,end,step,width,Lt,na,nb,wa,wb){
#   breaks <- base::seq(from = start,to = end, by = step)
#   width <- width/2
#   #t1 <- Lt/(na*nb)
#   t1 <- (na/Lt)
#   g_values <- sapply(breaks,function(dist){
#     d1 <- dist + width
#     d2 <- dist - width
#     int_mat <- (dist_mat <= d1 & dist_mat >= d2)
#     int_mat <- base::sweep(int_mat, MARGIN = 2, FUN = "*", wa)
#     int_mat <- base::sweep(int_mat,MARGIN = 1, FUN = "*", wb)
#     tot <- rowSums(int_mat)
#     k <- t1 * sum(tot)
#     return(k)
#   })
#   return(g_values)
# }



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### randomization functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Points on network randomization simplified
#'
#' @description Randomize location of points on a network.
#'
#' @param graph An graph object from igraph
#' @param edge_df A DataFrame describing the edges
#' @param n The number of point
#' @param resolution The maximum size of the network edges
#' @param nsim The number of distance matrices to generate
#' @param start_vert The vertices from which the distances will be calculated.
#' if null, then the distances are calculated from the generated locations.
#'
#' @return A numeric matrix with the distances between points
#' @keywords internal
#' @importFrom igraph E edge_attr
#' @importFrom stats runif
#' @examples
#' #This is an internal function, no example provided
randomize_distmatrix2 <- function(graph, edge_df, n, resolution, nsim, start_vert = NULL){
  ## step1 : generate all the candidate nodes on the graph
  ## a : select all the edges with a length superior to resolution
  edit_edge <- subset(edge_df, edge_df$weight > resolution)
  all_names <- names(igraph::V(graph))

  ## step2 : create the new needed vertices and edges
  vertices_and_distances <- lapply(1:nrow(edit_edge), function(i){
    this_edge <- edit_edge[i,]
    start_node <- all_names[this_edge$start_oid]
    end_node <- all_names[this_edge$end_oid]
    dists <- rep(resolution, floor(this_edge$weight/resolution))
    if(sum(dists) == this_edge$weight){
      dists <- dists[1:(length(dists)-1)]
    }
    names(dists) <- paste("fict",i,1:length(dists), sep="_")
    return(data.frame(starts = c(start_node, names(dists)),
                      ends = c(names(dists), end_node),
                      weight = c(dists, this_edge$weight - sum(dists))))
  })

  all_elements <- do.call(rbind, vertices_and_distances)
  all_elements <- subset(all_elements,all_elements$weight>0)

  ## step 3 : creating a new graph with the new vertices and edges
  new_graph <- igraph::graph_from_data_frame(all_elements, directed = FALSE)

  ## then merging it with the previous graph
  tot_graph <- igraph::union(graph,new_graph, byname = TRUE)

  ## step 4 correcting the weights
  ws <- igraph::E(tot_graph)
  df_tmp <- data.frame("w1" = ws$weight_1,
                       "w2" = ws$weight_2)

  df_tmp[is.na(df_tmp$w1),"w1"] <- 0
  df_tmp[is.na(df_tmp$w2),"w2"] <- 0
  tot_graph <- igraph::set_edge_attr(tot_graph, "weight",
                                     value = df_tmp$w1 + df_tmp$w2,
                                     index = igraph::E(tot_graph))


  ## and now, calculating the distance matrices
  verts <- as.numeric(igraph::V(tot_graph))
  dist_matrices <- lapply(1:nsim, function(i){
    new_vert <- sample(verts,size = n, replace = F)
    if (is.null(start_vert)){
      dist_mat <- igraph::distances(tot_graph,v = new_vert, to = new_vert)
    }else{
      dist_mat <- igraph::distances(tot_graph,v = start_vert, to = new_vert)

    }
  })
  return(dist_matrices)
}


#' @title Points on network randomization
#'
#' @description Randomize location of points on a network.
#'
#' @param graph An graph object from igraph
#' @param edge_df A DataFrame describing the edges
#' @param n The number of point
#' @param start_vert The vertices from which the distances will be calculated.
#' if null, then the distances are calculated from the generated locations.
#'
#' @return A numeric matrix with the distances between points
#' @keywords internal
#' @importFrom igraph E edge_attr
#' @importFrom stats runif
#' @examples
#' #This is an internal function, no example provided
randomize_distmatrix <- function(graph, edge_df, n, start_vert = NULL){

  vec_runif <- Vectorize(runif, vectorize.args = c("max"))

  ## prefered case where points are randomly located on edges

  #a. selecting the edges that will have points
  sel_edges_id <- sample(edge_df$edge_id,
                         size = n, replace = TRUE,
                         prob = 1/edge_df$weight * edge_df$probs)

  sel_edges <- edge_df[sel_edges_id,]
  sel_edges_len <- sel_edges$weight

  # preparing some variables for later
  start_oids <- sel_edges$start_oid
  end_oids <- sel_edges$end_oid
  all_names <- names(igraph::V(graph))

  # finding the start nodes and end nodes of the selected edges
  start_names <- all_names[start_oids]
  end_names <- all_names[end_oids]

  #b. calculating the position of the point on the edges
  # each edge will receive one point
  dists <- vec_runif(n=1,min = 0, max = sel_edges_len)
  # creating virtual names for the new vertices
  new_vert <- paste0(rep("virt_"),seq_len(length(dists)))

  #c. creating the new edges and nodes as another graph
  df <- data.frame(start = c(start_names,new_vert),
                   end = c(new_vert,end_names),
                   weight = c(dists,(sel_edges_len - dists)))

  #d. and then merging the graphs
  new_graph <- igraph::graph_from_data_frame(df, directed = FALSE)
  tot_graph <- igraph::union(graph,new_graph, byname = TRUE)

  # We just need to merge the weights of the graphs
  ws <- igraph::E(tot_graph)
  df_tmp <- data.frame("w1" = ws$weight_1,
                       "w2" = ws$weight_2)

  df_tmp[is.na(df_tmp$w1),"w1"] <- 0
  df_tmp[is.na(df_tmp$w2),"w2"] <- 0
  tot_graph <- igraph::set_edge_attr(tot_graph, "weight",
                                     value = df_tmp$w1 + df_tmp$w2,
                                     index = igraph::E(tot_graph))


  # calculating the distances
  if (is.null(start_vert)){
    dist_mat <- igraph::distances(tot_graph,v = new_vert, to = new_vert)


  }else{
    dist_mat <- igraph::distances(tot_graph,v = start_vert, to = new_vert)

  }

  return(dist_mat)
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
#' @return A list with the following values : \cr \itemize{ \item{plotk}{ A
#'   ggplot2 object representing the values of the k-function} \item{plotg}{ A
#'   ggplot2 object representing the values of the g-function} \item{values}{ A
#'   DataFrame with the values used to build the plots} }
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot geom_ribbon geom_path aes_string labs
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' main_network_mtl <- sf::st_read(networkgpkg,layer="main_network_mtl")
#' mtl_libraries <- sf::st_read(eventsgpkg,layer="mtl_libraries")
#' result <- kfunctions(main_network_mtl, mtl_libraries,
#'      start = 0, end = 2500, step = 10,
#'      width = 200, nsim = 50,
#'      conf_int = 0.05, tol = 0.1, agg = NULL,
#'      verbose = FALSE)
#' }
kfunctions <- function(lines, points,
                       start, end, step, width,
                       nsim, conf_int = 0.05,
                       digits = 2, tol = 0.1,
                       resolution = NULL, agg = NULL,
                       verbose = TRUE, return_sims = FALSE){

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
  ## step6 : calcualte the kfunction and the g function
  if (verbose){
    print("Calculating k and g functions ...")
  }
  k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = snapped_events$weight)
  g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = snapped_events$weight)

  ## step7 : generate the permutations
  if (verbose){
    print("Calculating the simulations ...")
  }
  w <- rep(1,times = n)
  if(verbose){
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }

  # the case where we can simplified the situation
  if (is.null(resolution)==FALSE){
    dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
                                           edge_df = graph_result$spedges,
                                           n = n,
                                           resolution = resolution,
                                           nsim = nsim)

    all_values <- lapply(1:nsim,function(i){
      dist_mat <- dist_matrices[[i]]
      k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
      g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
      if(verbose){
        setTxtProgressBar(pb, i)
      }
      return(cbind(k_vals,g_vals))
    })

  }else{
    # the case where we can not simplified the situation
    all_values <- lapply(1:nsim,function(i){
      dist_mat <- randomize_distmatrix(graph_result$graph,graph_result$spedges,n)
      k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
      g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
      if(verbose){
        setTxtProgressBar(pb, i)
      }
      return(cbind(k_vals,g_vals))
    })
  }

  ## step8 : extract the k_vals and g_vals matrices
  k_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,1])}))
  g_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,2])}))

  ## step9 : calculating the summary stats
  upper <- 1-conf_int / 2
  lower <- conf_int / 2
  k_stats <- apply(k_mat,MARGIN = 1, function(i){
    return(quantile(i,probs = c(lower,upper)))
  })
  g_stats <- apply(g_mat,MARGIN = 1, function(i){
    return(quantile(i,probs = c(lower,upper)))
  })

  plot_df <- data.frame(
    "obs_k" = k_vals,
    "lower_k" = k_stats[1,],
    "upper_k" = k_stats[2,],
    "obs_g" = g_vals,
    "lower_g" = g_stats[1,],
    "upper_g" = g_stats[2,],
    "distances" = seq(start,end,step)
  )

  plotk <- ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_k",ymax = "upper_k"),
                fill = grDevices::rgb(0.1,0.1,0.1), alpha=0.4)+
    geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
    labs(x = "distances",
         y = "empirical K-function")

  plotg <- ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_g", ymax = "upper_g"),
                fill = grDevices::rgb(0.1,0.1,0.1),alpha=0.4, )+
    geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
    labs(x = "distances",
         y = "empirical G-function")

  obj <- list(
    "plotk" = plotk,
    "plotg" = plotg,
    "values" = plot_df
  )
  if(return_sims){
    obj$sim_k_values <- k_mat
    obj$sim_g_values <- g_mat
  }

  return(obj)
}


#' @title Network k and g functions (multicore, maturing)
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
#'
#' @return A list with the following values : \cr \itemize{ \item{plotk}{ A
#'   ggplot2 object representing the values of the k-function} \item{plotg}{ A
#'   ggplot2 object representing the values of the g-function} \item{values}{ A
#'   DataFrame with the values used to build the plots} }
#'
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_path labs aes_string
#' @importFrom igraph E
#'
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' main_network_mtl <-  sf::st_read(networkgpkg,layer="main_network_mtl")
#' mtl_libraries <- sf::st_read(eventsgpkg,layer="mtl_libraries")
#' future::plan(future::multisession(workers=2))
#' result <- kfunctions.mc(main_network_mtl, mtl_libraries,
#'      start = 0, end = 2500, step = 10,
#'      width = 200, nsim = 50,
#'      conf_int = 0.05, tol = 0.1, agg = NULL,
#'      verbose = FALSE)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
kfunctions.mc <- function(lines, points, start, end, step, width, nsim, conf_int = 0.05,
                          digits = 2 ,tol = 0.1, resolution = 50, agg = NULL,
                          verbose = TRUE, return_sims = FALSE){

  ## step0 : clean the points
  if(verbose){
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
  if(verbose){
    print("Snapping points on lines ...")
  }
  snapped_events <- snapPointsToLines2(points,lines,idField = "oid")
  new_lines <- split_lines_at_vertex(lines, snapped_events,
                                     snapped_events$nearest_line_id, tol)
  #new_lines <- add_vertices_lines(lines,snapped_events,snapped_events$nearest_line_id,tol)

  ## step3 : splitting the lines
  if(verbose){
    print("Building graph ...")
  }
  #new_lines <- simple_lines(new_lines)
  new_lines$length <- as.numeric(st_length(new_lines))
  new_lines <- subset(new_lines,new_lines$length>0)
  new_lines <- remove_loop_lines(new_lines,digits)
  new_lines$oid <- seq_len(nrow(new_lines))
  new_lines <- new_lines[c("length","oid","probs")]
  Lt <- sum(as.numeric(st_length(new_lines)))

  ## step4 : building the graph for the real case
  graph_result <- build_graph(new_lines,digits = digits,line_weight = "length", attrs = TRUE)
  graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
                                                           name = "probs")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  snapped_events$vertex_id <- closest_points(snapped_events, nodes)

  ## step5 : calculating the distance matrix
  dist_mat <- igraph::distances(graph,v = snapped_events$vertex_id,
                                to = snapped_events$vertex_id)
  ## step6 : calcualte the kfunction and the g function
  if(verbose){
    print("Calculating k and g functions ...")
  }
  k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = snapped_events$weight)
  g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = snapped_events$weight)

  ## step7 : generate the permutations
  if(verbose){
    print("Calculating the simulations ...")
  }

  w <- rep(1,times = n)
  sim_seq <- 1:nsim
  graph <- graph_result$graph
  edgesdf <- st_drop_geometry(graph_result$spedges)

  # the classical way
  if (is.null(resolution)){
    if(verbose){
      progressr::with_progress({
        p <- progressr::progressor(along = sim_seq)
        all_values <- future.apply::future_lapply(sim_seq, function(i){
          dist_mat <- randomize_distmatrix(graph, edgesdf, n)
          k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
          g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
          return(cbind(k_vals,g_vals))
        },future.packages = c("igraph"))
      })
    }else{
      all_values <- future.apply::future_lapply(sim_seq, function(i){
        dist_mat <- randomize_distmatrix(graph, edgesdf, n)
        k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
        g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
        return(cbind(k_vals,g_vals))
      },future.packages = c("igraph"))

    }
  }else{
    # the simplified way
    ## first : generating the matrices
    if(verbose){
      print("generating the randomized distance matrices...")
    }
    dist_mats <- randomize_distmatrix2(graph, edgesdf,
                                      n = n, nsim = nsim,
                                      resolution = resolution)
    if(verbose){
      print("calculating the k and g functions for the randomized matrices..")
    }
    all_values <- future.apply::future_lapply(dist_mats, function(dist_mat){
      k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
      g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
      return(cbind(k_vals,g_vals))
    },future.packages = c("igraph"))

  }


  ## step8 : extract the k_vals and g_vals matrices
  k_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,1])}))
  g_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,2])}))

  ## step9 : calculating the summary stats
  lower <- 1-conf_int/2
  upper <- conf_int/2
  k_stats <- apply(k_mat,MARGIN = 1, function(i){
    return(quantile(i,probs = c(lower,upper)))
  })
  g_stats <- apply(g_mat,MARGIN = 1, function(i){
    return(quantile(i,probs = c(lower,upper)))
  })

  plot_df <- data.frame(
    "obs_k" = k_vals,
    "upper_k" = k_stats[1,],
    "lower_k" = k_stats[2,],
    "obs_g" = g_vals,
    "upper_g" = g_stats[1,],
    "lower_g" = g_stats[2,],
    "distances" = seq(start,end,step)
  )

  plotk <- ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_k", ymax = "upper_k"),
                fill = grDevices::rgb(0.1,0.1,0.1),alpha=0.4, )+
    geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
    labs(x = "distances",
         y = "empirical K-function")

  plotg <- ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_g", ymax = "upper_g"),
                fill = grDevices::rgb(0.1,0.1,0.1),alpha=0.4, )+
    geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
    labs(x = "distances",
         y = "empirical G-function")

  obj <- list(
    "plotk" = plotk,
    "plotg" = plotg,
    "values" = plot_df
  )

  if(return_sims){
    obj$sim_k_values <- k_mat
    obj$sim_g_values <- g_mat
  }

  return(obj)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### execution cross-k functions ####
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
#' @return A list with the following values : \cr \itemize{ \item{plotk}{ A
#'   ggplot2 object representing the values of the cross k-function}
#'   \item{plotg}{ A ggplot2 object representing the values of the cross
#'   g-function} \item{values}{ A DataFrame with the values used to build the
#'   plots} }
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot geom_ribbon geom_path aes_string labs
#' @importFrom grDevices rgb
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' main_network_mtl <- sf::st_read(networkgpkg,layer="main_network_mtl")
#' mtl_libraries <- sf::st_read(eventsgpkg,layer="mtl_libraries")
#' mtl_theatres <- sf::st_read(eventsgpkg,layer="mtl_theatres")
#' result <- cross_kfunctions(main_network_mtl, mtl_theatres, mtl_libraries,
#'                            start = 0, end = 2500, step = 10, width = 250,
#'                            nsim = 50, conf_int = 0.05, digits = 2,
#'                            tol = 0.1, agg = NULL, verbose = FALSE)
#' }
cross_kfunctions <- function(lines, pointsA, pointsB,
                             start, end, step, width,
                             nsim, conf_int = 0.05,
                             digits = 2, tol = 0.1,
                             resolution = NULL, agg = NULL,
                             verbose = TRUE, return_sims = FALSE){

  ## step0 : clean the points
  if(verbose){
    print("Preparing data ...")
  }
  na <- nrow(pointsA)
  nb <- nrow(pointsB)

  probs <- NULL

  pointsA$weight <- rep(1,nrow(pointsA))
  pointsA <- clean_events(pointsA,digits,agg)
  pointsA$goid <- seq_len(nrow(pointsA))

  pointsB$weight <- rep(1,nrow(pointsB))
  pointsB <- clean_events(pointsB,digits,agg)
  pointsB$goid <- seq_len(nrow(pointsB))

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
  # new_lines <- add_vertices_lines(lines, snapped_events,
  #                                 snapped_events$nearest_line_id, tol)

  ## step3 : splitting the lines
  if(verbose){
    print("Building graph ...")
  }
  #new_lines <- simple_lines(new_lines)
  new_lines$length <- as.numeric(st_length(new_lines))
  new_lines <- subset(new_lines, new_lines$length>0)
  new_lines <- remove_loop_lines(new_lines, digits)
  new_lines$oid <- seq_len(nrow(new_lines))
  new_lines <- new_lines[c("length","oid","probs")]
  Lt <- sum(as.numeric(st_length(new_lines)))

  ## step4 : building the graph for the real case
  graph_result <- build_graph(new_lines,digits = digits,
                              line_weight = "length", attrs = TRUE)
  graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
                                                           name = "probs")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  snapped_events$vertex_id <- closest_points(snapped_events, nodes)
  snappedA <- subset(snapped_events, snapped_events$type == "A")
  snappedB <- subset(snapped_events, snapped_events$type == "B")

  ## step5 : calculating the distance matrix
  dist_mat <- igraph::distances(graph,v = snappedB$vertex_id,
                                to = snappedA$vertex_id, mode = "out")
  ## step6 : calcualte the kfunction and the g function
  if(verbose){
    print("Calculating k and g functions ...")
  }
  k_vals <- cross_kfunc_cpp(dist_mat,start,end,step,Lt,na,nb,snappedA$weight,snappedB$weight)
  g_vals <- cross_gfunc_cpp(dist_mat,start,end,step,width,Lt,na,nb,snappedA$weight,snappedB$weight)

  ## step7 : generate the permutations
  if(verbose){
    print("Calculating the simulations ...")
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }
  w <- rep(1,times = na)

  # the case where we can simplified the situation
  if (is.null(resolution)==FALSE){
    dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
                                           edge_df = graph_result$spedges,
                                           n = na,
                                           start_vert = snappedB$vertex_id,
                                           resolution = resolution,
                                           nsim = nsim)

    all_values <- lapply(1:nsim,function(i){
      dist_mat <- dist_matrices[[i]]
      k_vals <- cross_kfunc_cpp(dist_mat,start,end,step,Lt,na,nb,w,snappedB$weight)
      g_vals <- cross_gfunc_cpp(dist_mat,start,end,step,width,Lt,na,nb,w,snappedB$weight)
      if(verbose){
        setTxtProgressBar(pb, i)
      }
      return(cbind(k_vals,g_vals))
    })

  }else{
    # the case where we can not simplified the situation
    all_values <- lapply(1:nsim,function(i){
      dist_mat <- randomize_distmatrix(graph_result$graph,graph_result$spedges,
                                       na,start_vert = snappedB$vertex_id)
      k_vals <- cross_kfunc_cpp(dist_mat,start,end,step,Lt,na,nb,w,snappedB$weight)
      g_vals <- cross_gfunc_cpp(dist_mat,start,end,step,width,Lt,na,nb,w,snappedB$weight)
      if(verbose){
        setTxtProgressBar(pb, i)
      }
      return(cbind(k_vals,g_vals))
    })
  }


  ## step8 : extract the k_vals and g_vals matrices
  k_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,1])}))
  g_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,2])}))

  ## step9 : calculating the summary stats
  upper <- 1-conf_int / 2
  lower <- conf_int / 2
  k_stats <- apply(k_mat,MARGIN = 1, function(i){
    return(quantile(i,probs = c(lower,upper)))
  })
  g_stats <- apply(g_mat,MARGIN = 1, function(i){
    return(quantile(i,probs = c(lower,upper)))
  })

  plot_df <- data.frame(
    "obs_k" = k_vals,
    "lower_k" = k_stats[1,],
    "upper_k" = k_stats[2,],
    "obs_g" = g_vals,
    "lower_g" = g_stats[1,],
    "upper_g" = g_stats[2,],
    "distances" = seq(start,end,step)
  )

  plotk <- ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin = "lower_k", ymax = "upper_k"),
                fill = rgb(0.1,0.1,0.1),alpha=0.4, )+
    geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
    labs(x = "distances",
         y = "empirical cross-K-function")

  plotg <- ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin = "lower_g", ymax = "upper_g"),
                fill = rgb(0.1,0.1,0.1),alpha=0.4, )+
    geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
    labs(x = "distances",
         y = "empirical cross-G-function")

  obj <- list(
    "plotk" = plotk,
    "plotg" = plotg,
    "values" = plot_df
  )
  if(return_sims){
    obj$sim_k_values <- k_mat
    obj$sim_g_values <- g_mat
  }
  return(obj)
}


#' @title Network cross k and g functions (multicore, maturing)
#'
#' @description Calculate the cross k and g functions for a set of points on a
#'   network with multicore support. (maturing)
#'
#' @template kross_kfunctions-arg
#' @template common_kfunctions-arg
#' @param return_sims a boolean indicating if the simulated k and g values must also
#' be returned as matrices
#'
#' @return A list with the following values : \cr \itemize{ \item{plotk}{ A
#'   ggplot2 object representing the values of the cross k-function}
#'   \item{plotg}{ A ggplot2 object representing the values of the cross
#'   g-function} \item{values}{ A DataFrame with the values used to build the
#'   plots} }
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot geom_ribbon geom_path aes_string labs
#' @importFrom grDevices rgb
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' main_network_mtl <- sf::st_read(networkgpkg,layer="main_network_mtl")
#' mtl_libraries <- sf::st_read(eventsgpkg,layer="mtl_libraries")
#' mtl_theatres <- sf::st_read(eventsgpkg,layer="mtl_theatres")
#' future::plan(future::multisession(workers=2))
#' result <- cross_kfunctions.mc(main_network_mtl, mtl_libraries, mtl_theatres,
#'                            start = 0, end = 2500, step = 10, width = 250,
#'                            nsim = 50, conf_int = 0.05, digits = 2,
#'                            tol = 0.1, agg = NULL, verbose = TRUE)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
cross_kfunctions.mc <- function(lines, pointsA, pointsB,
                                start, end, step, width,
                                nsim, conf_int = 0.05,
                                digits = 2, tol = 0.1,
                                resolution = NULL, agg = NULL,
                                verbose = TRUE, return_sims = FALSE){

  ## step0 : clean the points
  if(verbose){
    print("Preparing data ...")
  }
  na <- nrow(pointsA)
  nb <- nrow(pointsB)

  probs <- NULL

  pointsA$weight <- rep(1,nrow(pointsA))
  pointsA <- clean_events(pointsA,digits,agg)
  pointsA$goid <- seq_len(nrow(pointsA))

  pointsB$weight <- rep(1,nrow(pointsB))
  pointsB <- clean_events(pointsB,digits,agg)
  pointsB$goid <- seq_len(nrow(pointsB))

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
  # new_lines <- add_vertices_lines(lines, snapped_events,
  #                                 snapped_events$nearest_line_id, tol)

  ## step3 : splitting the lines
  if(verbose){
    print("Building graph ...")
  }
  #new_lines <- simple_lines(new_lines)
  new_lines$length <- as.numeric(st_length(new_lines))
  new_lines <- subset(new_lines,new_lines$length>0)
  new_lines <- remove_loop_lines(new_lines,digits)
  new_lines$oid <- seq_len(nrow(new_lines))
  new_lines <- new_lines[c("length","oid","probs")]
  Lt <- sum(as.numeric(st_length(new_lines)))

  ## step4 : building the graph for the real case
  graph_result <- build_graph(new_lines,digits = digits,
                              line_weight = "length", attrs = TRUE)
  graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
                                                           name = "probs")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  snapped_events$vertex_id <- closest_points(snapped_events, nodes)
  snappedA <- subset(snapped_events, snapped_events$type == "A")
  snappedB <- subset(snapped_events, snapped_events$type == "B")

  ## step5 : calculating the distance matrix
  dist_mat <- igraph::distances(graph,v = snappedB$vertex_id,
                                to = snappedA$vertex_id)
  ## step6 : calcualte the kfunction and the g function
  if(verbose){
    print("Calculating k and g functions ...")
  }
  k_vals <- cross_kfunc_cpp(dist_mat, start, end, step, Lt, na, nb,
                        snappedA$weight, snappedB$weight)
  g_vals <- cross_gfunc_cpp(dist_mat, start, end, step, width, Lt, na,
                       nb, snappedA$weight, snappedB$weight)

  ## step7 : generate the permutations
  if(verbose){
    print("Calculating the simulations ...")
  }
  w <- rep(1,times = na)
  sim_seq <- 1:nsim

  # classical approach
  if (is.null(resolution)){
    progressr::with_progress({
      p <- progressr::progressor(along = sim_seq)
      all_values <- future.apply::future_lapply(sim_seq,function(i){
        dist_mat <- randomize_distmatrix(graph_result$graph,
                                         graph_result$spedges,
                                         na, start_vert = snappedB$vertex_id)
        k_vals <- cross_kfunc_cpp(dist_mat, start, end, step, Lt, na, nb,
                              w, snappedB$weight)
        g_vals <- cross_gfunc_cpp(dist_mat, start, end, step, width, Lt, na, nb,
                             w,snappedB$weight)
        return(cbind(k_vals,g_vals))
      },future.packages = c("igraph","base"))
    })
  }else{
    #simplified approach
    if(verbose){
      print("calculating the randomized distance matrices...")
    }
    dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
                                           edge_df = graph_result$spedges,
                                           n = na,
                                           start_vert = snappedB$vertex_id,
                                           resolution = resolution,
                                           nsim = nsim)
    if(verbose){
      print("calculating the k and g functions for the randomized matrices..")
    }
    all_values <- future.apply::future_lapply(dist_matrices,function(dist_mat){
        k_vals <- cross_kfunc_cpp(dist_mat, start, end, step, Lt, na, nb,
                              w, snappedB$weight)
        g_vals <- cross_gfunc_cpp(dist_mat, start, end, step, width, Lt, na, nb,
                             w,snappedB$weight)
        return(cbind(k_vals,g_vals))
      },future.packages = c("igraph","base"))
  }


  ## step8 : extract the k_vals and g_vals matrices
  k_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,1])}))
  g_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,2])}))

  ## step9 : calculating the summary stats
  upper <- 1-conf_int / 2
  lower <- conf_int / 2
  k_stats <- apply(k_mat,MARGIN = 1, function(i){
    return(quantile(i,probs = c(lower,upper)))
  })
  g_stats <- apply(g_mat,MARGIN = 1, function(i){
    return(quantile(i,probs = c(lower,upper)))
  })

  plot_df <- data.frame(
    "obs_k" = k_vals,
    "lower_k" = k_stats[1,],
    "upper_k" = k_stats[2,],
    "obs_g" = g_vals,
    "lower_g" = g_stats[1,],
    "upper_g" = g_stats[2,],
    "distances" = seq(start,end,step)
  )

  plotk <- ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin = "lower_k", ymax = "upper_k"),
                fill = rgb(0.1,0.1,0.1),alpha=0.4, )+
    geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
    labs(x = "distances",
         y = "empirical cross-K-function")

  plotg <- ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin = "lower_g", ymax = "upper_g"),
                fill = rgb(0.1,0.1,0.1),alpha=0.4, )+
    geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
    labs(x = "distances",
         y = "empirical cross-G-function")

  obj <- list(
    "plotk" = plotk,
    "plotg" = plotg,
    "values" = plot_df
  )
  if(return_sims){
    obj$sim_k_values <- k_mat
    obj$sim_g_values <- g_mat
  }
  return(obj)
}

