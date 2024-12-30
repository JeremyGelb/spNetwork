# k_nt_functions <- function(lines, points, points_time,
#                            start_net, end_net, step_net, width_net,
#                            start_time, end_time, step_time, width_time,
#                            nsim, conf_int = 0.05, digits = 2, tol = 0.1,
#                            resolution = NULL, agg = NULL, verbose = TRUE){
#
#   ## step0 : clean the points
#   if (verbose){
#     print("Preparing data ...")
#   }
#   n <- nrow(points)
#   points$goid <- seq_len(nrow(points))
#   points$weight <- rep(1,nrow(points))
#
#   # we aggregate the points to have a fewer number of locations
#   # but because of the time dimension, we have to keep them separated in the end
#   agg_points <- clean_events(points,digits,agg)
#   agg_points$locid <- 1:nrow(agg_points)
#
#   # and we need now to find the new location of the original points
#   xy_origin <- st_coordinates(points)
#   xy_agg <- st_coordinates(agg_points)
#   ids <- dbscan::kNN(xy_agg, k = 1, query = xy_origin)
#   points$locid <- ids$id[,1]
#
#   probs <- NULL
#
#   ## step1 : clean the lines
#   if(is.null(probs)){
#     lines$probs <- 1
#   }else{
#     lines$probs <- probs
#   }
#
#   lines$length <- as.numeric(st_length(lines))
#   lines <- subset(lines, lines$length>0)
#   lines$oid <- seq_len(nrow(lines))
#
#
#   ## step2 : adding the points to the lines
#   if (verbose){
#     print("Snapping points on lines ...")
#   }
#   snapped_events <- snapPointsToLines2(agg_points, lines, idField = "oid")
#   new_lines <- split_lines_at_vertex(lines, snapped_events,
#                                      snapped_events$nearest_line_id, tol)
#
#   ## step3 : splitting the lines
#   if (verbose){
#     print("Building graph ...")
#   }
#   #new_lines <- simple_lines(new_lines)
#   new_lines$length <- as.numeric(st_length(new_lines))
#   new_lines <- subset(new_lines,new_lines$length>0)
#
#   new_lines <- remove_loop_lines(new_lines,digits)
#
#   new_lines$oid <- seq_len(nrow(new_lines))
#   new_lines <- new_lines[c("length","oid","probs")]
#
#   # calculating the extent of the study area
#   Lt <- sum(as.numeric(st_length(new_lines)))
#   Tt <- max(points_time) - min(points_time)
#
#   new_lines$weight <- as.numeric(st_length(new_lines))
#
#   ## step4 : building the graph for the real case
#   graph_result <- build_graph(new_lines,digits = digits,
#                               line_weight = "weight",
#                               attrs = TRUE)
#
#   graph <- graph_result$graph
#   nodes <- graph_result$spvertices
#   graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
#                                                            name = "probs")
#   snapped_events$vertex_id <- closest_points(snapped_events, nodes)
#
#   ## step5 : calculating the distance matrix
#   dist_mat <- igraph::distances(graph,v = snapped_events$vertex_id,
#                                 to = snapped_events$vertex_id)
#
#   ## step 5.5 I must now deal with the dupplicated ids
#   ## minus 1 because c++ indexing starts at 0
#   dist_mat_net <- extend_matrix_by_ids(dist_mat, points$goid, points$locid-1)
#
#   ## and generate a matrix with the time distances !
#   dist_mat_time <- as.matrix(stats::dist(points_time))
#
#   ## step6 : calcualte the kfunction and the g function
#   if (verbose){
#     print("Calculating k and g functions ...")
#   }
#   k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time,
#                           start_net, end_net, step_net,
#                           start_time, end_time, step_time,
#                           Lt, Tt, n, w = points$weight)
#
#   g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time,
#                           start_net, end_net, step_net, width_net,
#                           start_time, end_time, step_time, width_time,
#                           Lt, Tt, n,w = points$weight)
#
#   ## step7 : generate the permutations
#   if (verbose){
#     print("Calculating the simulations ...")
#   }
#   w <- rep(1,times = n)
#   if(verbose){
#     pb <- txtProgressBar(min = 0, max = nsim, style = 3)
#   }
#
#   dims_time <- dim(dist_mat_time)
#
#   # the case where we can simplified the situation
#   if (is.null(resolution)==FALSE){
#     dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
#                                            edge_df = graph_result$spedges,
#                                            n = n,
#                                            resolution = resolution,
#                                            nsim = nsim)
#
#     all_values <- lapply(1:nsim,function(i){
#       dist_mat_net <- dist_matrices[[i]]
#
#       dist_mat_time <- runif(n = prod(dims_time),
#                              min = min(points_time), max = max(points_time))
#       dim(dist_mat_time) <- dims_time
#
#       k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time,
#                               start_net, end_net, step_net,
#                               start_time, end_time, step_time,
#                               Lt, Tt, n, w = w)
#
#       g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time,
#                               start_net, end_net, step_net, width_net,
#                               start_time, end_time, step_time, width_time,
#                               Lt, Tt, n, w = w)
#       if(verbose){
#         setTxtProgressBar(pb, i)
#       }
#       return(list(k_vals,g_vals))
#     })
#
#   }else{
#     # the case where we can not simplified the situation
#     all_values <- lapply(1:nsim,function(i){
#       dist_mat_net <- randomize_distmatrix(graph_result$graph,graph_result$spedges,n)
#       dist_mat_time <- runif(n = prod(dims_time),
#                              min = min(points_time), max = max(points_time))
#       dim(dist_mat_time) <- dims_time
#
#       k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time,
#                               start_net, end_net, step_net,
#                               start_time, end_time, step_time,
#                               Lt, Tt, n, w = w)
#
#       g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time,
#                               start_net, end_net, step_net, width_net,
#                               start_time, end_time, step_time, width_time,
#                               Lt, Tt, n, w = w)
#       if(verbose){
#         setTxtProgressBar(pb, i)
#       }
#       return(list(k_vals,g_vals))
#     })
#   }
#
#   ## step8 : extract the k_vals and g_vals matrices
#   L1 <- lapply(all_values,function(i){return(i[[1]])})
#   L2 <- lapply(all_values,function(i){return(i[[2]])})
#   L1[["along"]] <- 3
#   L2[["along"]] <- 3
#   k_mat <- do.call(abind::abind, L1)
#   g_mat <- do.call(abind::abind, L2)
#
#   ## step9 : calculating the summary stats
#   upper <- 1 - conf_int / 2
#   lower <- conf_int / 2
#   k_stats_lower <- apply(k_mat,MARGIN = c(1,2), function(i){
#     return(quantile(i,probs = lower))
#   })
#   k_stats_upper <- apply(k_mat,MARGIN = c(1,2), function(i){
#     return(quantile(i,probs = upper))
#   })
#   g_stats_upper <- apply(g_mat,MARGIN = c(1,2), function(i){
#     return(quantile(i,probs = upper))
#   })
#   g_stats_lower <- apply(g_mat,MARGIN = c(1,2), function(i){
#     return(quantile(i,probs = lower))
#   })
#
#   seq_net <- seq_num2(start_net, end_net, step_net)
#   seq_time <- seq_num2(start_time, end_time, step_time)
#
#   row.names(k_vals) <- seq_net
#   row.names(k_stats_lower) <- seq_net
#   row.names(k_stats_upper) <- seq_net
#   row.names(g_vals) <- seq_net
#   row.names(g_stats_lower) <- seq_net
#   row.names(g_stats_upper) <- seq_net
#
#   colnames(k_vals) <- seq_time
#   colnames(k_stats_lower) <- seq_time
#   colnames(k_stats_upper) <- seq_time
#   colnames(g_vals) <- seq_time
#   colnames(g_stats_lower) <- seq_time
#   colnames(g_stats_upper) <- seq_time
#
#   results <- list(
#     "obs_k" = k_vals,
#     "lower_k" = k_stats_lower,
#     "upper_k" = k_stats_upper,
#     "obs_g" = g_vals,
#     "lower_g" = g_stats_lower,
#     "upper_g" = g_stats_upper,
#     "distances_net" = seq_num2(start_net, end_net, step_net),
#     "distances_time" = seq_num2(start_time, end_time, step_time)
#   )
#   # see this plot https://www.researchgate.net/publication/320760197_Spatiotemporal_Point_Pattern_Analysis_Using_Ripley%27s_K_Function
#   return(results)
# }



# k_nt_functions.mc <- function(lines, points, points_time,
#                            start_net, end_net, step_net, width_net,
#                            start_time, end_time, step_time, width_time,
#                            nsim, conf_int = 0.05, digits = 2, tol = 0.1,
#                            resolution = NULL, agg = NULL, verbose = TRUE){
#
#   ## step0 : clean the points
#   if (verbose){
#     print("Preparing data ...")
#   }
#   n <- nrow(points)
#   points$goid <- seq_len(nrow(points))
#   points$weight <- rep(1,nrow(points))
#
#   # we aggregate the points to have a fewer number of locations
#   # but because of the time dimension, we have to keep them separated in the end
#   agg_points <- clean_events(points,digits,agg)
#   agg_points$locid <- 1:nrow(agg_points)
#
#   # and we need now to find the new location of the original points
#   xy_origin <- st_coordinates(points)
#   xy_agg <- st_coordinates(agg_points)
#   ids <- dbscan::kNN(xy_agg, k = 1, query = xy_origin)
#   points$locid <- ids$id[,1]
#
#   probs <- NULL
#
#   ## step1 : clean the lines
#   if(is.null(probs)){
#     lines$probs <- 1
#   }else{
#     lines$probs <- probs
#   }
#
#   lines$length <- as.numeric(st_length(lines))
#   lines <- subset(lines, lines$length>0)
#   lines$oid <- seq_len(nrow(lines))
#
#
#   ## step2 : adding the points to the lines
#   if (verbose){
#     print("Snapping points on lines ...")
#   }
#   snapped_events <- snapPointsToLines2(agg_points, lines, idField = "oid")
#   new_lines <- split_lines_at_vertex(lines, snapped_events,
#                                      snapped_events$nearest_line_id, tol)
#
#   ## step3 : splitting the lines
#   if (verbose){
#     print("Building graph ...")
#   }
#   #new_lines <- simple_lines(new_lines)
#   new_lines$length <- as.numeric(st_length(new_lines))
#   new_lines <- subset(new_lines,new_lines$length>0)
#
#   new_lines <- remove_loop_lines(new_lines,digits)
#
#   new_lines$oid <- seq_len(nrow(new_lines))
#   new_lines <- new_lines[c("length","oid","probs")]
#   Lt <- sum(as.numeric(st_length(new_lines)))
#   Tt <- max(points_time) - min(points_time)
#
#   new_lines$weight <- as.numeric(st_length(new_lines))
#
#   ## step4 : building the graph for the real case
#   graph_result <- build_graph(new_lines,digits = digits,
#                               line_weight = "weight",
#                               attrs = TRUE)
#   graph <- graph_result$graph
#   nodes <- graph_result$spvertices
#   graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
#                                                            name = "probs")
#   snapped_events$vertex_id <- closest_points(snapped_events, nodes)
#
#   ## step5 : calculating the distance matrix
#   dist_mat <- igraph::distances(graph,v = snapped_events$vertex_id,
#                                 to = snapped_events$vertex_id)
#
#   ## step 5.5 I must now deal with the dupplicated ids
#   ## minus 1 because c++ indexing starts at 0
#   dist_mat_net <- extend_matrix_by_ids(dist_mat, points$goid, points$locid-1)
#
#   ## and generate a matrix with the time distances !
#   dist_mat_time <- as.matrix(stats::dist(points_time))
#
#   ## step6 : calcualte the kfunction and the g function
#   if (verbose){
#     print("Calculating k and g functions ...")
#   }
#
#   k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time,
#                           start_net, end_net, step_net,
#                           start_time, end_time, step_time,
#                           Lt, Tt, n, w = points$weight)
#
#   g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time,
#                           start_net, end_net, step_net, width_net,
#                           start_time, end_time, step_time, width_time,
#                           Lt, Tt, n,w = points$weight)
#
#   ## step7 : generate the permutations
#   if (verbose){
#     print("Calculating the simulations ...")
#   }
#   w <- rep(1,times = n)
#   if(verbose){
#     pb <- txtProgressBar(min = 0, max = nsim, style = 3)
#   }
#
#   dims_time <- dim(dist_mat_time)
#   sim_seq <- 1:nsim
#   graph <- graph_result$graph
#   edgesdf <- st_drop_geometry(graph_result$spedges)
#
#   # the classical way
#   if (is.null(resolution)){
#     if(verbose){
#       progressr::with_progress({
#         p <- progressr::progressor(along = sim_seq)
#         all_values <- future.apply::future_lapply(sim_seq, function(i){
#
#           dist_mat_net <- randomize_distmatrix(graph_result$graph,graph_result$spedges,n)
#           dist_mat_time2 <- runif(n = prod(dims_time),
#                                  min = min(points_time), max = max(points_time))
#           dim(dist_mat_time2) <- dims_time
#
#           k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time2,
#                                   start_net, end_net, step_net,
#                                   start_time, end_time, step_time,
#                                   Lt, Tt, n, w = w)
#
#           g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time2,
#                                   start_net, end_net, step_net, width_net,
#                                   start_time, end_time, step_time, width_time,
#                                   Lt, Tt, n, w = w)
#           return(list(k_vals,g_vals))
#         },future.packages = c("igraph"))
#       })
#     }else{
#       all_values <- future.apply::future_lapply(sim_seq, function(i){
#         dist_mat_net <- randomize_distmatrix(graph_result$graph,graph_result$spedges,n)
#         dist_mat_time2 <- runif(n = prod(dims_time),
#                                 min = min(points_time), max = max(points_time))
#         dim(dist_mat_time2) <- dims_time
#
#         k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time2,
#                                 start_net, end_net, step_net,
#                                 start_time, end_time, step_time,
#                                 Lt, Tt, n, w = w)
#
#         g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time2,
#                                 start_net, end_net, step_net, width_net,
#                                 start_time, end_time, step_time, width_time,
#                                 Lt, Tt, n, w = w)
#         p(sprintf("i=%g", i))
#         return(list(k_vals,g_vals))
#       },future.packages = c("igraph"))
#
#     }
#   }else{
#     # the simplified way
#     ## first : generating the matrices
#     if(verbose){
#       print("generating the randomized distance matrices...")
#     }
#     dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
#                                            edge_df = graph_result$spedges,
#                                            n = n,
#                                            resolution = resolution,
#                                            nsim = nsim)
#     if(verbose){
#       print("calculating the k and g functions for the randomized matrices..")
#     }
#     all_values <- future.apply::future_lapply(dist_matrices, function(dist_mat_net){
#       dist_mat_time2 <- runif(n = prod(dims_time),
#                              min = min(points_time), max = max(points_time))
#       dim(dist_mat_time2) <- dims_time
#
#       k_vals <- k_nt_func_cpp(dist_mat_net, dist_mat_time2,
#                               start_net, end_net, step_net,
#                               start_time, end_time, step_time,
#                               Lt, Tt, n, w = w)
#
#       g_vals <- g_nt_func_cpp(dist_mat_net, dist_mat_time2,
#                               start_net, end_net, step_net, width_net,
#                               start_time, end_time, step_time, width_time,
#                               Lt, Tt, n, w = w)
#       return(list(k_vals,g_vals))
#     },future.packages = c("igraph"))
#
#   }
#
#   ## step8 : extract the k_vals and g_vals matrices
#   L1 <- lapply(all_values,function(i){return(i[[1]])})
#   L2 <- lapply(all_values,function(i){return(i[[2]])})
#   L1[["along"]] <- 3
#   L2[["along"]] <- 3
#   k_mat <- do.call(abind::abind, L1)
#   g_mat <- do.call(abind::abind, L2)
#
#   ## step9 : calculating the summary stats
#   upper <- 1 - conf_int / 2
#   lower <- conf_int / 2
#   k_stats_lower <- apply(k_mat,MARGIN = c(1,2), function(i){
#     return(quantile(i,probs = lower))
#   })
#   k_stats_upper <- apply(k_mat,MARGIN = c(1,2), function(i){
#     return(quantile(i,probs = upper))
#   })
#   g_stats_upper <- apply(g_mat,MARGIN = c(1,2), function(i){
#     return(quantile(i,probs = upper))
#   })
#   g_stats_lower <- apply(g_mat,MARGIN = c(1,2), function(i){
#     return(quantile(i,probs = lower))
#   })
#
#   seq_net <- seq_num2(start_net, end_net, step_net)
#   seq_time <- seq_num2(start_time, end_time, step_time)
#
#   row.names(k_vals) <- seq_net
#   row.names(k_stats_lower) <- seq_net
#   row.names(k_stats_upper) <- seq_net
#   row.names(g_vals) <- seq_net
#   row.names(g_stats_lower) <- seq_net
#   row.names(g_stats_upper) <- seq_net
#
#   colnames(k_vals) <- seq_time
#   colnames(k_stats_lower) <- seq_time
#   colnames(k_stats_upper) <- seq_time
#   colnames(g_vals) <- seq_time
#   colnames(g_stats_lower) <- seq_time
#   colnames(g_stats_upper) <- seq_time
#
#   results <- list(
#     "obs_k" = k_vals,
#     "lower_k" = k_stats_lower,
#     "upper_k" = k_stats_upper,
#     "obs_g" = g_vals,
#     "lower_g" = g_stats_lower,
#     "upper_g" = g_stats_upper,
#     "distances_net" = seq_num2(start_net, end_net, step_net),
#     "distances_time" = seq_num2(start_time, end_time, step_time)
#   )
#   # see this plot https://www.researchgate.net/publication/320760197_Spatiotemporal_Point_Pattern_Analysis_Using_Ripley%27s_K_Function
#   return(results)
# }

