# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### randomization functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' randomize_distmatrix2 <- function(graph, edge_df, n, resolution, nsim, start_vert = NULL){
#'   ## step1 : generate all the candidate nodes on the graph
#'   ## a : select all the edges with a length superior to resolution
#'   edit_edge <- subset(edge_df, edge_df$weight > resolution)
#'   all_names <- names(igraph::V(graph))
#'
#'   ## step2 : create the new needed vertices and edges
#'   vertices_and_distances <- lapply(1:nrow(edit_edge), function(i){
#'     this_edge <- edit_edge[i,]
#'     start_node <- all_names[this_edge$start_oid]
#'     end_node <- all_names[this_edge$end_oid]
#'     dists <- rep(resolution, floor(this_edge$weight/resolution))
#'     if(sum(dists) == this_edge$weight){
#'       dists <- dists[1:(length(dists)-1)]
#'     }
#'     names(dists) <- paste("fict",i,1:length(dists), sep="_")
#'     return(data.frame(starts = c(start_node, names(dists)),
#'                       ends = c(names(dists), end_node),
#'                       weight = c(dists, this_edge$weight - sum(dists))))
#'   })
#'
#'   all_elements <- do.call(rbind, vertices_and_distances)
#'   all_elements <- subset(all_elements,all_elements$weight>0)
#'
#'   ## step 3 : creating a new graph with the new vertices and edges
#'   new_graph <- igraph::graph_from_data_frame(all_elements, directed = FALSE)
#'
#'   ## then merging it with the previous graph
#'   tot_graph <- igraph::union(graph,new_graph, byname = TRUE)
#'
#'   ## step 4 correcting the weights
#'   ws <- igraph::E(tot_graph)
#'   df_tmp <- data.frame("w1" = ws$weight_1,
#'                        "w2" = ws$weight_2)
#'
#'   df_tmp[is.na(df_tmp$w1),"w1"] <- 0
#'   df_tmp[is.na(df_tmp$w2),"w2"] <- 0
#'   tot_graph <- igraph::set_edge_attr(tot_graph, "weight",
#'                                      value = df_tmp$w1 + df_tmp$w2,
#'                                      index = igraph::E(tot_graph))
#'
#'
#'   ## and now, calculating the distance matrices
#'   verts <- as.numeric(igraph::V(tot_graph))
#'   dist_matrices <- lapply(1:nsim, function(i){
#'     new_vert <- sample(verts,size = n, replace = FALSE)
#'     if (is.null(start_vert)){
#'       dist_mat <- igraph::distances(tot_graph,v = new_vert, to = new_vert)
#'     }else{
#'       dist_mat <- igraph::distances(tot_graph,v = start_vert, to = new_vert)
#'
#'     }
#'   })
#'   return(dist_matrices)
#' }


# randomize_distmatrix <- function(graph, edge_df, n, start_vert = NULL){
#
#   vec_runif <- Vectorize(runif, vectorize.args = c("max"))
#
#   ## prefered case where points are randomly located on edges
#
#   #a. selecting the edges that will have points
#   sel_edges_id <- sample(edge_df$edge_id,
#                          size = n, replace = TRUE,
#                          prob = 1/edge_df$weight * edge_df$probs)
#
#   sel_edges <- edge_df[sel_edges_id,]
#   sel_edges_len <- sel_edges$weight
#
#   # preparing some variables for later
#   start_oids <- sel_edges$start_oid
#   end_oids <- sel_edges$end_oid
#   all_names <- names(igraph::V(graph))
#
#   # finding the start nodes and end nodes of the selected edges
#   start_names <- all_names[start_oids]
#   end_names <- all_names[end_oids]
#
#   #b. calculating the position of the point on the edges
#   # each edge will receive one point
#   dists <- vec_runif(n=1,min = 0, max = sel_edges_len)
#   # creating virtual names for the new vertices
#   new_vert <- paste0(rep("virt_"),seq_len(length(dists)))
#
#   #c. creating the new edges and nodes as another graph
#   df <- data.frame(start = c(start_names,new_vert),
#                    end = c(new_vert,end_names),
#                    weight = c(dists,(sel_edges_len - dists)))
#
#   #d. and then merging the graphs
#   new_graph <- igraph::graph_from_data_frame(df, directed = FALSE)
#   tot_graph <- igraph::union(graph,new_graph, byname = TRUE)
#
#   # We just need to merge the weights of the graphs
#   ws <- igraph::E(tot_graph)
#   df_tmp <- data.frame("w1" = ws$weight_1,
#                        "w2" = ws$weight_2)
#
#   df_tmp[is.na(df_tmp$w1),"w1"] <- 0
#   df_tmp[is.na(df_tmp$w2),"w2"] <- 0
#   tot_graph <- igraph::set_edge_attr(tot_graph, "weight",
#                                      value = df_tmp$w1 + df_tmp$w2,
#                                      index = igraph::E(tot_graph))
#
#
#   # calculating the distances
#   if (is.null(start_vert)){
#     dist_mat <- igraph::distances(tot_graph,v = new_vert, to = new_vert)
#
#
#   }else{
#     dist_mat <- igraph::distances(tot_graph,v = start_vert, to = new_vert)
#
#   }
#
#   return(dist_mat)
# }


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### execution k functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# kfunctions <- function(lines, points,
#                        start, end, step, width,
#                        nsim, conf_int = 0.05,
#                        digits = 2, tol = 0.1,
#                        resolution = NULL, agg = NULL,
#                        verbose = TRUE, return_sims = FALSE){
#
#   ## step0 : clean the points
#   if (verbose){
#     print("Preparing data ...")
#   }
#   n <- nrow(points)
#   points$goid <- seq_len(nrow(points))
#   points$weight <- rep(1,nrow(points))
#   points <- clean_events(points,digits,agg)
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
#   snapped_events <- snapPointsToLines2(points,lines,idField = "oid")
#   new_lines <- split_lines_at_vertex(lines, snapped_events,
#                                      snapped_events$nearest_line_id, tol)
#
#   # new_lines <- add_vertices_lines(lines,snapped_events,
#   #                                 snapped_events$nearest_line_id, tol)
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
#   ## HERE, I SHOULD CHECK THAT TWO POINTS DO NOT SHARE THE SAME vertex
#   if(max(table(snapped_events$vertex_id))>1){
#     stop("After snapping the points on the network, some of them share the same location.
#          To correct it, please consider setting or increasing the value of the parameter agg.
#          They will be merged and their weights added")
#   }
#
#   ## step5 : calculating the distance matrix
#   dist_mat <- igraph::distances(graph,v = snapped_events$vertex_id,
#                                 to = snapped_events$vertex_id)
#   ## step6 : calcualte the kfunction and the g function
#   if (verbose){
#     print("Calculating k and g functions ...")
#   }
#   k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = snapped_events$weight)
#   g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = snapped_events$weight)
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
#   # the case where we can simplified the situation
#   if (is.null(resolution)==FALSE){
#     dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
#                                            edge_df = graph_result$spedges,
#                                            n = n,
#                                            resolution = resolution,
#                                            nsim = nsim)
#
#     all_values <- lapply(1:nsim,function(i){
#       dist_mat <- dist_matrices[[i]]
#       k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
#       g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
#       if(verbose){
#         setTxtProgressBar(pb, i)
#       }
#       return(cbind(k_vals,g_vals))
#     })
#
#   }else{
#     # the case where we can not simplified the situation
#     all_values <- lapply(1:nsim,function(i){
#       dist_mat <- randomize_distmatrix(graph_result$graph,graph_result$spedges,n)
#       k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
#       g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
#       if(verbose){
#         setTxtProgressBar(pb, i)
#       }
#       return(cbind(k_vals,g_vals))
#     })
#   }
#
#   ## step8 : extract the k_vals and g_vals matrices
#   k_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,1])}))
#   g_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,2])}))
#
#   ## step9 : calculating the summary stats
#   upper <- 1-conf_int / 2
#   lower <- conf_int / 2
#   k_stats <- apply(k_mat,MARGIN = 1, function(i){
#     return(quantile(i,probs = c(lower,upper)))
#   })
#   g_stats <- apply(g_mat,MARGIN = 1, function(i){
#     return(quantile(i,probs = c(lower,upper)))
#   })
#
#   plot_df <- data.frame(
#     "obs_k" = k_vals,
#     "lower_k" = k_stats[1,],
#     "upper_k" = k_stats[2,],
#     "obs_g" = g_vals,
#     "lower_g" = g_stats[1,],
#     "upper_g" = g_stats[2,],
#     "distances" = seq(start,end,step)
#   )
#
#   plotk <- ggplot(plot_df)+
#     geom_ribbon(aes_string(x = "distances", ymin="lower_k",ymax = "upper_k"),
#                 fill = grDevices::rgb(0.1,0.1,0.1), alpha=0.4)+
#     geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
#     labs(x = "distances",
#          y = "empirical K-function")
#
#   plotg <- ggplot(plot_df)+
#     geom_ribbon(aes_string(x = "distances", ymin="lower_g", ymax = "upper_g"),
#                 fill = grDevices::rgb(0.1,0.1,0.1),alpha=0.4, )+
#     geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
#     labs(x = "distances",
#          y = "empirical G-function")
#
#   obj <- list(
#     "plotk" = plotk,
#     "plotg" = plotg,
#     "values" = plot_df
#   )
#   if(return_sims){
#     obj$sim_k_values <- k_mat
#     obj$sim_g_values <- g_mat
#   }
#
#   return(obj)
# }



# kfunctions.mc <- function(lines, points, start, end, step, width, nsim, conf_int = 0.05,
#                           digits = 2 ,tol = 0.1, resolution = 50, agg = NULL,
#                           verbose = TRUE, return_sims = FALSE){
#
#   ## step0 : clean the points
#   if(verbose){
#     print("Preparing data ...")
#   }
#   n <- nrow(points)
#   points$goid <- seq_len(nrow(points))
#   points$weight <- rep(1,nrow(points))
#   points <- clean_events(points,digits,agg)
#
#   probs <- NULL
#
#   ## step1 : clean the lines
#   if(is.null(probs)){
#     lines$probs <- 1
#   }else{
#     lines$probs <- probs
#   }
#   lines$length <- as.numeric(st_length(lines))
#   lines <- subset(lines, lines$length>0)
#   lines$oid <- seq_len(nrow(lines))
#
#   ## step2 : adding the points to the lines
#   if(verbose){
#     print("Snapping points on lines ...")
#   }
#   snapped_events <- snapPointsToLines2(points,lines,idField = "oid")
#   new_lines <- split_lines_at_vertex(lines, snapped_events,
#                                      snapped_events$nearest_line_id, tol)
#   #new_lines <- add_vertices_lines(lines,snapped_events,snapped_events$nearest_line_id,tol)
#
#   ## step3 : splitting the lines
#   if(verbose){
#     print("Building graph ...")
#   }
#   #new_lines <- simple_lines(new_lines)
#   new_lines$length <- as.numeric(st_length(new_lines))
#   new_lines <- subset(new_lines,new_lines$length>0)
#   new_lines <- remove_loop_lines(new_lines,digits)
#   new_lines$oid <- seq_len(nrow(new_lines))
#   new_lines <- new_lines[c("length","oid","probs")]
#   Lt <- sum(as.numeric(st_length(new_lines)))
#
#   ## step4 : building the graph for the real case
#   graph_result <- build_graph(new_lines,digits = digits,line_weight = "length", attrs = TRUE)
#   graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
#                                                            name = "probs")
#   graph <- graph_result$graph
#   nodes <- graph_result$spvertices
#   snapped_events$vertex_id <- closest_points(snapped_events, nodes)
#
#   ## HERE, I SHOULD CHECK THAT TWO POINTS DO NOT SHARE THE SAME vertex
#   if(max(table(snapped_events$vertex_id))>1){
#     warning("After snapping the points on the network, some of them share the same location.
#          To correct it, please consider setting or increasing the value of the parameter agg.
#          They will be merged and their weights added")
#   }
#
#   ## step5 : calculating the distance matrix
#   # dist_mat <- igraph::distances(graph,v = snapped_events$vertex_id,
#   #                               to = snapped_events$vertex_id)
#   dist_mat <- dist_mat_dupl(graph, snapped_events$vertex_id, snapped_events$vertex_id)
#
#   ## step6 : calcualte the kfunction and the g function
#   if(verbose){
#     print("Calculating k and g functions ...")
#   }
#   k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = snapped_events$weight)
#   g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = snapped_events$weight)
#
#   ## step7 : generate the permutations
#   if(verbose){
#     print("Calculating the simulations ...")
#   }
#
#   w <- rep(1,times = n)
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
#           dist_mat <- randomize_distmatrix(graph, edgesdf, n)
#           k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
#           g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
#           return(cbind(k_vals,g_vals))
#         },future.packages = c("igraph"))
#       })
#     }else{
#       all_values <- future.apply::future_lapply(sim_seq, function(i){
#         dist_mat <- randomize_distmatrix(graph, edgesdf, n)
#         k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
#         g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
#         return(cbind(k_vals,g_vals))
#       },future.packages = c("igraph"))
#
#     }
#   }else{
#     # the simplified way
#     ## first : generating the matrices
#     if(verbose){
#       print("generating the randomized distance matrices...")
#     }
#     dist_mats <- randomize_distmatrix2(graph, edgesdf,
#                                       n = n, nsim = nsim,
#                                       resolution = resolution)
#     if(verbose){
#       print("calculating the k and g functions for the randomized matrices..")
#     }
#     all_values <- future.apply::future_lapply(dist_mats, function(dist_mat){
#       k_vals <- kfunc_cpp(dist_mat,start,end,step,Lt,n,w = w)
#       g_vals <- gfunc_cpp(dist_mat,start,end,step,width,Lt,n,w = w)
#       return(cbind(k_vals,g_vals))
#     },future.packages = c("igraph"))
#
#   }
#
#
#   ## step8 : extract the k_vals and g_vals matrices
#   k_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,1])}))
#   g_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,2])}))
#
#   ## step9 : calculating the summary stats
#   lower <- 1-conf_int/2
#   upper <- conf_int/2
#   k_stats <- apply(k_mat,MARGIN = 1, function(i){
#     return(quantile(i,probs = c(lower,upper)))
#   })
#   g_stats <- apply(g_mat,MARGIN = 1, function(i){
#     return(quantile(i,probs = c(lower,upper)))
#   })
#
#   plot_df <- data.frame(
#     "obs_k" = k_vals,
#     "upper_k" = k_stats[1,],
#     "lower_k" = k_stats[2,],
#     "obs_g" = g_vals,
#     "upper_g" = g_stats[1,],
#     "lower_g" = g_stats[2,],
#     "distances" = seq(start,end,step)
#   )
#
#   plotk <- ggplot(plot_df)+
#     geom_ribbon(aes_string(x = "distances", ymin="lower_k", ymax = "upper_k"),
#                 fill = grDevices::rgb(0.1,0.1,0.1),alpha=0.4, )+
#     geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
#     labs(x = "distances",
#          y = "empirical K-function")
#
#   plotg <- ggplot(plot_df)+
#     geom_ribbon(aes_string(x = "distances", ymin="lower_g", ymax = "upper_g"),
#                 fill = grDevices::rgb(0.1,0.1,0.1),alpha=0.4, )+
#     geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
#     labs(x = "distances",
#          y = "empirical G-function")
#
#   obj <- list(
#     "plotk" = plotk,
#     "plotg" = plotg,
#     "values" = plot_df
#   )
#
#   if(return_sims){
#     obj$sim_k_values <- k_mat
#     obj$sim_g_values <- g_mat
#   }
#
#   return(obj)
# }


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### execution cross-k functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# cross_kfunctions <- function(lines, pointsA, pointsB,
#                              start, end, step, width,
#                              nsim, conf_int = 0.05,
#                              digits = 2, tol = 0.1,
#                              resolution = NULL, agg = NULL,
#                              verbose = TRUE, return_sims = FALSE){
#
#   ## step0 : clean the points
#   if(verbose){
#     print("Preparing data ...")
#   }
#   na <- nrow(pointsA)
#   nb <- nrow(pointsB)
#
#   probs <- NULL
#
#   st_geometry(pointsA) <- 'geometry'
#   st_geometry(pointsB) <- 'geometry'
#
#   pointsA$weight <- rep(1,nrow(pointsA))
#   pointsA <- clean_events(pointsA,digits,agg)
#   pointsA$goid <- seq_len(nrow(pointsA))
#
#   pointsB$weight <- rep(1,nrow(pointsB))
#   pointsB <- clean_events(pointsB,digits,agg)
#   pointsB$goid <- seq_len(nrow(pointsB))
#
#   ## step1 : clean the lines
#   if(is.null(probs)){
#     lines$probs <- 1
#   }else{
#     lines$probs <- probs
#   }
#   lines$length <- as.numeric(st_length(lines))
#   lines <- subset(lines, lines$length>0)
#   lines$oid <- seq_len(nrow(lines))
#
#   ## step2 : adding the points to the lines
#   if(verbose){
#     print("Snapping points on lines ...")
#   }
#   pointsA$type <- "A"
#   pointsB$type <- "B"
#   all_events <- rbind(pointsA[c("type","goid","weight")],
#                       pointsB[c("type","goid","weight")])
#
#   snapped_events <- snapPointsToLines2(all_events, lines, idField = "oid")
#
#
#   new_lines <- split_lines_at_vertex(lines, snapped_events,
#                                      snapped_events$nearest_line_id, tol)
#
#   # HERE, I should ensure that snapped points have unique locations
#   # otherwise, I must add them
#
#
#
#   # new_lines <- add_vertices_lines(lines, snapped_events,
#   #                                 snapped_events$nearest_line_id, tol)
#
#   ## step3 : splitting the lines
#   if(verbose){
#     print("Building graph ...")
#   }
#   #new_lines <- simple_lines(new_lines)
#   new_lines$length <- as.numeric(st_length(new_lines))
#   new_lines <- subset(new_lines, new_lines$length>0)
#   new_lines <- remove_loop_lines(new_lines, digits)
#   new_lines$oid <- seq_len(nrow(new_lines))
#   new_lines <- new_lines[c("length","oid","probs")]
#   Lt <- sum(as.numeric(st_length(new_lines)))
#
#   ## step4 : building the graph for the real case
#   graph_result <- build_graph(new_lines,digits = digits,
#                               line_weight = "length", attrs = TRUE)
#   graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
#                                                            name = "probs")
#   graph <- graph_result$graph
#   nodes <- graph_result$spvertices
#   snapped_events$vertex_id <- closest_points(snapped_events, nodes)
#
#   ## HERE, I SHOULD CHECK THAT TWO POINTS DO NOT SHARE THE SAME vertex
#   if(max(table(snapped_events$vertex_id))>1){
#     warning("After snapping the points on the network, some of them share the same location.
#          To correct it, please consider setting or increasing the value of the parameter agg.
#          They will be merged and their weights added")
#   }
#
#   snappedA <- subset(snapped_events, snapped_events$type == "A")
#   snappedB <- subset(snapped_events, snapped_events$type == "B")
#
#   ## step5 : calculating the distance matrix
#
#   dist_mat <- dist_mat_dupl(graph, snappedB$vertex_id, snappedA$vertex_id, mode = "out")
#
#
#   ## step6 : calcualte the kfunction and the g function
#   if(verbose){
#     print("Calculating k and g functions ...")
#   }
#   k_vals <- cross_kfunc_cpp(dist_mat,start,end,step,Lt,na,nb,snappedA$weight,snappedB$weight)
#   g_vals <- cross_gfunc_cpp(dist_mat,start,end,step,width,Lt,na,nb,snappedA$weight,snappedB$weight)
#
#   ## step7 : generate the permutations
#   if(verbose){
#     print("Calculating the simulations ...")
#     pb <- txtProgressBar(min = 0, max = nsim, style = 3)
#   }
#   w <- rep(1,times = na)
#
#   # the case where we can simplified the situation
#   if (is.null(resolution)==FALSE){
#     dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
#                                            edge_df = graph_result$spedges,
#                                            n = na,
#                                            start_vert = snappedB$vertex_id,
#                                            resolution = resolution,
#                                            nsim = nsim)
#
#     all_values <- lapply(1:nsim,function(i){
#       dist_mat <- dist_matrices[[i]]
#       k_vals <- cross_kfunc_cpp(dist_mat,start,end,step,Lt,na,nb,w,snappedB$weight)
#       g_vals <- cross_gfunc_cpp(dist_mat,start,end,step,width,Lt,na,nb,w,snappedB$weight)
#       if(verbose){
#         setTxtProgressBar(pb, i)
#       }
#       return(cbind(k_vals,g_vals))
#     })
#
#   }else{
#     # the case where we can not simplified the situation
#     all_values <- lapply(1:nsim,function(i){
#       dist_mat <- randomize_distmatrix(graph_result$graph,graph_result$spedges,
#                                        na,start_vert = snappedB$vertex_id)
#       k_vals <- cross_kfunc_cpp(dist_mat,start,end,step,Lt,na,nb,w,snappedB$weight)
#       g_vals <- cross_gfunc_cpp(dist_mat,start,end,step,width,Lt,na,nb,w,snappedB$weight)
#       if(verbose){
#         setTxtProgressBar(pb, i)
#       }
#       return(cbind(k_vals,g_vals))
#     })
#   }
#
#
#   ## step8 : extract the k_vals and g_vals matrices
#   k_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,1])}))
#   g_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,2])}))
#
#   ## step9 : calculating the summary stats
#   upper <- 1-conf_int / 2
#   lower <- conf_int / 2
#   k_stats <- apply(k_mat,MARGIN = 1, function(i){
#     return(quantile(i,probs = c(lower,upper)))
#   })
#   g_stats <- apply(g_mat,MARGIN = 1, function(i){
#     return(quantile(i,probs = c(lower,upper)))
#   })
#
#   plot_df <- data.frame(
#     "obs_k" = k_vals,
#     "lower_k" = k_stats[1,],
#     "upper_k" = k_stats[2,],
#     "obs_g" = g_vals,
#     "lower_g" = g_stats[1,],
#     "upper_g" = g_stats[2,],
#     "distances" = seq(start,end,step)
#   )
#
#   plotk <- ggplot(plot_df)+
#     geom_ribbon(aes_string(x = "distances", ymin = "lower_k", ymax = "upper_k"),
#                 fill = rgb(0.1,0.1,0.1),alpha=0.4, )+
#     geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
#     labs(x = "distances",
#          y = "empirical cross-K-function")
#
#   plotg <- ggplot(plot_df)+
#     geom_ribbon(aes_string(x = "distances", ymin = "lower_g", ymax = "upper_g"),
#                 fill = rgb(0.1,0.1,0.1),alpha=0.4, )+
#     geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
#     labs(x = "distances",
#          y = "empirical cross-G-function")
#
#   obj <- list(
#     "plotk" = plotk,
#     "plotg" = plotg,
#     "values" = plot_df
#   )
#   if(return_sims){
#     obj$sim_k_values <- k_mat
#     obj$sim_g_values <- g_mat
#   }
#   return(obj)
# }



# cross_kfunctions.mc <- function(lines, pointsA, pointsB,
#                                 start, end, step, width,
#                                 nsim, conf_int = 0.05,
#                                 digits = 2, tol = 0.1,
#                                 resolution = NULL, agg = NULL,
#                                 verbose = TRUE, return_sims = FALSE){
#
#   ## step0 : clean the points
#   if(verbose){
#     print("Preparing data ...")
#   }
#   na <- nrow(pointsA)
#   nb <- nrow(pointsB)
#
#   probs <- NULL
#
#   pointsA$weight <- rep(1,nrow(pointsA))
#   pointsA <- clean_events(pointsA,digits,agg)
#   pointsA$goid <- seq_len(nrow(pointsA))
#
#   pointsB$weight <- rep(1,nrow(pointsB))
#   pointsB <- clean_events(pointsB,digits,agg)
#   pointsB$goid <- seq_len(nrow(pointsB))
#
#   ## step1 : clean the lines
#   if(is.null(probs)){
#     lines$probs <- 1
#   }else{
#     lines$probs <- probs
#   }
#   lines$length <- as.numeric(st_length(lines))
#   lines <- subset(lines, lines$length>0)
#   lines$oid <- seq_len(nrow(lines))
#
#   ## step2 : adding the points to the lines
#   if(verbose){
#     print("Snapping points on lines ...")
#   }
#   pointsA$type <- "A"
#   pointsB$type <- "B"
#   all_events <- rbind(pointsA[c("type","goid","weight")],
#                       pointsB[c("type","goid","weight")])
#
#   snapped_events <- snapPointsToLines2(all_events, lines, idField = "oid")
#   new_lines <- split_lines_at_vertex(lines, snapped_events,
#                                   snapped_events$nearest_line_id, tol)
#   # new_lines <- add_vertices_lines(lines, snapped_events,
#   #                                 snapped_events$nearest_line_id, tol)
#
#   ## step3 : splitting the lines
#   if(verbose){
#     print("Building graph ...")
#   }
#   #new_lines <- simple_lines(new_lines)
#   new_lines$length <- as.numeric(st_length(new_lines))
#   new_lines <- subset(new_lines,new_lines$length>0)
#   new_lines <- remove_loop_lines(new_lines,digits)
#   new_lines$oid <- seq_len(nrow(new_lines))
#   new_lines <- new_lines[c("length","oid","probs")]
#   Lt <- sum(as.numeric(st_length(new_lines)))
#
#   ## step4 : building the graph for the real case
#   graph_result <- build_graph(new_lines,digits = digits,
#                               line_weight = "length", attrs = TRUE)
#   graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
#                                                            name = "probs")
#   graph <- graph_result$graph
#   nodes <- graph_result$spvertices
#   snapped_events$vertex_id <- closest_points(snapped_events, nodes)
#
#   ## HERE, I SHOULD CHECK THAT TWO POINTS DO NOT SHARE THE SAME vertex
#   if(max(table(snapped_events$vertex_id))>1){
#     stop("After snapping the points on the network, some of them share the same location.
#          To correct it, please consider setting or increasing the value of the parameter agg.
#          They will be merged and their weights added")
#   }
#
#   snappedA <- subset(snapped_events, snapped_events$type == "A")
#   snappedB <- subset(snapped_events, snapped_events$type == "B")
#
#   ## step5 : calculating the distance matrix
#   dist_mat <- igraph::distances(graph,v = snappedB$vertex_id,
#                                 to = snappedA$vertex_id)
#   ## step6 : calcualte the kfunction and the g function
#   if(verbose){
#     print("Calculating k and g functions ...")
#   }
#   k_vals <- cross_kfunc_cpp(dist_mat, start, end, step, Lt, na, nb,
#                         snappedA$weight, snappedB$weight)
#   g_vals <- cross_gfunc_cpp(dist_mat, start, end, step, width, Lt, na,
#                        nb, snappedA$weight, snappedB$weight)
#
#   ## step7 : generate the permutations
#   if(verbose){
#     print("Calculating the simulations ...")
#   }
#   w <- rep(1,times = na)
#   sim_seq <- 1:nsim
#
#   # classical approach
#   if (is.null(resolution)){
#     progressr::with_progress({
#       p <- progressr::progressor(along = sim_seq)
#       all_values <- future.apply::future_lapply(sim_seq,function(i){
#         dist_mat <- randomize_distmatrix(graph_result$graph,
#                                          graph_result$spedges,
#                                          na, start_vert = snappedB$vertex_id)
#         k_vals <- cross_kfunc_cpp(dist_mat, start, end, step, Lt, na, nb,
#                               w, snappedB$weight)
#         g_vals <- cross_gfunc_cpp(dist_mat, start, end, step, width, Lt, na, nb,
#                              w,snappedB$weight)
#         return(cbind(k_vals,g_vals))
#       },future.packages = c("igraph","base"))
#     })
#   }else{
#     #simplified approach
#     if(verbose){
#       print("calculating the randomized distance matrices...")
#     }
#     dist_matrices <- randomize_distmatrix2(graph = graph_result$graph,
#                                            edge_df = graph_result$spedges,
#                                            n = na,
#                                            start_vert = snappedB$vertex_id,
#                                            resolution = resolution,
#                                            nsim = nsim)
#     if(verbose){
#       print("calculating the k and g functions for the randomized matrices..")
#     }
#     all_values <- future.apply::future_lapply(dist_matrices,function(dist_mat){
#         k_vals <- cross_kfunc_cpp(dist_mat, start, end, step, Lt, na, nb,
#                               w, snappedB$weight)
#         g_vals <- cross_gfunc_cpp(dist_mat, start, end, step, width, Lt, na, nb,
#                              w,snappedB$weight)
#         return(cbind(k_vals,g_vals))
#       },future.packages = c("igraph","base"))
#   }
#
#
#   ## step8 : extract the k_vals and g_vals matrices
#   k_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,1])}))
#   g_mat <- do.call(cbind,lapply(all_values,function(i){return(i[,2])}))
#
#   ## step9 : calculating the summary stats
#   upper <- 1-conf_int / 2
#   lower <- conf_int / 2
#   k_stats <- apply(k_mat,MARGIN = 1, function(i){
#     return(quantile(i,probs = c(lower,upper)))
#   })
#   g_stats <- apply(g_mat,MARGIN = 1, function(i){
#     return(quantile(i,probs = c(lower,upper)))
#   })
#
#   plot_df <- data.frame(
#     "obs_k" = k_vals,
#     "lower_k" = k_stats[1,],
#     "upper_k" = k_stats[2,],
#     "obs_g" = g_vals,
#     "lower_g" = g_stats[1,],
#     "upper_g" = g_stats[2,],
#     "distances" = seq(start,end,step)
#   )
#
#   plotk <- ggplot(plot_df)+
#     geom_ribbon(aes_string(x = "distances", ymin = "lower_k", ymax = "upper_k"),
#                 fill = rgb(0.1,0.1,0.1),alpha=0.4, )+
#     geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
#     labs(x = "distances",
#          y = "empirical cross-K-function")
#
#   plotg <- ggplot(plot_df)+
#     geom_ribbon(aes_string(x = "distances", ymin = "lower_g", ymax = "upper_g"),
#                 fill = rgb(0.1,0.1,0.1),alpha=0.4, )+
#     geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
#               linetype="dashed")+
#     geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
#     labs(x = "distances",
#          y = "empirical cross-G-function")
#
#   obj <- list(
#     "plotk" = plotk,
#     "plotg" = plotg,
#     "values" = plot_df
#   )
#   if(return_sims){
#     obj$sim_k_values <- k_mat
#     obj$sim_g_values <- g_mat
#   }
#   return(obj)
# }
#
