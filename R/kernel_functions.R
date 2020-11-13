# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### available kernels ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title triangle kernel
#'
#' @description Function implementing the triangle kernel.
#'
#' @param d The distance from the event
#' @param bw The bandwidth used for the kernel
#' @return The estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
triangle_kernel <- function(d, bw){
  u <- d/bw
  k <- 1 - abs(u)
  k <- k/bw
  k <- ifelse(abs(d)>bw,0,k)
  return(k)
}

#' @title Uniform kernel
#'
#' @description Function implementing the uniform kernel.
#'
#' @param d The distance from the event
#' @param bw The bandwidth used for the kernel
#' @return The estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
uniform_kernel <- function(d, bw){
  k <- 1/(2*bw)
  k <- ifelse(abs(d)>bw,0,k)
  return(k)
}

#' @title Epanechnikov kernel
#'
#' @description Function implementing the epanechnikov kernel.
#'
#' @param d The distance from the event
#' @param bw The bandwidth used for the kernel
#' @return The estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
epanechnikov_kernel <- function(d, bw){
  u <- d/bw
  k <- (3/4) * (1-u**2)
  k <- k/bw
  k <- ifelse(abs(d)>bw,0,k)
  return(k)
}

#' @title Quartic kernel
#'
#' @description Function implementing the quartic kernel.
#'
#' @param d The distance from the event
#' @param bw The bandwidth used for the kernel
#' @return The estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
quartic_kernel <- function(d, bw){
  u <- d/bw
  k <- (15/16)*(1-u**2)**2
  k <- k / bw
  k <- ifelse(abs(d)>bw,0,k)
  return(k)
}


#' @title Triweight kernel
#'
#' @description Function implementing the triweight kernel.
#'
#' @param d The distance from the event
#' @param bw The bandwidth used for the kernel
#' @return The estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
triweight_kernel <- function(d, bw){
  u <- d/bw
  k <- (35/32)*(1-u**2)**3
  k <- k / bw
  k <- ifelse(abs(d)>bw,0,k)
  return(k)
}

#' @title Tricube kernel
#'
#' @description Function implementing the tricube kernel.
#'
#' @param d The distance from the event
#' @param bw The bandwidth used for the kernel
#' @return The estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
tricube_kernel <- function(d, bw){
  u <- d/bw
  k <- (70/81)*(1-abs(u)**3)**3
  k <- k / bw
  k <- ifelse(abs(d)>bw,0,k)
  return(k)
}

#' @title Cosine kernel
#'
#' @description Function implementing the cosine kernel.
#'
#' @param d The distance from the event
#' @param bw The bandwidth used for the kernel
#' @return The estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
cosine_kernel <- function(d, bw){
  u <- abs(d/bw)
  k <- (pi/4) * cos((pi/2)*u)
  k <- k / bw
  k <- ifelse(abs(d)>bw,0,k)
  return(k)
}

#' @title Gaussian kernel
#'
#' @description Function implementing the gaussian kernel.
#'
#' @param d The distance from the event
#' @param bw The bandwidth used for the kernel
#' @return The estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
gaussian_kernel <- function(d, bw){
  u <- d/bw
  t1 <- 1/(sqrt(2*pi))
  t2 <- exp(-1 * (1/2) * u**2  )
  k <- t1*t2
  k <- k/bw
  k <- ifelse(abs(d)>bw,0,k)
  return(k)
}


#' @title Scaled gaussian kernel
#'
#' @description Function implementing the scaled gaussian kernel.
#'
#' @param d The distance from the event
#' @param bw The bandwidth used for the kernel
#' @return The estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
gaussian_kernel_scaled <- function(d, bw){
  bw2 <- bw/3
  t1 <- 1/(bw2 * sqrt(2*pi))
  t2 <- exp(-1 * (d**2 / (2*bw2**2)))
  k <- t1*t2
  k <- ifelse(abs(d)>bw,0,k)
  return(k)
}



#' @title Select kernel function
#'
#' @description select the kernel function with its name.
#'
#' @param name The name of the kernel to use
#' @return A kernel function
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
select_kernel <- function(name){
  if((name %in% c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov','uniform'))==FALSE){
    stop('the name of the kernel function must be one of c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov" ,"uniform")')
  }
  if(name=="gaussian"){
    return(gaussian_kernel)
  }
  if(name == "scaled gaussian"){
    return(gaussian_kernel_scaled)
  }
  if(name=="triangle"){
    return(triangle_kernel)
  }
  if(name=="tricube"){
    return(tricube_kernel)
  }
  if(name=="triweight"){
    return(triweight_kernel)
  }
  if(name=="quartic"){
    return(quartic_kernel)
  }
  if(name=="cosine"){
    return(cosine_kernel)
  }
  if(name=="epanechnikov"){
    return(epanechnikov_kernel)
  }
  if(name=="uniform"){
    return(uniform_kernel)
  }

}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Functions for adaptative bandwidth calculation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Geometric mean
#'
#' @description Function to calculate the geometric mean.
#'
#' @param x A vector of numeric values
#' @param na.rm A boolean indicating if we filter the NA values
#' @return The geometric mean of x
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#' @title Gamma parameter for Abramson’s adaptive bandwidth
#'
#' @description Function to calculate the gamma parameter in Abramson’s smoothing regimen.
#'
#' @param k a vector of numeric values (the estimated kernel densities)
#' @return the gamma parameter in Abramson’s smoothing regimen
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
calc_gamma <- function(k){
  return(gm_mean(k**(-1/2)))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Functions to perform the simple NKDE ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Simple NKDE algorithm
#'
#' @description Function to perform the simple nkde.
#'
#' @param graph a graph object from igraph representing the network
#' @param events a SpatialPointsDataFrame representing the events. It must be
#' snapped on the network, and be nodes of the network. A column vertex_id
#' must indicate for each event its corresponding node
#' @param samples a SpatialPointsDataFrame representing the sampling points.
#' The samples must be snapped on the network. A column edge_id must indicate
#' for each sample on which edge it is snapped.
#' @param bws a vector indicating the kernel bandwidth (in meters) for each
#' event
#' @param kernel_func a function obtained with the function select_kernel
#' @param nodes a SpatialPointsDataFrame representing the nodes of the network
#' @param edges a SpatialLinesDataFrame representing the edges of the network
#' @return a dataframe with two columns. sum_k is the sum for each sample point
#'  of the kernel values. n is the number of events influencing each sample
#' point
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
simple_nkde <- function(graph, events, samples, bws, kernel_func, nodes, edges){

  ## step 1 set all values to 0
  base_k <- rep(0,nrow(samples))
  base_count <- rep(0,nrow(samples))

  #if the number of event is 0, then return only 0 values
  if(nrow(events)==0){
    return(data.frame("sum_k"=base_k,
                      "n"=base_count))
  }

  sample_tree <- build_quadtree(samples)
  edges_tree <- build_quadtree(edges)
  ## step2 iterate over each event
  pb <- txtProgressBar(min = 0, max = nrow(events), style = 3)
  for(i in 1:nrow(events)){
    setTxtProgressBar(pb, i)
    bw <- bws[[i]]
    #extracting the starting values
    e <- events[i,]
    y <- e$vertex_id
    w <- e$weight
    #calculating the kernel values
    samples_k <- ess_kernel(graph,y,bw, kernel_func, samples, nodes, edges, sample_tree, edges_tree)
    samples_count <- ifelse(samples_k>0,1,0)
    base_k <- samples_k * w + base_k
    base_count <- base_count + (samples_count * w)
  }
  return(data.frame("sum_k"=base_k,
              "n"=base_count))
}


#' @title Worker for simple NKDE algorithm
#'
#' @description The worker function to perform the simple nkde.
#'
#' @param graph a graph object from igraph representing the network
#' @param y the index of the actual event
#' @param bw a float indicating the kernel bandwidth (in meters)
#' @param kernel_func a function obtained with the function select_kernel
#' @param samples a SpatialPointsDataFrame representing the sampling points.
#' The samples must be snapped on the network. A column edge_id must indicate
#' for each sample on which edge it is snapped.
#' @param nodes a SpatialPointsDataFrame representing the nodes of the network
#' @param edges a SpatialLinesDataFrame representing the edges of the network
#' @param sample_tree a quadtree object, the spatial index of the samples
#' @param edges_tree a quadtree object, the spatial index of the edges
#' @importFrom igraph get.edge.attribute
#' @importFrom rgeos gBuffer
#' @importFrom igraph ends distances
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
ess_kernel <- function(graph, y, bw, kernel_func, samples, nodes, edges, sample_tree, edges_tree){
  samples_k <- rep(0,nrow(samples))
  event_node <- nodes[y,]
  buff <- gBuffer(event_node,width = bw)
  ## step1 find all the samples in the radius
  ok_samples <- spatial_request(buff,sample_tree,samples)
  if(nrow(ok_samples)==0){
    return(samples_k)
  }
  ## step2 find all the edges in the radius
  buff2 <- gBuffer(event_node,width = bw+0.1*bw)
  ok_edges <- spatial_request(buff2,edges_tree,edges)$edge_id
  ## Step3 for each edge, find its two vertices
  vertices <- ends(graph,ok_edges,names=F)
  ## step4 calculate the the distance between the start node and each edge vertex
  un_vertices <- unique(c(vertices[,1],vertices[,2]))
  dist1 <- as.numeric(distances(graph,y,to=un_vertices,mode="out"))

  dist_table <- data.frame("vertex"=un_vertices,
                           "distance" = dist1)
  ## step5 aggregate all the data
  df_edges <- data.frame("edge_id" = ok_edges,
                         "node1" = vertices[,1],
                         "node2" = vertices[,2]
                         )
  A1 <- data.table(df_edges)
  A2 <- data.table(df_edges)
  B <- data.table(dist_table)
  df_edges$d1 <- A1[B, on = c("node1" = "vertex"),
                   names(B) := mget(paste0("i.", names(B)))]$distance
  df_edges$d2 <- A2[B, on = c("node2" = "vertex"),
                   names(B) := mget(paste0("i.", names(B)))]$distance
  #df_edges$d1 <- left_join(df_edges,dist_table,by=c("node1"="vertex"))$distance
  #df_edges$d2 <- left_join(df_edges,dist_table,by=c("node2"="vertex"))$distance

  ## now, we will calculate for each sample, the minimum distance for the two distances
  df1 <- df_edges[c("edge_id","node1","d1")]
  df2 <- df_edges[c("edge_id","node2","d2")]

  ## getting the first part of the distances
  df_samples1 <- data.table(ok_samples@data)
  df_samples2 <- data.table(ok_samples@data)
  B <- data.table(df1)
  C <- data.table(df2)

  df_samples1[B, on = "edge_id", names(B) := mget(paste0("i.", names(B)))]
  df_samples2[C, on = "edge_id", names(C) := mget(paste0("i.", names(C)))]
  #df_samples1 <-  left_join(ok_samples@data,df1,by="edge_id")
  #df_samples2 <-  left_join(ok_samples@data,df2,by="edge_id")

  ## getting the second part of the distance
  start_nodes1 <- nodes@data[df_samples1$node1,]
  start_nodes2 <- nodes@data[df_samples2$node2,]

  dist1_b <- sqrt((df_samples1$X_coords - start_nodes1$X_coords)**2 +
                    (df_samples1$Y_coords - start_nodes1$Y_coords)**2)

  dist2_b <- sqrt((df_samples2$X_coords - start_nodes2$X_coords)**2 +
                    (df_samples2$Y_coords - start_nodes2$Y_coords)**2)

  ## getting the full distances
  full_dist1 <- dist1_b + df_samples1$d1
  full_dist2 <- dist2_b + df_samples2$d2

  ok_distances <-ifelse(full_dist1 < full_dist2, full_dist1, full_dist2)


  ##and then, the final distance

  k <- kernel_func(ok_distances,bw)
  samples_k[df_samples1$oid] <- k
  return(samples_k)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Functions to perform the discontinuous NKDE ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Function to perform the discontinuous nkde
# this function is deprecated and repalced by a Rcpp one
# the code is keeped for debugging purpose
# discontinuous_nkde2 <-  function(edge_list,neighbour_list, v_events, weights,
#                                  samples, bw, kernel_func, nodes, linelist, max_depth, verbose){
#
#   ##step 1 : mettre toutes les valeurs a 0
#   base_k <- rep(0,nrow(samples))
#   base_count <- rep(0,nrow(samples))
#
#
#   if(length(v_events)==0){
#     return(data.frame("sum_k"=base_k,
#                       "n"=base_count))
#   }
#   lines_weight <- linelist$weight
#
#   pb <- txtProgressBar(min = 0, max = length(v_events), style = 3)
#   ##step2 : iterer sur chaque event
#   for(i in 1:length(v_events)){
#     #preparer les differentes valeurs de departs pour l'event y
#     setTxtProgressBar(pb, i)
#     y <- v_events[[i]]
#     w <- weights[[i]]
#
#     samples_k <- esd_kernel2(y, edge_list,neighbour_list,
#                             samples, bw, kernel_func, nodes, lines_weight, linelist, max_depth, verbose)
#
#     samples_count <- ifelse(samples_k>0,1,0)
#     base_k <- samples_k * w + base_k
#     base_count <- base_count + (samples_count * w)
#   }
#   return(data.frame("sum_k"=base_k,
#                     "n"=base_count))
# }


#' Worker function for the discontinuous nkde
# this function is deprecated and repalced by a Rcpp one
# the code is keeped for debugging purpose

# esd_kernel2 <- function(y,edge_list,neighbour_list,
#                         samples, bw, kernel_func, nodes, lines_weight, linelist, max_depth, verbose){
#
#   #definit les premiere valeurs a 0
#   samples_k <- rep(0,nrow(samples))
#   #definir la premiere serie de parametres
#   all_parameters <- list(
#     list("v"=y,
#          "d"=0,
#          "alpha" = 1,
#          "depth" = 0,
#          "prev_node"=-999)
#   )
#
#   #lancement des iterations
#   while(length(all_parameters)>0){
#
#     #on cree une petite liste vide
#     new_parameters <- list()
#
#     #on itere sur les cas en cours
#     for(params in all_parameters){
#       #step1 : unpacking the values
#       v <- params[["v"]]
#       alpha <- params[["alpha"]]
#       prev_node <- params[["prev_node"]]
#       d <- params[["d"]]
#       depth <- params[["depth"]]
#
#       #step2 : on trouve les voisins de v (et on enleve le precedent node)
#       v_neighbours <- neighbour_list[[v]]
#       test <- v_neighbours != prev_node
#       v_neighbours <- v_neighbours[test]
#
#       #avec ces voisins, on peut setter le new_alpha
#       if(prev_node!=-999){
#         new_alpha <- (1/(length(v_neighbours))) * alpha
#       }else{
#         new_alpha <- 1
#       }
#
#
#       if(length(v_neighbours)>0){
#
#         #step3 : on trouve les edges entre v et ses voisins
#         edges_id <- paste(v,v_neighbours,sep="_")
#         edges <- sapply(edges_id,function(x){return(edge_list[[x]])})
#         node <- nodes[v,]
#
#         #step4 : on va iterer sur chacune de ces lignes
#         for(i in 1:length(edges)){
#           li <- edges[[i]]
#           vi <- v_neighbours[[i]]
#           #il faut trouver les echantillons concernes
#           test <- samples$edge_id == li
#           sub_samples <- samples[test,]
#           #il faut calculer les distances entre le point de depart et cet echantillon
#           d1 <- sqrt((node$X_coords - sub_samples$X_coords)**2 + (node$Y_coords - sub_samples$Y_coords)**2)
#           d2 <- d1 + d
#           #on calcule maintenant la valeur kernel
#           k1 <- kernel_func(d2,bw) * new_alpha
#           #et on l'ajoute a la valeur precedente
#           old_k <- samples_k[test]
#           samples_k[test] <- old_k + k1
#
#           #il ne reste plus que a voir si on peut continuer sur le prochain noeud
#           d3 <- d+lines_weight[li]
#           if(d3<bw & depth<max_depth){
#
#             if(new_alpha == alpha){
#               new_depth <- depth
#             }else{
#               new_depth <- depth+1
#             }
#
#             new_params <- list("v" = vi,
#                                "prev_node"=v,
#                                "d" = d3,
#                                "depth" = new_depth,
#                                "alpha" = new_alpha
#                                )
#             new_parameters[[length(new_parameters)+1]] <- new_params
#           }
#
#         }
#
#       }
#
#     }
#
#     #on reset les nouveaux parameters
#     all_parameters <- new_parameters
#   }
#
#   return(samples_k)
# }


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Functions to perform the continuous NKDE ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Function to perform the continuous nkde
#
# This function is deprecated because replaced by a Rcpp one

# continuous_nkde<- function(graph,events,samples,bw,kernel_func,nodes, line_list, max_depth){
#   ##step 1 : mettre toutes les valeurs a 0
#   base_k <- rep(0,nrow(samples))
#   base_count <- rep(0,nrow(samples))
#   ##step2 : iterer sur chaque event
#   for(i in 1:nrow(events)){
#     #preparer les differentes valeurs de departs pour l'event y
#     e <- events[i,]
#     y <- e$vertex_id
#     w <- e$weight
#     #on veut trouver toutes les lignes emannant de y (Ly)
#     y_neighbours <- as.numeric(igraph::neighbors(graph,y,mode="out"))
#     eids <- cbind(rep(y,length(y_neighbours)),y_neighbours)
#     eids <- c(t(eids))
#     Ly <- as.numeric(igraph::get.edge.ids(graph,eids))
#     for(j in 1:length(Ly)){
#       li <- Ly[j]
#       vi <- y_neighbours[j]
#       d <- 0
#       depth <- 0
#       alpha <- 1
#       samples_k <<- rep(0,nrow(samples))
#       esc_kernel(graph,y,vi,li,d,alpha,bw,kernel_func, samples, nodes,line_list,depth,max_depth)
#       base_k <- samples_k * w + base_k
#       samples_count <- ifelse(samples_k>0,1,0)
#       base_count <- base_count + (samples_count * w)
#     }
#
#   }
#   return(data.frame("sum_k" = base_k,
#               "n" = base_count))
# }


# Worker function for the continuous nkde
#
# This function is deprecated because replaced by a Rcpp one.

# esc_kernel <- function(graph,v, v1, l1, d,alpha, bw,kernel_func, samples,nodes,line_list,depth,max_depth){
#   x_samples <- subset(samples, samples$edge_id==l1)
#   x_dists <- rgeos::gDistance(nodes[v,],x_samples,byid = T)[,1] + d
#   samples_k[x_samples$oid] <<- samples_k[x_samples$oid] + alpha * (kernel_func(x_dists, bw))
#   d2 <-  get.edge.attribute(graph,"weight",l1)
#   d2 <- d+d2
#   new_depth <- depth+1
#   if(bw>d2 & new_depth < max_depth){
#     #on veut trouver toutes les lignes emannant de v (Lv)
#     v_neighbours <- as.numeric(igraph::neighbors(graph,v1,mode="out"))
#     eids <- cbind(rep(v1,length(v_neighbours)),v_neighbours)
#     eids <- c(t(eids))
#     Lv <- as.numeric(igraph::get.edge.ids(graph,eids))
#     n <- length(v_neighbours)
#     if(n>1){
#       for(j in 1:length(Lv)){
#         li <- Lv[j]
#         vi <- v_neighbours[j]
#         if(li==l1){
#           n_alpha <- -1*alpha * ((n-2)/n)
#         }else{
#           n_alpha <- alpha * (2/n)
#         }
#         esc_kernel(graph,v1,vi, li, d2, n_alpha, bw,kernel_func, samples,nodes,line_list, new_depth, max_depth)
#       }
#     }
#   }
# }
