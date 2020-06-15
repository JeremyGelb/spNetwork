# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### available kernels ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' triangle kernel
#'
#' @param d the distance from the event
#' @param bw the bandwidth used for the kernel
#' @return the estimated density
#' @export
#' @examples
#' #This is an internal function, no example provided
triangle_kernel <- function(d, bw){
  u <- d/bw
  k <- 1 - abs(u)
  k <- k/bw
  k <- ifelse(d>bw,0,k)
  return(k)
}

#' epanechnikov kernel
#'
#' @param d the distance from the event
#' @param bw the bandwidth used for the kernel
#' @return the estimated density
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

#' quartic kernel
#'
#' @param d the distance from the event
#' @param bw the bandwidth used for the kernel
#' @return the estimated density
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


#' triweight kernel
#'
#' @param d the distance from the event
#' @param bw the bandwidth used for the kernel
#' @return the estimated density
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

#' tricube kernel
#'
#' @param d the distance from the event
#' @param bw the bandwidth used for the kernel
#' @return the estimated density
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

#' cosine kernel
#'
#' @param d the distance from the event
#' @param bw the bandwidth used for the kernel
#' @return the estimated density
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

#' gaussian kernel
#'
#' @param d the distance from the event
#' @param bw the bandwidth used for the kernel
#' @return the estimated density
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


#' select the right kernel function with its name
#'
#' @param name The name of the kernel to use
#' @return A kernel function
#' @examples
#' #This is an internal function, no example provided
select_kernel <- function(name){
  if((name %in% c("triangle", "gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov'))==FALSE){
    stop('the name of the kernel function must be one of c("triangle", "gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov")')
  }
  if(name=="gaussian"){
    return(gaussian_kernel)
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

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Functions to perform the simple NKDE ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Function to perform the simple nkde
#'
#' @param graph a graph object from igraph representing the network
#' @param events a SpatialPointsDataFrame representing the events. It must be
#' snapped on the network, and be nodes of the network. A column vertex_id
#' must indicate for each event its corresponding node
#' @param samples a SpatialPointsDataFrame representing the sampling points.
#' The samples must be snapped on the network. A column edge_id must indicate
#' for each sample on which edge it is snapped.
#' @param bw a float indicating the kernel bandwidth (in meters)
#' @param kernel_func a function obtained with the function select_kernel
#' @param nodes a SpatialPointsDataFrame representing the nodes of the network
#' @param edges a SpatialLinesDataFrame representing the edges of the network
#' @return a dataframe with two columns. sum_k is the sum for each sample point
#'  of the kernel values. n is the number of events influencing each sample
#' point
#' @examples
#' #This is an internal function, no example provided
simple_nkde <- function(graph, events, samples, bw, kernel_func, nodes, edges){
  ##step 1 : mettre toutes les valeurs a 0
  base_k <- rep(0,nrow(samples))
  base_count <- rep(0,nrow(samples))

  if(nrow(events)==0){
    return(data.frame("sum_k"=base_k,
                      "n"=base_count))
  }

  sample_tree <- build_quadtree(samples)
  edges_tree <- build_quadtree(edges)
  ##step2 : iterer sur chaque event
  pb <- txtProgressBar(min = 0, max = nrow(events), style = 3)
  for(i in 1:nrow(events)){
    #preparer les differentes valeurs de departs pour l'event y
    setTxtProgressBar(pb, i)
    e <- events[i,]
    y <- e$vertex_id
    w <- e$weight
    samples_k <- ess_kernel(graph,y,bw, kernel_func, samples, nodes, edges, sample_tree, edges_tree)
    samples_count <- ifelse(samples_k>0,1,0)
    base_k <- samples_k * w + base_k
    base_count <- base_count + (samples_count * w)
  }
  return(data.frame("sum_k"=base_k,
              "n"=base_count))
}


#' Worker function for the simple nkde
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
#' @importFrom dplyr left_join
#' @examples
#' #This is an internal function, no example provided
ess_kernel <- function(graph, y, bw, kernel_func, samples, nodes, edges, sample_tree, edges_tree){
  samples_k <- rep(0,nrow(samples))
  event_node <- nodes[y,]
  buff <- gBuffer(event_node,width = bw)
  ##step1 : find all the samples in the radius
  ok_samples <- spatial_request(buff,sample_tree,samples)
  if(nrow(ok_samples)==0){
    return(samples_k)
  }
  ##step2 : find all the edges
  ok_edges <- spatial_request(buff,edges_tree,edges)$edge_id
  ##Step3 : for each edge, find the two vertices
  vertices <- ends(graph,ok_edges,names=F)
  ##calculate the two distances (one for each vertex)
  un_vertices <- unique(c(vertices[,1],vertices[,2]))
  dist1 <- as.numeric(distances(graph,y,to=un_vertices,mode="out"))

  dist_table <- data.frame("vertex"=un_vertices,
                           "distance" = dist1)
  ##aggregate all the data
  df_edges <- data.frame("edge_id" = ok_edges,
                         "node1" = vertices[,1],
                         "node2" = vertices[,2]
                         )

  df_edges$d1 <- left_join(df_edges,dist_table,by=c("node1"="vertex"))$distance
  df_edges$d2 <- left_join(df_edges,dist_table,by=c("node2"="vertex"))$distance

  ## now, we will calculate for each sample, the minimum distance for the two distances
  df1 <- df_edges[c("edge_id","node1","d1")]
  df2 <- df_edges[c("edge_id","node2","d2")]

  ## getting the first part of the distances
  df_samples1 <-  left_join(ok_samples@data,df1,by="edge_id")
  df_samples2 <-  left_join(ok_samples@data,df2,by="edge_id")

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

#' Function to perform the discontinuous nkde
#'
#' @param graph a graph object from igraph representing the network
#' @param events a SpatialPointsDataFrame representing the events. It must be
#' snapped on the network, and be nodes of the network. A column vertex_id
#' must indicate for each event its corresponding node
#' @param samples a SpatialPointsDataFrame representing the sampling points.
#' The samples must be snapped on the network. A column edge_id must indicate
#' for each sample on which edge it is snapped.
#' @param bw a float indicating the kernel bandwidth (in meters)
#' @param kernel_func a function obtained with the function select_kernel
#' @param nodes a SpatialPointsDataFrame representing the nodes of the network
#' @return a dataframe with two columns. sum_k is the sum for each sample point
#'  of the kernel values. n is the number of events influencing each sample
#' point
#' @param edges a SpatialLinesDataFrame representing the edges of the network
#' @examples
#' #This is an internal function, no example provided
discontinuous_nkde <-  function(graph, events, samples, bw, kernel_func, nodes, edges){
  ##step 1 : mettre toutes les valeurs a 0
  base_k <- rep(0,nrow(samples))
  base_count <- rep(0,nrow(samples))

  if(nrow(events)==0){
    return(data.frame("sum_k"=base_k,
                      "n"=base_count))
  }

  sample_tree <- build_quadtree(samples)
  edges_tree <- build_quadtree(edges)

  ##step2 : iterer sur chaque event
  pb <- txtProgressBar(min = 0, max = nrow(events), style = 3)
  for(i in 1:nrow(events)){
    setTxtProgressBar(pb, i)
    #preparer les differentes valeurs de departs pour l'event y
    e <- events[i,]
    y <- e$vertex_id
    w <- e$weight
    samples_k <- esd_kernel(graph,y,bw, kernel_func, samples, nodes, edges, sample_tree, edges_tree)
    samples_count <- ifelse(samples_k>0,1,0)
    base_k <- samples_k * w + base_k
    base_count <- base_count + (samples_count * w)
  }
  return(data.frame("sum_k"=base_k,
              "n"=base_count))
}



#' Worker function for the discontinuous nkde
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
#' @importFrom igraph ends distances shortest_paths adjacent_vertices
#' @importFrom dplyr left_join
#' @examples
#' #This is an internal function, no example provided
esd_kernel <- function(graph, y, bw, kernel_func, samples, nodes, edges, sample_tree, edges_tree){
  samples_k <- rep(0,nrow(samples))
  event_node <- nodes[y,]
  buff <- gBuffer(event_node,width = bw)

  ##step1 : find all the samples in the radius
  ok_samples <- spatial_request(buff,sample_tree,samples)
  if(nrow(ok_samples)==0){
    return(samples_k)
  }

  ##step2 : find all the edges
  ok_edges <- spatial_request(buff,edges_tree,edges)$edge_id
  ##Step3 : for each edge, find the two vertices
  vertices <- igraph::ends(graph,ok_edges,names=F)

  ##calculate the two paths for each vertex
  un_vertices <- unique(c(vertices[,1],vertices[,2]))
  paths <- shortest_paths(graph,from=y,to=un_vertices)
  dist1 <- as.numeric(distances(graph,y,to=un_vertices))


  ##calculating alpha !
  alphas <- sapply(paths$vpath,function(a_path){
    path_vert <- as.numeric(a_path)
    vert_neighbours <- adjacent_vertices(graph,path_vert,mode = "out")
    n_neighbours <- as.numeric(sapply(vert_neighbours,length))-1
    return(prod(n_neighbours))
  })

  alphas <- ifelse(alphas == 0, 1 ,alphas)


  dist_table <- data.frame("vertex"=un_vertices,
                           "distance" = dist1,
                           "alpha" = alphas)
  ##aggregate all the data
  df_edges <- data.frame("edge_id" = ok_edges,
                         "node1" = vertices[,1],
                         "node2" = vertices[,2]
  )

  merge1 <- left_join(df_edges,dist_table,by=c("node1"="vertex"))
  df_edges$d1 <- merge1$distance
  df_edges$alpha1 <- merge1$alpha

  merge2 <-  left_join(df_edges,dist_table,by=c("node2"="vertex"))
  df_edges$d2 <- merge2$distance
  df_edges$alpha2 <- merge2$alpha


  ## now, we will calculate for each sample, the minimum distance for the two distances
  df1 <- df_edges[c("edge_id","node1","d1","alpha1")]
  df2 <- df_edges[c("edge_id","node2","d2","alpha2")]

  ## getting the first part of the distances
  df_samples1 <-  left_join(ok_samples@data,df1,by="edge_id")
  df_samples2 <-  left_join(ok_samples@data,df2,by="edge_id")

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

  ok_distances <- ifelse(full_dist1 < full_dist2, full_dist1, full_dist2)
  ok_alphas <- ifelse(full_dist1 < full_dist2, df_samples1$alpha1, df_samples2$alpha2)


  ##and then, the final distance

  k <- kernel_func(ok_distances,bw)
  k <- k * (1/ok_alphas)
  samples_k[df_samples1$oid] <- k
  return(samples_k)
}




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
