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
  list_kernels <- c("triangle", "gaussian", "scaled gaussian", "tricube",
                    "cosine" ,"triweight", "quartic", 'epanechnikov',
                    'uniform')
  if((name %in% list_kernels)==FALSE){
    allkernels <- paste(list_kernels,collapse = ", ")
    stop(paste("the name of the kernel function must be one of",
               allkernels, sep=" "))
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
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#' @title Gamma parameter for Abramson’s adaptive bandwidth
#'
#' @description Function to calculate the gamma parameter in Abramson’s
#'   smoothing regimen.
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
    samples_k <- ess_kernel(graph,y,bw, kernel_func, samples,
                            nodes, edges, sample_tree, edges_tree)
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
#' @param samples a SpatialPointsDataFrame representing the sampling points. The
#'   samples must be snapped on the network. A column edge_id must indicate for
#'   each sample on which edge it is snapped.
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
  vertices <- ends(graph,ok_edges,names = FALSE)
  ## step4 calculate the distance between the start node and each edge vertex
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
