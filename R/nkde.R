#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### available kernels ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' calcualte the density at a point with a quartic function given
#' the distances of the points around, their weights and the kernel range
#'
#' @param d A vector of distances
#' @param w A vector of weights
#' @param r A float representing the range
#' @return The kernel density
#' @examples
#' #This is an internal function, no example provided
quartic_kernel <- function(d,w,r){
  w <- w[d<r]
  d <- d[d<r]
  u <- (d/r)
  K <- (3/pi) * w * (1-u**2)**2
  tot <- (1/r**2) * sum(K)
  return(sum(w)*tot)
}

#' calcualte the density at a point with a gaussian function given
#' the distances of the points around, their weights and the kernel range
#'
#' @param d A vector of distances
#' @param w A vector of weights
#' @param r A float representing the range
#' @return The kernel density
#' @examples
#' #This is an internal function, no example provided
gaussian_kernel <- function(d,w,r){
  w <- w[d<r]
  d <- d[d<r]
  u <- (d/r)
  K <- (1/(sqrt(2*pi)))*exp(-(1/2)*u**2)
  tot <- (1/r**2) * sum(K)
  return(sum(w)*tot)
}

#' calcualte the density at a point with a epanechnikov function given
#' the distances of the points around, their weights and the kernel range
#'
#' @param d A vector of distances
#' @param w A vector of weights
#' @param r A float representing the range
#' @return The kernel density
#' @examples
#' #This is an internal function, no example provided
epanechnikov_kernel <- function(d,w,r){
  w <- w[d<r]
  d <- d[d<r]
  u <- (d/r)
  K <- (3/4)*(1-u**2)
  tot <- (1/r**2) * sum(K)
  return(sum(w)*tot)
}

#' select the right kernel function with its name
#'
#' @param name The name of the kernel to use
#' @return A kernel function
#' @examples
#' #This is an internal function, no example provided
select_kernel<-function(name){
  if(name=="quartic"){
    return(quartic_kernel)
  }else if (name=="gaussian"){
    return(gaussian_kernel)
  }else if (name=="epanechnikov"){
    return(epanechnikov_kernel)
  }else{
    stop("the specified kernel is not implemented. The kernel name must be in c('quartic','gaussian','epanechnikov')")
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### NKDE calculation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate the kernel density of a set of vertices in a graph
#'
#' @param graph The graph to use to calculate the densities
#' @param origins The indices of the vertices of the graph for which the
#' densities will be calculated
#' @param destinations The indices of the vertices to use as events for
#' the kernel density calculation
#' @param range The range value of the kernel to use
#' @param kernel The name of the kernel function to use (must be one of
#' quartic, gaussian or epanechnikov)
#' @param weights A vector of weights for the destination. Default is NULL
#' if NULL, all the destinations have a weight of 1.
#' @return Vector of kernel densities (one for each origin)
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' #This is an internal function, no example provided
calc_NKDE <- function(graph,origins,destinations,range,kernel="quartic",weights=NULL){
  #selecting the kernel function
  kernel_fun <- select_kernel(kernel)
  #adjusting the weights if needed
  if (is.null(weights)){
    weights <- rep(1,length(destinations))
  }
  df <- data.frame(dest = destinations, weights = weights)
  counts <- df %>% dplyr::group_by(dest) %>% dplyr::summarise_all(sum)
  destinations <- counts$dest
  weights <- counts$weights
  i<-0
  Values <- sapply(origins,function(o){
    i<i+1
    alldistances <- igraph::distances(graph,v=o,to = destinations, mode="out")
    d <- alldistances[alldistances<range]
    w <- weights[alldistances<range]
    v <- kernel_fun(d,w,range)
    return(v)
  })
  return(Values)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### perform nkde analysis ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Lixelize a set of SpatiaLinesDataFrame, snap events (SpatialPointsDataFrame)
#' to the lixels, generate a graph and then calculate the kernel density
#' estimate for each lixels centroids
#'
#' @param lines The SpatialLinesDataFrame to use as a network
#' @param points The events to use in the kernel density estimation
#' @param snap_dist The maximum distance snapping between points and lines
#' @param lx_length The normal length of a lixel
#' @param kernel_range The range of the kernel function
#' @param kernel The name of the kernel function to use (must be one of
#' quartic, gaussian or epanechnikov). Default is Quartic.
#' @param tol The tolerence for topological operations
#' @param digits The number of digits to keep in the coordinates of the
#' geometries. This simplification is used to reduce chances of topological
#' errors
#' @param mindist The minimum length of a lixel. Defaut is NULL. If NULL,
#' mindist is set as lx_length/10
#' @param weights The name of a column of the SpatialPointsDataFrame
#' containing the weight of each point. Default is NULL. If NULL then each
#' point has a weight of 1.
#' @param show_progress A boolean indicating if a progress bar must be displayed
#' @return Vector of kernel densities (one for each origin)
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid
#' @export
#' @examples
#' data(mtl_network)
#' lixels <- nkde(mtl_network, snap_dist = 150,
#'       lx_length = 150, kernel_range = 800, mindist=50)
nkde <- function(lines,points,snap_dist,lx_length,kernel_range,kernel="quartic",weights=NULL,tol=0.1,digits=3,mindist=NULL, show_progress=T){

  #adjusting the weights if needed
  if (is.null(weights)){
    W <- rep(1,nrow(points))
  }else{
    W <- points[[weights]]
  }
  points$tmpweight <- W

  print("generating the lixels...")
  lixels <- lixelize_lines(lines,lx_length = lx_length,mindist=mindist, show_progress=show_progress)

  print("extracting the centroids of the lixels...")
  lixelscenters <- gCentroid(lixels,byid = T)
  lixelscenters <- SpatialPointsDataFrame(lixelscenters,lixels@data)

  print("combining lixels and events ...")
  lixelscenters$type <- "lixel"
  lixelscenters$OID <- 1:nrow(lixelscenters)
  points$type <- "event"
  points$OID <- 1:nrow(points)
  allpts <- rbind(lixelscenters[c("type","OID")],points[c("type","OID")])

  print("snapping points on lines...")
  snapped_points <- maptools::snapPointsToLines(allpts, lines,maxDist = snap_dist)
  snapped_points$spOID <- sp_char_index(coordinates(snapped_points),digits = digits)

  print("edditing vertices of lines...")
  newlines <- add_vertices_lines(lines,snapped_points,tol=tol, show_progress=show_progress)

  print("converting polylines to simple lines...")
  lines2 <- simple_lines(newlines)

  print("building the graph...")
  ListNet <- build_graph(lines2,digits = digits)
  Graph <- ListNet$graph

  print("getting the starting points...")
  start_pts <- subset(snapped_points,snapped_points$type=="lixel")
  origins <- find_vertices(ListNet$spvertices,start_pts,tol = tol,digits = digits)

  print("getting the destination points...")
  end_pts <- subset(snapped_points,snapped_points$type=="event")
  destinations <- find_vertices(ListNet$spvertices,end_pts,tol = tol,digits = digits)

  print("calulating the density values...")
  Values <- calc_NKDE(Graph,origins,destinations, kernel=kernel ,range = kernel_range, weights = points$tmpweight)
  start_pts$density <- Values

  lixels$density <- Values

  return(list("lixel_lines"=lixels,"lixels_points"=start_pts))
}


#' Gridded version of nkde function.
#'
#' The dataset is spatialy splitted in squares convering the extent of the
#' SpatiaLine. This greatly improve performance and memory usage. To avoid
#' frontier effect, a buffer is used (size of the kernel range parameter)
#' around the square when building local graphs. As a consequence, too small
#' squares will lead to longer calculation.
#'
#' @param lines The SpatialLinesDataFrame to use as a network
#' @param points The events to use in the kernel density estimation
#' @param snap_dist The maximum distance snapping between points and lines
#' @param lx_length The expected length of a lixel
#' @param kernel_range The range of the kernel function
#' @param kernel The name of the kernel function to use (must be one of
#' quartic, gaussian or epanechnikov). Default is Quartic
#' @param tol The tolerence for topological operations
#' @param digits The number of digits to keep in the coordinates of the
#' geometries
#' @param mindist The minimum length of a lixel. Defaut is NULL. If NULL,
#' mindist is set as lx_length/10
#' @param weights The name of a column of the SpatialPointsDataFrame
#' containing the weight of each point. Default is NULL. If NULL then each
#' point has a weight of 1.
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(2,2)
#' @param show_progress A boolean indicating if a plot should be displayed to
#' track the progressing of the algorithm.
#' @return Vector of kernel densities (one for each origin)
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @export
#' @examples
#' data(mtl_network)
#' data(bike_accidents)
#' lixels_nkde <- nkde_grided(mtl_network, bike_accident,
#'       snap_dist = 150,
#'       lx_length = 150,
#'       mindist = 50,
#'       kernel_range = 800,
#'       kernel="quartic",
#'       grid_shape=c(5,5)
#' )
nkde_grided <- function(lines,points,snap_dist,lx_length,kernel_range, kernel="quartic",tol=0.1,digits=3,mindist=NULL,weights=NULL,grid_shape=c(2,2), show_progress=TRUE){
  #adjusting the weights if needed
  if (is.null(weights)){
    W <- rep(1,nrow(points))
  }else{
    W <- points[[weights]]
  }
  points$tmpweight <- W
  print("generating the grid...")
  grid <- build_grid(grid_shape,spatial=lines)
  print("generating the lixels...")
  lixels <- lixelize_lines(lines,lx_length = lx_length, mindist=mindist)
  lixels$lxid <- 1:nrow(lixels)
  print("starting the calculation on the grid")
  remaininglixels <- lixels
  nquadra <- length(grid)

  if(show_progress){
    plot(grid)
    plot(lines,add=T)
  }

  values <- lapply(1:length(grid),function(i){
    print(paste("-------------------Iterating on quadra : ",i,"/",nquadra,"----------------",sep=""))
    #extracting the analyzed pixels
    quadra <- grid[i]
    if (show_progress){
      plot(quadra,col="red",add=T)
    }
    test_lixels <- as.vector(gIntersects(remaininglixels,quadra,byid=T))
    if(any(test_lixels)==FALSE){
      print("---- passing, empty quadra----")
      return()
    }else{
      print(paste("---quadra with values : -----",i,sep=""))
      selected_lixels <- subset(remaininglixels,test_lixels)
      if (nrow(selected_lixels)>0){
        remaininglixels <<- subset(remaininglixels,test_lixels==F)
        #extracting the needed lines
        ext <- raster::extent(selected_lixels)
        poly <- as(ext, 'SpatialPolygons')
        raster::crs(poly) <- raster::crs(selected_lixels)
        #selecting the lines close to build the network
        buff <- gBuffer(poly,width=kernel_range)
        test_lines <- as.vector(gIntersects(lines,buff,byid=T))
        selected_lines <- subset(lines,test_lines)
        #extracting the needed points
        test_points <- as.vector(gIntersects(points,buff,byid=T))
        selected_points <- subset(points,test_points)
        ##performing the real job
        print("extracting the centroids of the lixels...")
        #step3 : extraire le centroid ce ces lixels
        lixelscenters <- gCentroid(selected_lixels,byid = T)
        lixelscenters <- SpatialPointsDataFrame(lixelscenters,selected_lixels@data)
        print("combining lixels and events ...")
        #step4 : combiner les evenements et les centres de lixels
        lixelscenters$type <- "lixel"
        lixelscenters$OID <- 1:nrow(lixelscenters)

        if(any(test_points)==FALSE){
          selected_lixels$density<-0
          return(selected_lixels)
        }

        selected_points$type <- "event"
        selected_points$OID <- 1:nrow(selected_points)
        allpts <- rbind(lixelscenters[c("type","OID")],selected_points[c("type","OID")])

        #step 5 : snapper ces points sur les lignes
        print("snapping points on lines...")
        snapped_points <- maptools::snapPointsToLines(allpts, selected_lines,maxDist = snap_dist)
        snapped_points$spOID <- sp_char_index(coordinates(snapped_points),digits = digits)

        #step6 : ajouter ces points comme vertices
        print("edditing vertices of lines...")
        newlines <- add_vertices_lines(selected_lines,snapped_points,tol=tol)


        #step7 : couversion vers des lignes simple
        print("converting polylines to simple lines...")
        lines2 <- simple_lines(newlines)

        print("building the graph...")
        #step8 : generation et dessin du graph
        ListNet <- build_graph(lines2,digits = digits)
        Graph <- ListNet$graph

        print("getting the starting points...")
        #step9 : retrouver les points de depart
        start_pts <- subset(snapped_points,snapped_points$type=="lixel")
        origins <- find_vertices(ListNet$spvertices,start_pts,tol = tol,digits = digits)

        print("getting the destination points...")
        end_pts <- subset(snapped_points,snapped_points$type=="event")
        destinations <- find_vertices(ListNet$spvertices,end_pts,tol = tol,digits = digits)

        print("calulating the density values...")
        Values <- calc_NKDE(Graph,origins,destinations,range = kernel_range, kernel=kernel, weights=selected_points$tmpweight)
        selected_lixels$density <- Values

        return(selected_lixels)
      }
    }
  })

  okvalues <- values[lengths(values) != 0]
  alllixels <- do.call("rbind",okvalues)
  return(alllixels)
}




