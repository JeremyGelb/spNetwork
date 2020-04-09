#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### worker functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Worker function for the nkde_grided.mc function.
#' @param i The row index of the spatiale grid to use
#' @param grid The SpatialPolygons object reprsenting the grid for splitted
#' calculation
#' @param lines The SpatialLinesDataFrame to use as a network
#' @param points The events to use in the kernel density estimation
#' @param kernel_range The range of the kernel function
#' @param kernel The name of the kernel function to use (must be one of
#' quartic, gaussian or epanechnikov). Default is Quartic
#' @param snap_dist The maximum distance snapping between points and lines
#' @param tol The tolerence for topological operations
#' @param digits The number of digits to keep in the coordinates of the
#' geometries
#' @param weights A vector of weights for the events. Default is NULL
#' if NULL, all the events have a weight of 1.
#' track the progressing of the calculation.
#' @return A SpatialLinesDataFrame of the lines and their kernel densities
#' @importFrom sp coordinates  SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @examples
#' #This is an internal function, no example provided
exe_nkde <- function(i,grid,lines,lixels,points,kernel,kernel_range,snap_dist,digits,tol){
  #extracting the analyzed pixels
  quadra <- grid[i]
  test_lixels <- as.vector(gIntersects(lixels,quadra,byid=T))
  if(any(test_lixels)==FALSE){
    return()
  }else{
    selected_lixels <- subset(lixels,test_lixels)
    if (nrow(selected_lixels)>0){
      #extracting the needed lines
      ext <- raster::extent(selected_lixels)
      poly <- as(ext, 'SpatialPolygons')
      raster::crs(poly) <- raster::crs(selected_lixels)
      buff <- gBuffer(poly,width=kernel_range)
      test_lines <- as.vector(gIntersects(lines,buff,byid=T))
      selected_lines <- subset(lines,test_lines)
      #extracting the needed points
      test_points <- as.vector(gIntersects(points,buff,byid=T))
      selected_points <- subset(points,test_points)
      ##performing the real job
      #step3 : extraire le centroid ce ces lixels
      lixelscenters <- gCentroid(selected_lixels,byid = T)
      lixelscenters <- SpatialPointsDataFrame(lixelscenters,selected_lixels@data)
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
      snapped_points <- maptools::snapPointsToLines(allpts, selected_lines,maxDist = snap_dist)
      snapped_points$spOID <- sp_char_index(coordinates(snapped_points),digits = digits)

      #step6 : ajouter ces points comme vertices
      newlines <- add_vertices_lines(selected_lines,snapped_points,tol=tol, show_progress=F)


      #step7 : couversion vers des lignes simple
      lines2 <- simple_lines(newlines)

      #step8 : generation et dessin du graph
      ListNet <- build_graph(lines2,digits = digits)
      Graph <- ListNet$graph

      #step9 : retrouver les points de depart
      start_pts <- subset(snapped_points,snapped_points$type=="lixel")
      origins <- find_vertices(ListNet$spvertices,start_pts,tol = tol,digits = digits)

      end_pts <- subset(snapped_points,snapped_points$type=="event")
      destinations <- find_vertices(ListNet$spvertices,end_pts,tol = tol,digits = digits)

      Values <- calc_NKDE(Graph,origins,destinations,kernel=kernel,range = kernel_range, weights=selected_points$tmpweight)
      selected_lixels$density <- Values
      return(selected_lixels)
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### launching functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Gridded version of nkde with multicore support.
#'
#' The dataset is spatialy splitted in squares convering
#' the extent of the SpatiaLine. This greatly improve performance and memory
#' usage. To avoid frontier effect, a buffer is used (twice the size of the
#' range parameter) around the square when building local graphs. As a
#' consequence, too small squares will lead to longer calculation.
#' This version of the function can use multicore if a plan is defined
#' with the future package. If you do not use multicore, the regular function
#' nkde_grided will be slightely more efficient.
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
#' mindist is set as lixel_length/10
#' @param weights A vector of weights for the events. Default is NULL
#' if NULL, all the events have a weight of 1.
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(2,2)
#' @param show_progress A boolean indicating if a plot should be displayed to
#' track the progressing of the calculation.
#' @return Vector of kernel densities (one for each origin)
#' @importFrom sp coordinates  SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @export
#' @examples
#' data(mtl_network)
#' data(bike_accidents)
#' future::plan(future::multiprocess(workers=2))
#' lixels_nkde <- nkde_grided.mc(mtl_network, bike_accident,
#'       snap_dist = 150,
#'       lx_length = 150,
#'       mindist = 50,
#'       kernel_range = 800,
#'       kernel="quartic",
#'       grid_shape=c(5,5)
#' )
nkde_grided.mc <- function(lines,points,snap_dist,lx_length,kernel_range, kernel="quartic",tol=0.1,digits=3,mindist=NULL,weights=NULL,grid_shape=c(2,2), show_progress=TRUE){
  #adjusting the weights if needed
  if (is.null(weights)){
    W <- rep(1,nrow(points))
  }else{
    W <- points[[weights]]
  }
  points$tmpweight <- W
  print("generating the grid...")
  grid <- build_grid(grid_shape,spatial=lines)
  ##step2 : generer les lixels
  print("generating the lixels...")
  lixels <- lixelize_lines(lines,lx_length = lx_length, mindist=mindist)
  lixels$lxid <- 1:nrow(lixels)
  ##step3 : lancer les iterations
  print("starting the calculation on the grid")
  remaininglixels <- lixels
  iseq <- 1:length(grid)
  if (show_progress){
    progressr::with_progress({
      p <- progressr::progressor(along = iseq)
      values <- future.apply::future_lapply(iseq,function(i){
        p(sprintf("i=%g", i))
        return(exe_nkde(i,grid,lines,lixels,points,kernel,kernel_range,snap_dist,digits,tol))
      })
    })
  }else{
    values <-  future.apply::future_lapply(iseq,exe_nkde,grid=grid,lines=lines,
                     lixels=lixels,points=points,
                     kernel_range=kernel_range, kernel=kernel,
                     snap_dist=snap_dist,digits=digits,tol=tol)
  }
  print("Combining the results from processes ...")
  okvalues <- values[lengths(values) != 0]
  alllixels <- do.call("rbind",okvalues)
  alllixels <- alllixels[!duplicated(alllixels$lxid), ]
  return(alllixels)
}

