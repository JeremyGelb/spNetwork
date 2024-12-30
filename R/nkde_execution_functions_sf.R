#' @importFrom Rdpack reprompt

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### helper functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Projection test
#' @description Check if a feature collection is in a projected CRS
#' @param obj A feature collection
#' @return A boolean
#' @importFrom sf st_is_longlat st_crs st_crs<-
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
is_projected <- function(obj) {
  !is.na(st_crs(obj)) &
    !st_is_longlat(obj) &
    is.null(st_crs(obj)$to_meter)
}

#' @title Geometry sanity check
#'
#' @description Function to check if the geometries given by the user are valid.
#'
#' @param lines A feature collection of lines
#' @param samples A feature collection of points (the samples)
#' @param events A feature collection of points (the events)
#' @param study_area A feature collection of polygons (the study_area)
#' @return TRUE if all the checks are passed
#' @importFrom sf st_geometry_type st_is_valid
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
check_geometries <- function(lines,samples,events, study_area){ # nocov start

  # checking if geometries are all valid, simple and planar
  obj_names <- c("lines","samples","events")
  objs <- list(lines,samples,events)
  for(i in 1:length(obj_names)){
    obj <- objs[[i]]
    obj_name <- obj_names[[i]]
    types <- st_geometry_type(obj)
    test <- grepl("MULTI", types, fixed = TRUE)
    if(any(test)){
      stop(paste("the ",obj_name," must be simple geometries,
                 considere using the function sf::st_cast to convert them to simple geometries",sep=""))
    }
    if(any(st_is_valid(obj) == FALSE)){
      stop(paste("the ",obj_name," must be valid geometries",sep=""))
    }
    if(is_projected(obj)==FALSE){
      stop(paste("the ",obj_name," must be projected (planar coordinates)",sep=""))
    }
  }
  # checking if the geometries type are good
  if(class(events)[[1]]!="sf" | unique(st_geometry_type(events))!="POINT"){
    stop("the events must be given as a feature class of points")
  }
  if(class(samples)[[1]]!="sf" | unique(st_geometry_type(samples))!="POINT"){
    stop("the samples must be given as a feature class of points")
  }
  if(class(lines)[[1]]!="sf" | unique(st_geometry_type(lines))!="LINESTRING"){
    stop("the lines must be given as a feature class of linestrings")
  }

  # checking if the CRS are good
  if(is.null(study_area)){
    comp <- c(st_crs(samples) == st_crs(events),
              st_crs(lines) == st_crs(events),
              st_crs(lines) == st_crs(samples))
  }else{
    comp <- c(st_crs(study_area) == st_crs(events),
              st_crs(study_area) == st_crs(samples),
              st_crs(study_area) == st_crs(lines),
              st_crs(samples) == st_crs(events),
              st_crs(lines) == st_crs(events),
              st_crs(lines) == st_crs(samples))
  }

  if(any(comp==FALSE)){
    stop("the lines, events and samples must have the same Coordinates Reference System (crs)")
  }
  return(TRUE)
}

#defining some global variables (weird felx but ok)
utils::globalVariables(c("spid", "weight", ".", "bws",'X','Y','wid','m_X','m_Y')) # nocov end


#' @title Clean events geometries
#'
#' @description Function to avoid having events at the same location.
#'
#' @param events The feature collection of points to contract (must have a weight column)
#' @param digits The number of digits to keep
#' @param agg A double indicating if the points must be aggregated within a distance.
#' if NULL, then the points are aggregated by rounding the coordinates.
#' @return A new feature collection of points
#' @importFrom data.table tstrsplit setDF
#' @importFrom sf st_coordinates
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
clean_events <- function(events,digits=5,agg=NULL){

  base_crs <- st_crs(events)

  if(is.null(agg)){
    events$spid <- sp_char_index(st_coordinates(events),digits)

    if("time" %in% names(events)){
      events$spid <- paste0(events$spid,"_",events$time)
    }

    edf <- st_drop_geometry(events)
    new_events <- data.table(edf)
    n1 <- names(new_events)[names(new_events) %in% c("weight","spid") == FALSE]


    agg_events <- new_events[, lapply(.SD, max) , by = .(spid),  .SDcols = n1]
    agg_events$weight <- new_events[, .(sum(weight)), by = .(spid)][,2]
    if("time" %in% names(events)){
      agg_events[,  c("X", "Y","time") := tstrsplit(spid, "_", fixed = TRUE)]
      agg_events$time <- as.numeric(agg_events$time)
    }else{
      agg_events[,  c("X", "Y") := tstrsplit(spid, "_", fixed = TRUE)]
    }


    agg_events$X <- as.numeric(agg_events$X)
    agg_events$Y <- as.numeric(agg_events$Y)
    new_events <- setDF(agg_events)
    new_events <- st_as_sf(new_events,
                           coords = c("X", "Y"),
                           crs = st_crs(events))
  }else{
    if("time" %in% names(events)){
      XY <- st_coordinates(events)
      events$X <- XY[,1]
      events$Y <- XY[,2]
      events$gp <- aggregate_points(events,agg, return_ids = TRUE)

      if(sum(events$gp) == 0){
        new_events <- events
        new_events$gp <- NULL
        new_events$spid <- sp_char_index(st_coordinates(new_events),digits)
      }else{
        events <- st_drop_geometry(events)
        part1 <- subset(events, events$gp == 0)
        part2 <- setDT(subset(events, events$gp != 0))
        mean_coords <- part2[,.(m_X = mean(X),
                              m_Y = mean(Y),
                              m_wid = data.table::first(wid)
                              ), by = "gp"]
        events2 <- setDT(merge(part2, mean_coords, by = "gp", all.x = TRUE))
        events2$m_X <- ifelse(is.na(events2$m_X), events2$X,events2$m_X)
        events2$m_Y <- ifelse(is.na(events2$m_Y), events2$X,events2$m_Y)
        new_events <- events2[, .(X = data.table::first(m_X),
                               Y = data.table::first(m_Y),
                               wid = data.table::first(wid),
                               weight = sum(weight)
        ), by = c("gp","time")]
        new_events <- rbind(new_events, part1[names(new_events)])
        new_events$gp <- NULL
        new_events <- st_as_sf(new_events, coords = c('X','Y'), crs = base_crs)

      }



    }else{
      new_events <- aggregate_points(events,agg)
      new_events$spid <- sp_char_index(st_coordinates(new_events),digits)
    }
  }

  # we can now add the goid : unique id of each location
  coords <- st_coordinates(new_events)
  new_events$goid <- as.numeric(as.factor(sp_char_index(coords,digits)))

  return(new_events)

}

#' @title Events aggregation
#'
#' @description Function to aggregate points within a radius.
#'
#' @details This function can be used to aggregate points within a radius. This
#'   is done by using the dbscan algorithm. This process is
#'   repeated until no more modification is applied.
#'
#' @param points The feature collection of points to contract (must have a weight column)
#' @param maxdist The distance to use
#' @param weight The name of the column to use as weight (default is "weight").
#' The values of the aggregated points for this column will be summed. For all
#' the other columns, only the max value is retained.
#' @param return_ids A boolean (default is FALSE), if TRUE, then an index indicating
#' for each point the group it belongs to is returned. If FALSE, then a spatial
#' point features is returned with the points already aggregated.
#' @return A new feature collection of points
#' @export
#' @examples
#' data(bike_accidents)
#' bike_accidents$weight <- 1
#' agg_points <- aggregate_points(bike_accidents, 5)
aggregate_points <- function(points, maxdist, weight = "weight", return_ids = FALSE){

  num_cols <- unlist(lapply(points, is.numeric))
  num_names <- names(points)[num_cols]

  points <- points[num_names]

  ## reimplementation avec dbscan
  coords <- st_coordinates(points)
  points$X <- coords[,1]
  points$Y <- coords[,2]
  result <- dbscan::dbscan(coords, eps = maxdist, minPts = 2)

  if(return_ids){
    return(result$cluster)
  }

  pt_list <- split(points, f = result$cluster)
  if(length(pt_list) == 1){
    all_pts <- points
  }else{
    old_pts <- pt_list[[1]]
    #mat1 <- cbind(st_coordinates(old_pts), st_drop_geometry(old_pts))
    mat1 <- st_drop_geometry(old_pts)
    new_pts <- data.frame(t(sapply(2:length(pt_list), function(i){
      pts <- pt_list[[i]]
      if(nrow(pts) == 1 ){
        return(c(st_drop_geometry(pts)))
      }else{
        coords <- colMeans(st_coordinates(pts))
        feat_max <- sapply(num_names, function(name){max(pts[[name]])})
        feat_max[names(feat_max) == weight] <- sum(pts[[weight]])
        #feat_max <- apply(st_drop_geometry(pts), 2, max)
        #feat_max[names(feat_max) == weight] <- sum(pts[[weight]])
        #feat_max <- feat_max[names(feat_max) %in% num_names]
        return(c(coords, feat_max))
      }
    })))
    new_pts <- rbind(mat1, new_pts)
    for (coln in c("X","Y", num_names)){
      new_pts[[coln]] <- as.numeric(new_pts[[coln]])
    }
    all_pts <- st_as_sf(new_pts, coords = c("X","Y"))
    all_pts[c("X","Y")] <- st_coordinates(all_pts)
    st_crs(all_pts) <- st_crs(points)
  }

  return(all_pts)
}



#' @title Prior data preparation
#'
#' @description A simple function to prepare data before the NKDE calculation.
#'
#' @param samples A feature collection of points representing the samples points
#' @param lines A feature collection of Linestrings representing the network
#' @param events A feature collection of points representing the events points
#' @param w A numeric vector representing the weight of the events
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @param agg A double indicating if the points must be aggregated within a distance.
#' if NULL, then the points are aggregated by rounding the coordinates.
#' @return the data prepared for the rest of the operations
#' @importFrom data.table tstrsplit setDF
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
prepare_data <- function(samples,lines,events, w ,digits,tol, agg){

  ## step1 cleaning the events
  events$weight <- w

  # the goid for the events is added here
  # two events at the same location will have the same goid
  # but the events with different time stamp must not be aggregated
  # geographically
  events <- clean_events(events,digits,agg)

  ## step2 defining the global IDS
  #events$goid <- 1:nrow(events)
  samples$goid <- 1:nrow(samples)
  samples <- samples[c("goid")]

  # removing existing X and Y coordinates
  if("X" %in% names(samples)){
    samples$X <- NULL
  }
  if("Y" %in% names(samples)){
    samples$Y <- NULL
  }
  if("X" %in% names(events)){
    events$X <- NULL
  }
  if("Y" %in% names(events)){
    events$Y <- NULL
  }

  ## step3 remove lines with no length
  lines$length <- as.numeric(st_length(lines))
  lines <- subset(lines, lines$length>0)
  lines$oid <- 1:nrow(lines)

  return(list("samples" = samples,
              "lines" = lines[c("length","oid")],
              "events" = events))

}





#' @title Split data with a grid
#'
#' @description Function to split the dataset according to a grid.
#'
#' @param grid A spatial grid to split the data within
#' @param samples A feature collection of points representing the samples points
#' @param events A feature collection of points representing the events points
#' @param lines A feature collection of linestrings representing the network
#' @param bw The kernel bandwidth (used to avoid edge effect)
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @param split_all A boolean indicating if we must split the lines at each vertex
#' (TRUE) or only at event vertices (FALSE)
#' @return A list with the split dataset
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
split_by_grid <- function(grid,samples,events,lines,bw,tol, digits, split_all = TRUE){

  grid$grid_id <- as.character(grid$oid)
  grid$oid <- NULL

  ## step 1 : preparing the spatial intersections
  inter_samples <- sf::st_join(samples, grid)
  inter_samples <- split(inter_samples,inter_samples$grid_id)
  inter_events <- sf::st_join(events, st_buffer(grid, dist=bw))
  inter_events <- split(inter_events,inter_events$grid_id)

  # pour les lignes on doit souffrir un peu plus
  inter_lines <- sf::st_join(lines, st_buffer(grid, dist=(bw+0.5*bw)))
  inter_lines <- split(inter_lines,inter_lines$grid_id)

  # inter_lines <- sf::st_intersects(lines, st_buffer(grid, dist=(bw+0.5*bw)))
  # idxs_lines <- rep(1:nrow(inter_lines), lengths(inter_lines))
  # idxs_gid <- unlist(inter_lines)
  # inter_lines <- lines[idxs_lines,]
  # inter_lines$grid_id <- idxs_gid

  ## step2 : split the datasets

  selections <- lapply(1:nrow(grid),function(i){

    gid <- grid$grid_id[[i]]

    # selecting the samples in the grid
    sel_samples <- inter_samples[[gid]]
    sel_samples$grid_id <- NULL

    # if there is no sampling points in the rectangle, then return NULL
    if(is.null(sel_samples)){
      return(NULL)
    }

    # selecting the events
    sel_events <- inter_events[[gid]]
    sel_events$grid_id <- NULL

    # selecting the lines in a buffer
    sel_lines <- inter_lines[[gid]]
    sel_lines$grid_id <- NULL

    sel_lines$oid <- 1:nrow(sel_lines)

    # snapping the events on the lines
    if(is.null(sel_events)){
      new_lines <- sel_lines
    }else{
      a <- nrow(sel_events)
      b <- nrow(sel_lines)
      x <-  a*b
      # if(is.na(x)){
      #   stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      # }
      # if(x >= 2*10^9){
      #   stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      # }
      sel_events <- snapPointsToLines2(sel_events,sel_lines,idField = "oid", snap_dist = bw)

      if(split_all){
        new_lines <- add_vertices_lines(sel_lines,sel_events,sel_events$nearest_line_id,tol)
      }else{
        new_lines <- split_lines_at_vertex(sel_lines,sel_events,sel_events$nearest_line_id,tol)
      }

    }

    # split lines at events
    if(split_all){
      new_lines <- simple_lines(new_lines)
    }
    new_lines$length <- as.numeric(st_length(new_lines))
    new_lines <- subset(new_lines,new_lines$length>0)

    # remove lines that are loops
    new_lines <- remove_loop_lines(new_lines,digits)
    new_lines$oid <- 1:nrow(new_lines)
    new_lines <- new_lines[c("length","oid")]


    return(list("samples" = sel_samples,
                "events" = sel_events,
                "lines" = new_lines))
  })
  #let us remove all the empty quadra
  selections <- selections[lengths(selections) != 0]

  for(i in 1:length(selections)){
    selections[[i]]$index <- i
  }

  return(selections)

}



#' @title Split data with a grid for the adaptive bw function
#'
#' @description Function to split the dataset according to a grid for the
#' adaptive bw function.
#'
#' @param grid A spatial grid to split the data within
#' @param events A feature collection of points representing the events
#' @param lines A feature collection of lines representing the network
#' @param bw The kernel bandwidth (used to avoid edge effect)
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @return A list with the split dataset
#' @importFrom sf st_buffer
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
split_by_grid_abw <- function(grid,events,lines,bw,tol,digits){

  grid$grid_id <- as.character(grid$oid)
  grid$oid <- NULL

  ## step1 : calculating the spatial joins
  inter_events_rect <- st_join(events, grid)
  inter_events_rect <- split(inter_events_rect,inter_events_rect$grid_id)

  inter_events <-  st_join(events, st_buffer(grid, dist=bw))
  inter_events <- split(inter_events,inter_events$grid_id)

  inter_lines <- st_join(lines, st_buffer(grid,dist=(bw+0.5*bw)))
  inter_lines <- split(inter_lines,inter_lines$grid_id)

  ## step2 : split the datasets

  selections <- lapply(1:nrow(grid),function(i){

    gid <- grid$grid_id[[i]]

    # selecting the events in the grid
    sel_events_rect <- inter_events_rect[[gid]]

    # if there is no event points in the rectangle, then return NULL
    if(is.null(sel_events_rect)){
      return(NULL)
    }
    sel_events <- inter_events[[gid]]

    # selecting the lines in a buffer
    sel_lines <- inter_lines[[gid]]
    sel_lines$oid <- 1:nrow(sel_lines)

    # snapping the events on the lines
    if(is.null(sel_events)){
      new_lines <- sel_lines
    }else{
      a <- nrow(sel_events)
      b <- nrow(sel_lines)
      x <-  a*b
      if(is.na(x)){
        stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      }
      if(x >= 2*10^9){
        stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      }
      sel_events <- snapPointsToLines2(sel_events,sel_lines,idField = "oid")
      new_lines <- add_vertices_lines(sel_lines,sel_events,sel_events$nearest_line_id,tol)
    }

    # split lines at events
    new_lines <- simple_lines(new_lines)
    new_lines$length <- as.numeric(st_length(new_lines))
    new_lines <- subset(new_lines,new_lines$length>0)

    # remove lines that are loops
    new_lines <- remove_loop_lines(new_lines,digits)
    new_lines$oid <- 1:nrow(new_lines)
    new_lines <- new_lines[c("length","oid")]

    sel_events_rect <- subset(sel_events, sel_events$goid %in% sel_events_rect$goid)

    return(list("samples" = sel_events_rect,
                "events" = sel_events,
                "lines" = new_lines))
  })
  #let us remove all the empty quadra
  selections <- selections[lengths(selections) != 0]

  for(i in 1:length(selections)){
    selections[[i]]$index <- i
  }

  return(selections)

}



#' @title Split data with a grid
#'
#' @description Function to split the dataset according to a grid.
#'
#' @param grid A spatial grid to split the data within
#' @param samples A feature collection of points representing the samples points
#' @param events A feature collection of points representing the events points
#' @param lines A feature collection of linestrings representing the network
#' @param bw The kernel bandwidth (used to avoid edge effect)
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @param split_all A boolean indicating if we must split the lines at each vertex
#' (TRUE) or only at event vertices (FALSE)
#' @return A list with the split dataset
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
split_by_grid.mc <- function(grid,samples,events,lines,bw,tol, digits, split_all = TRUE){


  grid$grid_id <- as.character(grid$oid)
  grid$oid <- NULL

  ## step 1 : preparing the spatial intersections (done on a single core for the moment)
  inter_samples <- sf::st_join(samples, grid)
  inter_samples <- split(inter_samples,inter_samples$grid_id)
  inter_events <- sf::st_join(events, st_buffer(grid, dist=bw))
  inter_events <- split(inter_events,inter_events$grid_id)

  inter_lines <- sf::st_join(lines, st_buffer(grid, dist=(bw+0.5*bw)))
  inter_lines <- split(inter_lines,inter_lines$grid_id)

  dfs <- lapply(1:nrow(grid), function(i){

    gid <- grid$grid_id[[i]]
    # selecting the samples in the grid
    sel_samples <- inter_samples[[gid]]
    sel_samples$grid_id <- NULL
    # selecting the events
    sel_events <- inter_events[[gid]]
    sel_events$grid_id <- NULL
    # selecting the lines in a buffer
    sel_lines <- inter_lines[[gid]]
    sel_lines$grid_id <- NULL
    sel_lines$oid <- 1:nrow(sel_lines)
    return(
      list(i, sel_samples, sel_events, sel_lines)
    )

  })

  ## step2 : split the datasets (done on multiple cores)

  selections <- future.apply::future_lapply(dfs,function(df){

    i <- df[[1]]
    sel_samples <- df[[2]]
    sel_events <- df[[3]]
    sel_lines <- df[[4]]

    # if there is no sampling points in the rectangle, then return NULL
    if(is.null(sel_samples)){
      return(NULL)
    }

    # snapping the events on the lines
    if(is.null(sel_events)){
      new_lines <- sel_lines
    }else{
      a <- nrow(sel_events)
      b <- nrow(sel_lines)
      x <-  a*b
      # if(is.na(x)){
      #   stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      # }
      # if(x >= 2*10^9){
      #   stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      # }
      sel_events <- snapPointsToLines2(sel_events,sel_lines,idField = "oid", snap_dist = bw)

      if(split_all){
        new_lines <- add_vertices_lines(sel_lines,sel_events,sel_events$nearest_line_id,tol)
      }else{
        new_lines <- split_lines_at_vertex(sel_lines,sel_events,sel_events$nearest_line_id,tol)
      }

    }

    # split lines at events
    if(split_all){
      new_lines <- simple_lines(new_lines)
    }
    new_lines$length <- as.numeric(st_length(new_lines))
    new_lines <- subset(new_lines,new_lines$length>0)

    # remove lines that are loops
    new_lines <- remove_loop_lines(new_lines,digits)
    new_lines$oid <- 1:nrow(new_lines)
    new_lines <- new_lines[c("length","oid")]


    return(list("samples" = sel_samples,
                "events" = sel_events,
                "lines" = new_lines))

  }, future.packages = c("spNetwork"))

  #let us remove all the empty quadra
  selections <- selections[lengths(selections) != 0]

  for(i in 1:length(selections)){
    selections[[i]]$index <- i
  }

  return(selections)

}



#' @title Split data with a grid
#'
#' @description Function to split the dataset according to a grid.
#'
#' @param grid A spatial grid to split the data within
#' @param samples A feature collection of points representing the samples points
#' @param events A feature collection of points representing the events points
#' @param lines A feature collection of linestrings representing the network
#' @param bw The kernel bandwidth (used to avoid edge effect)
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @param split_all A boolean indicating if we must split the lines at each vertex
#' (TRUE) or only at event vertices (FALSE)
#' @return A list with the split dataset
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
split_by_grid.mc <- function(grid,samples,events,lines,bw,tol, digits, split_all = TRUE){


  grid$grid_id <- as.character(grid$oid)
  grid$oid <- NULL

  ## step 1 : preparing the spatial intersections (done on a single core for the moment)
  inter_samples <- sf::st_join(samples, grid)
  inter_samples <- split(inter_samples,inter_samples$grid_id)
  inter_events <- sf::st_join(events, st_buffer(grid, dist=bw))
  inter_events <- split(inter_events,inter_events$grid_id)

  inter_lines <- sf::st_join(lines, st_buffer(grid, dist=(bw+0.5*bw)))
  inter_lines <- split(inter_lines,inter_lines$grid_id)

  dfs <- lapply(1:nrow(grid), function(i){

    gid <- grid$grid_id[[i]]
    # selecting the samples in the grid
    sel_samples <- inter_samples[[gid]]
    sel_samples$grid_id <- NULL
    # selecting the events
    sel_events <- inter_events[[gid]]
    sel_events$grid_id <- NULL
    # selecting the lines in a buffer
    sel_lines <- inter_lines[[gid]]
    sel_lines$grid_id <- NULL
    sel_lines$oid <- 1:nrow(sel_lines)
    return(
      list(i, sel_samples, sel_events, sel_lines)
    )

  })

  ## step2 : split the datasets (done on multiple cores)

  selections <- future.apply::future_lapply(dfs,function(df){

    i <- df[[1]]
    sel_samples <- df[[2]]
    sel_events <- df[[3]]
    sel_lines <- df[[4]]

    # if there is no sampling points in the rectangle, then return NULL
    if(is.null(sel_samples)){
      return(NULL)
    }

    # snapping the events on the lines
    if(is.null(sel_events)){
      new_lines <- sel_lines
    }else{
      a <- nrow(sel_events)
      b <- nrow(sel_lines)
      x <-  a*b
      # if(is.na(x)){
      #   stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      # }
      # if(x >= 2*10^9){
      #   stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      # }
      sel_events <- snapPointsToLines2(sel_events,sel_lines,idField = "oid", snap_dist = bw)

      if(split_all){
        new_lines <- add_vertices_lines(sel_lines,sel_events,sel_events$nearest_line_id,tol)
      }else{
        new_lines <- split_lines_at_vertex(sel_lines,sel_events,sel_events$nearest_line_id,tol)
      }

    }

    # split lines at events
    if(split_all){
      new_lines <- simple_lines(new_lines)
    }
    new_lines$length <- as.numeric(st_length(new_lines))
    new_lines <- subset(new_lines,new_lines$length>0)

    # remove lines that are loops
    new_lines <- remove_loop_lines(new_lines,digits)
    new_lines$oid <- 1:nrow(new_lines)
    new_lines <- new_lines[c("length","oid")]


    return(list("samples" = sel_samples,
                "events" = sel_events,
                "lines" = new_lines))

  }, future.packages = c("spNetwork"))

  #let us remove all the empty quadra
  selections <- selections[lengths(selections) != 0]

  for(i in 1:length(selections)){
    selections[[i]]$index <- i
  }

  return(selections)

}




#' @title Split data with a grid for the adaptive bw function (multicore)
#'
#' @description Function to split the dataset according to a grid for the
#' adaptive bw function with multicore support
#'
#' @param grid A spatial grid to split the data within
#' @param events A feature collection of points representing the events points
#' @param lines A feature collection of lines representing the network
#' @param bw The kernel bandwidth (used to avoid edge effect)
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @return A list with the split dataset
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
split_by_grid_abw.mc <- function(grid,events,lines,bw,tol,digits){

  grid$grid_id <- as.character(grid$oid)
  grid$oid <- NULL

  ## step1 : calculating the spatial joins (single core for the moment)
  inter_events_rect <- st_join(events, grid)
  inter_events_rect <- split(inter_events_rect,inter_events_rect$grid_id)

  inter_events <-  st_join(events, st_buffer(grid, dist=bw))
  inter_events <- split(inter_events,inter_events$grid_id)

  inter_lines <- st_join(lines, st_buffer(grid,dist=(bw+0.5*bw)))
  inter_lines <- split(inter_lines,inter_lines$grid_id)

  ## step2 split the datasets

  sub_samples <- lapply(1:nrow(grid),function(i){

    gid <- grid$grid_id[[i]]

    #selecting the samples in the grid
    sel_events_rect <- inter_events_rect[[gid]]

    ##return NULL if there is no sampling point in the rectangle
    if(is.null(sel_events_rect)){
      return(NULL)
    }

    #selecting the events in a buffer
    sel_events <- inter_events[[gid]]

    #selecting the lines in a buffer
    sel_lines <- inter_lines[[gid]]
    sel_lines$oid <- 1:nrow(sel_lines)

    return(list("sel_lines" = sel_lines,
                "sel_events"=sel_events,
                "sel_events_rect"=sel_events_rect))

  })


  selections <- future.apply::future_lapply(sub_samples,function(sub){
    if(is.null(sub)){
      return (NULL)
    }
    sel_lines <- sub$sel_lines
    sel_events <- sub$sel_events
    sel_events_rect <- sub$sel_events_rect

    #snapping the events on the lines
    if(is.null(sel_events)){
      new_lines <- sel_lines
    }else{
      a <- nrow(sel_events)
      b <- nrow(sel_lines)
      x <-  a*b
      if(is.na(x)){
        stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      }
      if(x >= 2*10^9){
        stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
      }
      sel_events <- snapPointsToLines2(sel_events,sel_lines,idField = "oid")
      invisible(capture.output(new_lines <- add_vertices_lines(sel_lines,sel_events,sel_events$nearest_line_id,tol)))
    }

    #split lines at events
    new_lines <- simple_lines(new_lines)
    new_lines$length <- as.numeric(st_length(new_lines))
    new_lines <- subset(new_lines,new_lines$length>0)

    # remove lines that are loops
    new_lines <- remove_loop_lines(new_lines,digits)
    new_lines$oid <- 1:nrow(new_lines)
    new_lines <- new_lines[c("length","oid")]

    sel_events_rect <- subset(sel_events, sel_events$goid %in% sel_events_rect$goid)


    return(list("samples" = sel_events_rect,
                "events" = sel_events,
                "lines" = new_lines))
  }, future.packages = c("spNetwork"))
  #let us remove the empty regions
  selections <- selections[lengths(selections) != 0]

  #randomize the elements to minimize calculation time
  selections <- sample(selections)

  for(i in 1:length(selections)){
    selections[[i]]$index <- i
  }

  return(selections)

}


#' @title Adaptive bandwidth
#'
#' @description Function to calculate Adaptive bandwidths according to Abramson’s smoothing regimen.
#'
#' @param grid A spatial grid to split the data within
#' @param events A feature collection of points representing the events points
#' @param lines A feature collection of linestrings representing the network
#' @param bw The fixed kernel bandwidth (can also be a vector, the value returned will be a matrix in that case)
#' @param trim_bw The maximum size of local bandwidths (can also be a vector, must match bw)
#' @param method The method to use when calculating the NKDE
#' @param kernel_name The name of the kernel to use
#' @param max_depth The maximum recursion depth
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @param sparse A Boolean indicating if sparse matrix should be used
#' @param verbose A Boolean indicating if update messages should be printed
#' @return A vector with the local bandwidths
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
adaptive_bw <- function(grid,events,lines,bw,trim_bw,method,kernel_name,max_depth,tol,digits,sparse,verbose){

  ##step 1 split the datas !
  selections <- split_by_grid_abw(grid,events,lines, bw = max(trim_bw), tol, digits)
  ## step 2 calculating the temp NKDE values
  if(verbose){
    print("start calculating the kernel values for the adaptive bandwidth...")
  }
  n_quadra <- length(selections)

  dfs <- lapply(1:n_quadra,function(i){
    sel <- selections[[i]]
    # bws should be a vector if we only have on global bandwidth
    if(length(bw) == 1){
      bws <- rep(bw,nrow(sel$events))
    }else{
      # and it will be a matrix if we have several bandwidths
      bws <- sapply(bw, function(x){
        rep(x, nrow(sel$events))
      })
    }

    if(verbose){
      print(paste("    quadra ",i,"/",n_quadra,sep=""))
    }
    values <- nkde_worker(lines =  sel$lines,
                          events = sel$events,
                          samples = sel$samples,
                          kernel_name = kernel_name,
                          bw = bw,
                          bws = bws,
                          method = method,
                          div = "none",
                          digits = digits,
                          tol = tol,
                          sparse = sparse,
                          max_depth = max_depth,
                          verbose = verbose)
    if(any(values==0)){
      # Adding a better way to explain the bug would be great !
      stop("This is embarassing, we should not get 0 values here (bug code 0001). Please report the bug at https://github.com/JeremyGelb/spNetwork/issues and provide the dataset used.")
    }
    # df <- data.frame("goid"=sel$samples$goid,
    #                  "k" = values)
    if(nrow(sel$samples) == 1){
      df <- c(sel$samples$goid, values)
    }else{
      df <- cbind(sel$samples$goid, values)
    }
    return(df)
  })

   ## step 3  combining the results
  tot_df <- do.call(rbind,dfs)
  #tot_df <- tot_df[order(tot_df[,1]),]
  #tot_df <- tot_df[match(tot_df[,1], events$goid),]
  tot_df <- tot_df[match(events$goid,tot_df[,1]),]
  ## step 4 calculating the new bandwidth !
  if(ncol(tot_df) == 2){
    k <- tot_df[,2]
    delta <- calc_gamma(k)
    new_bw <- bw * (k**(-1/2) * delta**(-1))
    new_bw <- ifelse(new_bw<trim_bw,new_bw,trim_bw)
  }else{
    new_bw <- sapply(2:ncol(tot_df), function(i){
      k <- tot_df[,i]
      delta <- calc_gamma(k)
      nbw <- bw[[i-1]] * (k**(-1/2) * delta**(-1))
      nbw <- ifelse(nbw<trim_bw[[i-1]],nbw,trim_bw[[i-1]])
      return(nbw)
    })
  }

  return(new_bw)
}


#' @title Adaptive bandwidth (multicore)
#'
#' @description Function to calculate Adaptive bandwidths according to
#'   Abramson’s smoothing regimen with multicore support
#'
#' @param grid A spatial grid to split the data within
#' @param events A feature collection of points representing the events
#' @param lines A feature collection of linestrings representing the network
#' @param bw The fixed kernel bandwidth (can also be a vector, the value returned will be a matrix in that case)
#' @param trim_bw The maximum size of local bandwidths (can also be a vector, must match bw)
#' @param method The method to use when calculating the NKDE
#' @param kernel_name The name of the kernel to use
#' @param max_depth The maximum recursion depth
#' @param digits The number of digits to keep
#' @param tol A float indicating the spatial tolerance when snapping events on
#' lines
#' @param sparse A Boolean indicating if sparse matrix should be used
#' @param verbose A Boolean indicating if update messages should be printed
#' @return A vector with the local bandwidths
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
adaptive_bw.mc <- function(grid,events,lines,bw,trim_bw,method,kernel_name,max_depth,tol,digits,sparse,verbose){

  ##step 1 split the datas !
  selections <- split_by_grid_abw.mc(grid,events,lines,trim_bw, tol, digits)
  ## step 2 calculating the temp NKDE values
  if(verbose){
    print("start calculating the kernel values for the adaptive bandwidth...")
  }

  if(verbose){
    print("starting the kernel calculations ...")
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {

        # bws should be a vector if we only have on global bandwidth
        if(length(bw) == 1){
          bws <- rep(bw,nrow(sel$events))
        }else{
          # and it will be a matrix if we have several bandwidths
          bws <- sapply(bw, function(x){
            rep(x, nrow(sel$events))
          })
        }
        invisible(capture.output(values <- spNetwork::nkde_worker(lines =  sel$lines,
                                                       events = sel$events,
                                                       samples = sel$samples,
                                                       kernel_name = kernel_name,
                                                       bw = bw,
                                                       bws = bws,
                                                       method = method,
                                                       div = "none",
                                                       digits = digits,
                                                       tol = tol,
                                                       sparse = sparse,
                                                       max_depth = max_depth,
                                                       verbose = verbose)

                                 ))

        # df <- data.frame("goid"=sel$samples$goid,
        #                  "k" = values)

        if(nrow(sel$samples) == 1){
          df <- c(sel$samples$goid, values)
        }else{
          df <- cbind(sel$samples$goid, values)
        }

        p(sprintf("i=%g", sel$index))
        return(df)
      }, future.packages = c("spNetwork"))
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {

      # bws should be a vector if we only have on global bandwidth
      if(length(bw) == 1){
        bws <- rep(bw,nrow(sel$events))
      }else{
        # and it will be a matrix if we have several bandwidths
        bws <- sapply(bw, function(x){
          rep(x, nrow(sel$events))
        })
      }

      values <- spNetwork::nkde_worker(lines =  sel$lines,
                            events = sel$events,
                            samples = sel$samples,
                            kernel_name = kernel_name,
                            bw = bw,
                            bws = bws,
                            method = method,
                            div = "none",
                            digits = digits,
                            tol = tol,
                            sparse = sparse,
                            max_depth = max_depth,
                            verbose = verbose)

      # df <- data.frame("goid"=sel$samples$goid,
      #                  "k" = values)
      if(nrow(sel$samples) == 1){
        df <- c(sel$samples$goid, values)
      }else{
        df <- cbind(sel$samples$goid, values)
      }
      return(df)
    }, future.packages = c("spNetwork"))
  }

  ## step 3  combining the results
  tot_df <- do.call(rbind,dfs)
  #tot_df <- tot_df[order(tot_df[,1]),]
  tot_df <- tot_df[match(events$goid,tot_df[,1]),]
  ## step 4 calculating the new bandwidth !
  if(ncol(tot_df) == 2){
    k <- tot_df[,2]
    delta <- calc_gamma(k)
    new_bw <- bw * (k**(-1/2) * delta**(-1))
    new_bw <- ifelse(new_bw<trim_bw,new_bw,trim_bw)
  }else{
    new_bw <- sapply(2:ncol(tot_df), function(i){
      k <- tot_df[,i]
      delta <- calc_gamma(k)
      nbw <- bw[[i-1]] * (k**(-1/2) * delta**(-1))
      nbw <- ifelse(nbw<trim_bw[[i-1]],nbw,trim_bw[[i-1]])
      return(nbw)
    })
  }
  return(new_bw)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### worker functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title NKDE worker
#'
#' @description The worker function for nkde and nkde.mc
#'
#' @param lines A feature collection of linestrings representing the network. The
#' geometries must be simple lines (may crash if some geometries
#'  are invalid)
#' @param events A feature collection of points representing the events on the
#' network. The points will be snapped on the network.
#' @param samples A feature collection of points representing the locations for
#' which the densities will be estimated.
#' @param kernel_name The name of the kernel to use
#' @param bw The global kernel bandwidth
#' @param bws The kernel bandwidth (in meters) for each event. Is usually a vector
#' but could also be a matrix if several global bandwidths were used. In this
#' case, the output value is also a matrix.
#' @param method The method to use when calculating the NKDE, must be one of
#' simple / discontinuous / continuous (see details for more information)
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#' events within the radius around each sampling point), "bw" (the bandwidth)
#' "none" (the simple sum).
#' @param digits The number of digits to keep in the spatial coordinates. It
#' ensures that topology is good when building the network. Default is 3
#' @param tol When adding the events and the sampling points to the network,
#' the minimum distance between these points and the lines extremities. When
#' points are closer, they are added at the extremity of the lines.
#' @param sparse A Boolean indicating if sparse or regular matrices should be
#' used by the Rcpp functions. Regular matrices are faster, but require more
#' memory and could lead to error, in particular with multiprocessing. Sparse
#' matrices are slower, but require much less memory.
#' @param max_depth When using the continuous and discontinuous methods, the
#' calculation time and memory use can go wild  if the network has a lot of
#' small edges (area with a lot of intersections and a lot of events). To
#' avoid it, it is possible to set here a maximum depth. Considering that the
#' kernel is divided at intersections, a value of 8 should yield good
#' estimates. A larger value can be used without problem for the discontinuous
#' method. For the continuous method, a larger value will strongly impact
#' calculation speed.
#' @param verbose A Boolean, indicating if the function should print messages
#' about the process.
#' @importFrom igraph adjacent_vertices get.edge.ids
#' @return A numeric vector with the nkde values
#' @keywords internal
#' @export
#' @examples
#' #This is an internal function, no example provided
nkde_worker <- function(lines, events, samples, kernel_name, bw, bws, method, div, digits, tol, sparse, max_depth, verbose = FALSE){

  # if we do not have event in that space, just return 0 values
  if(nrow(events)==0){
    values <- rep(0,nrow(samples))
    return(values)
  }

  ## step1 creating the graph
  if(verbose){
    print("    build graph ...")
  }
  graph_result <- build_graph(lines,digits = digits,line_weight = "length")

  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  edges <- graph_result$spedges

  ## step2 for each sample, find its belonging line
  a <- nrow(events)
  b <- nrow(lines)
  x <-  a*b
  if(is.na(x)){
    stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
  }
  if(x >= 2*10^9){
    stop(paste("The matrix size will be exceeded (",a," x ",b,"), please consider using a finer grid to split the study area",sep=""))
  }
  #snapped_samples <- maptools::snapPointsToLines(samples,edges,idField = "edge_id")
  snapped_samples <- snapPointsToLines2(samples,edges, snap_dist = max(bws), idField = "edge_id")
  samples$edge_id <- snapped_samples$nearest_line_id


  ## step3 finding for each event, its corresponding node
  events$vertex_id <- closest_points(events, nodes)

  ## step4 adding the spatial coordinates to samples and nodes
  XY_nodes <- st_coordinates(nodes)
  nodes$X_coords <- XY_nodes[,1]
  nodes$Y_coords <- XY_nodes[,2]

  XY_samples <- st_coordinates(samples)
  samples$X_coords <- XY_samples[,1]
  samples$Y_coords <- XY_samples[,2]

  ## step5 adding a local oid for samples and events
  events$oid <- 1:nrow(events)
  samples$oid <- 1:nrow(samples)

  ## step6 starting the calculations !

  if(verbose){
    print("        calculating NKDE values ...")
  }

  if(is.null(dim(bws))){
    # in this case, bws is a simple vector
    values <- calc_density(method, kernel_name, graph_result, graph, events, samples, bws, nodes, edges, div, max_depth, verbose, sparse)
  }else{
    values <- apply(bws, MARGIN = 2, FUN = function(x){
      vals <- calc_density(method, kernel_name, graph_result, graph, events, samples, x, nodes, edges, div, max_depth, verbose = FALSE, sparse)
      return(vals)
    })
  }
  return(values)
}


# A worker function to calculate densities. It has been added to facilitate bw selection in case of
# adaptive bandwidth
calc_density <- function(method, kernel_name, graph_result, graph, events, samples, bws, nodes, edges, div, max_depth, verbose, sparse){

  if(method == "simple"){
    # the cas of the simple method, no c++ here
    kernel_func <- select_kernel(kernel_name)
    if(verbose){
      values <- simple_nkde(graph, events, samples, bws, kernel_func, nodes, edges)
    }else{
      invisible(capture.output(values <- simple_nkde(graph, events, samples, bws, kernel_func, nodes, edges, div)))
    }

  }else{
    # we have to call the package with the Rcpp functions
    # this is necessary because this function can be used in a multicore context

    ## step6.5 preparing the complementary object for the rcpp function
    neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")
    neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})

    if(method=="continuous"){
      ##and finally calculating the values
      if (sparse){
        values <- spNetwork::continuous_nkde_cpp_arma_sparse(neighbour_list, events$vertex_id, events$weight,
                                                             st_drop_geometry(samples), bws, kernel_name, st_drop_geometry(nodes),
                                                             graph_result$linelist, max_depth, verbose, div)
      }else{
        values <- spNetwork::continuous_nkde_cpp_arma(neighbour_list, events$vertex_id, events$weight,
                                                      st_drop_geometry(samples), bws, kernel_name, st_drop_geometry(nodes),
                                                      graph_result$linelist, max_depth, verbose, div)
      }

    }

    if(method == "discontinuous"){
      #let this commented here to debug and test sessions
      # invisible(capture.output(values <- discontinuous_nkde2(edge_list,neighbour_list, events$vertex_id, events$weight,
      #                                                       samples@data, bw, kernel_func, nodes@data, graph_result$linelist, max_depth, verbose)))
      if(sparse){
        values <- spNetwork::discontinuous_nkde_cpp_arma_sparse(neighbour_list, events$vertex_id, events$weight,
                                                                st_drop_geometry(samples), bws, kernel_name, st_drop_geometry(nodes), graph_result$linelist, max_depth, verbose)
      }else{
        values <- spNetwork::discontinuous_nkde_cpp_arma(neighbour_list, events$vertex_id, events$weight,
                                                         st_drop_geometry(samples), bws, kernel_name, st_drop_geometry(nodes), graph_result$linelist, max_depth, verbose)
      }


    }

  }

  ## step7 adjusting the kernel values !
  if(div == "n"){
    return(values$sum_k / values$n)
  }else if (div == "bw"){
    #return(values$sum_k * (1/bw))
    # NOTE : if div == "bw", then the scaling is done during the estimation
    # because we might have adaptive BW
    return(values$sum_k)
  }else if (div == "none"){
    return(values$sum_k)
  }
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### main functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Network Kernel density estimate
#'
#' @description Calculate the Network Kernel Density Estimate based on a network of lines,
#' sampling points, and events
#' @md
#' @details
#' **The three NKDE methods**\cr
#' Estimating the density of a point process is commonly done by using an
#' ordinary two-dimensional kernel density function. However, there are numerous
#' cases for which the events do not occur in a two-dimensional space but on a
#' network (like car crashes, outdoor crimes, leaks in pipelines, etc.). New
#' methods were developed to adapt the methodology to networks, three of them
#' are available in this package.
#'
#' * The simple method: This first method was presented by \insertCite{xie2008kernel}{spNetwork} and proposes an
#' intuitive solution. The distances between events and sampling points are
#' replaced by network distances, and the formula of the kernel is adapted to
#' calculate the density over a linear unit instead of an areal unit.
#' * The discontinuous method: The previous method has been criticized by
#' \insertCite{okabe2009kernel}{spNetwork}, arguing that the estimator proposed
#' is biased, leading to an overestimation of density in events hot-spots. More
#' specifically, the simple method does not conserve mass and the induced kernel
#' is not a probability density along the network. They thus proposed a
#' discontinuous version of the kernel function on network, which equally
#' "divides" the mass density of an event at intersections.
#' * The continuous method: If the discontinuous method is unbiased, it leads
#' to a discontinuous kernel function which is a bit counter-intuitive.
#' \insertCite{okabe2009kernel;textual}{spNetwork} proposed another version of
#' the kernel, which divides the mass of the density at intersections but adjusts
#' the density before the intersection to make the function continuous.
#'
#' The three methods are available because, even though that the simple method is
#' less precise statistically speaking, it might be more intuitive. From a
#' purely geographical view, it might be seen as a sort of distance decay
#' function as used in Geographically Weighted Regression.\cr \cr\cr **adaptive bandwidth**\cr
#' It is possible to use adaptive bandwidth instead of fixed
#' bandwidth. Adaptive bandwidths are calculated using the Abramson’s smoothing
#' regimen \insertCite{abramson1982bandwidth}{spNetwork}. To do so, an original
#' fixed bandwidth must be specified (bw parameter), and is used to estimate the
#' priory densitiy at event locations. These densities are then used to
#' calculate local bandwidth. The maximum size of the local bandwidth can be
#' limited with the parameter trim_bw. For more details, see the vignettes.
#' \cr\cr **Optimization parameters**\cr The grid_shape parameter allows to
#' split the calculus of the NKDE according to a grid dividing the study area.
#' It might be necessary for big dataset to reduce the memory used. If the
#' grid_shape is c(1,1), then a full network is built for the area. If the
#' grid_shape is c(2,2), then the area is split in 4 rectangles. For each
#' rectangle, the sample points falling in the rectangle are used, the events
#' and the lines in a radius of the bandwidth length are used. The results are
#' combined at the end and ordered to match the original order of the samples.
#' \cr\cr The geographical coordinates of the start and end of lines are used to
#' build the network. To avoid troubles with digits, we truncate the coordinates
#' according to the digit parameter. A minimal loss of precision is expected but
#' results in a fast construction of the network. \cr\cr To calculate the
#' distances on the network, all the events are added as vertices. To reduce the
#' size of the network, it is possible to reduce the number of vertices by
#' adding the events at the extremity of the lines if they are close to them.
#' This is controlled by the parameter tol. \cr\cr In the same way, it is
#' possible to limit the number of vertices by aggregating the events that are
#' close to each other. In that case, the weights of the aggregated events are
#' summed. According to an aggregation distance, a buffer is drawn around the
#' fist event, all events falling in that buffer are aggregated to the first
#' event, forming a new event. The coordinates of this new event are the means of
#' the original events coordinates. This procedure is repeated until no events
#' are aggregated. The aggregation distance can be fixed with the parameter agg.
#' \cr\cr When using the continuous and discontinuous kernel, the density is
#' reduced at each intersection crossed. In the discontinuous case, after 5
#' intersections with four directions each, the density value is divided by 243
#' leading to very small values. In the same situation but with the continuous
#' NKDE, the density value is divided by approximately 7.6. The max_depth
#' parameters allows the user to control the maximum depth of these two NKDE.
#' The base value is 15, but a value of 10 would yield very close estimates. A
#' lower value might have a critical impact on speed when the bandwidth is large.
#' \cr\cr When using the continuous and discontinuous kernel, the connections
#' between graph nodes are stored in a matrix. This matrix is typically sparse,
#' and so a sparse matrix object is used to limit memory use. If the network is
#' small (typically when the grid used to split the data has small rectangles)
#' then a classical matrix could be used instead of a sparse one. It
#' significantly increases speed, but could lead to memory issues.
#'
#' @references{ \insertAllCited{} }
#'
#' @template nkde_params-arg
#' @template diggle_corr-arg
#' @param samples A feature collection of points representing the locations for
#'   which the densities will be estimated.
#' @template nkde_params2-arg
#' @template nkde_geoms-args
#' @template sparse-arg
#' @template grid_shape-arg
#' @template verbose-arg
#' @template check-arg
#' @return A vector of values, they are the density estimates at sampling points
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' data(bike_accidents)
#' lixels <- lixelize_lines(mtl_network,200,mindist = 50)
#' samples <- lines_center(lixels)
#' densities <- nkde(mtl_network,
#'                   events = bike_accidents,
#'                   w = rep(1,nrow(bike_accidents)),
#'                   samples = samples,
#'                   kernel_name = "quartic",
#'                   bw = 300, div= "bw",
#'                   adaptive = FALSE,
#'                   method = "discontinuous", digits = 1, tol = 1,
#'                   agg = 15,
#'                   grid_shape = c(1,1),
#'                   verbose=FALSE)
#' }
nkde <- function(lines, events, w, samples, kernel_name, bw, adaptive=FALSE, trim_bw=NULL, method, div="bw", diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE, grid_shape=c(1,1), verbose=TRUE, check=TRUE){

  ## step0 basic checks REPRENDRE ICI
  if(verbose){
    print("checking inputs ...")
  }
  if((kernel_name %in% c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov','uniform'))==FALSE){
    stop('the name of the kernel function must be one of c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov" ,"uniform")')
  }

  if(method %in% c("simple","continuous","discontinuous") == FALSE){
    stop('the method must be one of c("simple","continuous","discontinuous"')
  }
  if(method == "continuous" & kernel_name == "gaussian"){
    stop("using the continuous NKDE and the gaussian kernel function can yield negative values for densities because the gaussian kernel does not integrate to 1 within the bandiwdth, please consider using the quartic kernel instead")
  }

  if(min(bw)<=0){
    stop("the bandwidth for the kernel must be superior to 0")
  }

  if(adaptive & (length(bw) > 1 )){
    stop("When adaptive is TRUE, only a global bandwidth must be given")
  }

  if(adaptive & is.null(trim_bw)){
    stop("if adaptive is TRUE, a value for trim_bw must be supplied")
  }
  if(diggle_correction & is.null(study_area)){
    stop("the study_area must be defined if the Diggle correction factor is used")
  }
  if(check){
    checked <- check_geometries(lines,samples,events, study_area)
  }


  ## step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
  }
  events$bws <- bw
  data <- prepare_data(samples, lines, events, w, digits, tol, agg)
  lines <- data$lines
  samples <- data$samples
  events <- data$events

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))
  grid <- sf::st_as_sf(grid)

  ## adaptive bandwidth !
  if(adaptive == FALSE){
    # if (length(bw) == 1){
    #   bws <- rep(bw,nrow(events))
    # }else{
    #   bws <- bw
    # }
    bws <- events$bws

  }else{
    ## we want to use an adaptive bw
    # bws <- adaptive_bw.mc(grid, events, lines, bw, trim_bw, method,
    #                       kernel_name, max_depth, tol, digits, sparse, verbose)
    bws <- adaptive_bw(grid, events, lines, bw, trim_bw, method,
                          kernel_name, max_depth, tol, digits, sparse, verbose)
  }

  ## calculating the correction factor
  if(diggle_correction){
    if(verbose){
      print("Calculating the correction factor")
    }
    corr_factor <- correction_factor(study_area,events,lines,method,bws, kernel_name, tol, digits, max_depth, sparse)
  }else{
    corr_factor <- rep(1,nrow(events))
  }
  events$weight <- events$weight * corr_factor

  events$bw <- bws
  max_bw <- max(bws)

  ## step3 splitting the dataset with each rectangle
  if(verbose){
    print("Splitting the data with the spatial grid ...")
  }
  selections <- split_by_grid(grid,samples,events,lines,max_bw, tol, digits)

  ## step 4 calculating the values
  if(verbose){
    print("start calculating the kernel values ...")
  }
  n_quadra <- length(selections)

  dfs <- lapply(1:n_quadra,function(i){

    sel <- selections[[i]]

    if(verbose){
      print(paste("    quadra ",i,"/",n_quadra,sep=""))
    }

    values <- nkde_worker(lines = sel$lines,
                          events = sel$events,
                          samples =  sel$samples,
                          kernel_name = kernel_name,
                          bw = bw,
                          bws = sel$events$bw,
                          method = method,
                          div = div,
                          digits = digits,
                          tol = tol,
                          sparse = sparse,
                          max_depth = max_depth,
                          verbose = verbose)

    df <- data.frame("goid"=sel$samples$goid,
                     "k" = values)
    return(df)
  })

  if(verbose){
    print("combining the results ...")
  }

  ## step5  combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]
  if(adaptive){
    return(list("events" = events,
                "k" = tot_df$k))
  }else{
    return(tot_df$k)
  }
}


#' @title Network Kernel density estimate (multicore)
#'
#' @description Calculate the Network Kernel Density Estimate based on a network of lines,
#' sampling points, and events with multicore support.
#'
#' @details For more details, see help(nkde)
#'
#' @template nkde_params-arg
#' @template diggle_corr-arg
#' @param samples A feature collection of points representing the locations for
#' which the densities will be estimated.
#' @template nkde_params2-arg
#' @template nkde_geoms-args
#' @template sparse-arg
#' @template grid_shape-arg
#' @param verbose A Boolean, indicating if the function should print messages
#' about the process.
#' @template check-arg
#' @return A vector of values, they are the density estimates at sampling
#' points
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' data(bike_accidents)
#' future::plan(future::multisession(workers=1))
#' lixels <- lixelize_lines(mtl_network,200,mindist = 50)
#' samples <- lines_center(lixels)
#' densities <- nkde.mc(mtl_network,
#'                   events = bike_accidents,
#'                   w = rep(1,nrow(bike_accidents)),
#'                   samples = samples,
#'                   kernel_name = "quartic",
#'                   bw = 300, div= "bw",
#'                   adaptive = FALSE, agg = 15,
#'                   method = "discontinuous", digits = 1, tol = 1,
#'                   grid_shape = c(3,3),
#'                   verbose=TRUE)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
nkde.mc <- function(lines, events, w, samples, kernel_name, bw, adaptive=FALSE, trim_bw=NULL, method, div="bw", diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1,agg=NULL, sparse=TRUE, grid_shape=c(1,1), verbose=TRUE, check=TRUE){

  ## step0 basic checks
  if(verbose){
    print("checking inputs ...")
  }

  if((kernel_name %in% c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov','uniform'))==FALSE){
    stop('the name of the kernel function must be one of c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov" ,"uniform")')
  }

  if(method %in% c("simple","continuous","discontinuous") == FALSE){
    stop('The method must be one of c("simple","continuous","discontinuous"')
  }
  if(method == "continuous" & kernel_name == "gaussian"){
    stop("using the continuous NKDE and the gaussian kernel function can yield negative values for densities because the gaussian kernel does not integrate to 1 within the bandiwdth, please consider using the quartic kernel instead")
  }
  if(min(bw)<=0){
    stop("the bandwidth for the kernel must be superior to 0")
  }
  if(adaptive & (length(bw) > 1 )){
    stop("When adaptive is TRUE, only a global bandwidth must be given")
  }
  if(adaptive & is.null(trim_bw)){
    stop("if adaptive is TRUE, a value for trim_bw must be supplied")
  }
  if(diggle_correction & is.null(study_area)){
    stop("the study_area must be defined if the Diggle correction factor is used")
  }

  if(check){
    check_geometries(lines,samples,events,study_area)
  }

  ## step1 preparing the data
  if(verbose){
    print("prior data preparation ...")
  }

  events$bws <- bw
  data <- prepare_data(samples, lines, events,w,digits,tol,agg)
  lines <- data$lines
  samples <- data$samples
  events <- data$events

  ## step2 creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))

  ## adaptive bandwidth !
  if(adaptive==FALSE){
    bws <- events$bws
  }else{
    if(verbose){
      print("calculating the local bandwidth ...")
    }
    ## we want to use an adaptive bw
    bws <- adaptive_bw.mc(grid,events,lines,bw,trim_bw,method,kernel_name,max_depth,tol,digits,sparse,verbose)
  }

  ## calculating the correction factor
  if(diggle_correction){
    if(verbose){
      print("Calculating the correction factor")
    }
    corr_factor <- correction_factor(study_area,events,lines,method,bws, kernel_name, tol, digits, max_depth, sparse)
  }else{
    corr_factor <- rep(1,nrow(events))
  }
  events$weight <- events$weight * corr_factor

  events$bw <- bws
  max_bw <- max(bws)

  ## step3 splitting the dataset with each rectangle
  selections <- split_by_grid(grid,samples,events,lines,max_bw, digits,tol)

  ## step4 calculating the values

  if(verbose){
    print("starting the kernel calculations ...")
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {

        invisible(capture.output(values <- nkde_worker(sel$lines, sel$events,
                              sel$samples, kernel_name,bw,
                              sel$events$bw, method, div, digits,
                              tol,sparse, max_depth, verbose)))

        df <- data.frame("goid"=sel$samples$goid,
                         "k" = values)
        p(sprintf("i=%g", sel$index))
        return(df)
      }, future.packages = c("spNetwork"))
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {

      values <- nkde_worker(sel$lines, sel$events,
                            sel$samples, kernel_name,bw,
                            sel$events$bw, method, div, digits,
                            tol,sparse, max_depth, verbose)

      df <- data.frame("goid"=sel$samples$goid,
                       "k" = values)
      return(df)
    }, future.packages = c("spNetwork"))
  }


  if(verbose){
    print("combining the results ...")
  }
  ## step5 combining the results
  tot_df <- do.call(rbind,dfs)
  tot_df <- tot_df[order(tot_df$goid),]
  events$bws <- NULL
  if(adaptive){
    return(list("events" = events,
                "k" = tot_df$k))
  }else{
    return(tot_df$k)
  }
}


