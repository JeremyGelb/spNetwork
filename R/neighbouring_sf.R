#' @title Select the distance to weight function
#'
#' @description Select a function to convert distance to weights
#' if a function is provided, this function will be vectorized.
#'
#' @param dist_func Could be a name in c('inverse', 'identity',
#' 'squared inverse') or a function with only one parameter x
#' @return A vectorized function used to convert distance into spatial weights
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
select_dist_function <- function(dist_func = "inverse") {
    if (class(dist_func) == "character") {
        if (dist_func == "identity") {
            dist_func <- function(x) {
                return(x)
            }
        } else if (dist_func == "inverse") {
            dist_func <- function(x) {
                return(1 / x)
            }
        } else if (dist_func == "squared inverse") {
            dist_func <- function(x) {
                return(1 / x^2)
            }
        }
    }
    vdist_func <- Vectorize(dist_func, "x")
    return(vdist_func)
}

#defining some global variables (weird felx but ok)
utils::globalVariables(c("origin", "fid"))

#' @title network_listw worker
#'
#' @description The worker function of network_listw.
#'
#' @param points A feature collection of points corresponding to start and end
#' points. It must have a column fid, grouping the points if necessary.
#' @param lines A feature collection of lines representing the network
#' @param maxdistance The maximum distance between two observation to
#' consider them as neighbours.
#' @param dist_func A vectorized function converting spatial distances into
#' weights.
#' @param mindist The minimum distance between two different observations.
#' It is important for it to be different from 0 when a W style is used.
#' @param mindist The minimum distance between two different observations.
#' It is important for it to be different from 0 when a W style is used.
#' @param direction Indicate a field giving information about authorized
#' travelling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details
#' @param verbose A Boolean indicating if the function should print its
#' progress
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A list of neihbours as weights.
#' @importFrom sf st_drop_geometry
#' @importFrom utils capture.output
#' @importFrom data.table data.table .SD transpose
#' @keywords internal
#' @examples
#' #no example provided, this is an internal function
network_listw_worker<-function(points,lines,maxdistance,dist_func, direction=NULL, mindist=10, matrice_type = "B", verbose = FALSE, digits = 3, tol=0.1){
    ## step1 : adding the points to the lines
    lines$worker_id <- 1:nrow(lines)
    points$nearest_line_id <- as.numeric(as.character(points$nearest_line_id))
    joined <- data.table(st_drop_geometry(points))
    B <- data.table(st_drop_geometry(lines))
    joined[B,  on = c("nearest_line_id" = "tmpid"),
           names(B) := mget(paste0("i.", names(B)))]

    if(verbose){
        print("adding the points as vertices to nearest lines")
    }
    # in the previous version, we split all the lines at their vertex
    # however, this is not required here, we could just split them
    # at points added in the network
    # new_lines <- add_vertices_lines(lines, points, joined$worker_id, 1)
    graph_lines <- split_lines_at_vertex(lines, points, joined$worker_id, 1)

    ## step2 : splitting the lines on vertices and adjusting weights
    if(verbose){
        print("splitting the lines for the network")
    }
    #graph_lines <- simple_lines(new_lines)
    graph_lines$lx_length <- as.numeric(sf::st_length(graph_lines))
    graph_lines$lx_weight <- (graph_lines$lx_length / graph_lines$line_length) * graph_lines$line_weight

    ## step3 building the network
    if(verbose){
        print("generating the network")
    }
    if (is.null(direction)){
        result_graph <- build_graph(graph_lines, digits = digits,
                                    attrs = TRUE, line_weight = "lx_weight")
    }else{
        #dir <- ifelse(graph_lines[[direction]]=="Both",0,1)
        result_graph <- build_graph_directed(graph_lines, digits = digits,
                                    attrs = TRUE, line_weight='lx_weight',
                                    direction = direction)
    }
    ## step4 finding for each point its corresponding vertex
    points$vertex <- closest_points(points,result_graph$spvertices)
    if(verbose){
        print("calculating the distances on the network")
    }
    starts <- subset(points,points$pttype=="start")
    u <- unique(points$vertex)

    ## step 4 calculating the distances between the start and end points
    base_distances <- igraph::distances(result_graph$graph,
                                        v = starts$vertex, to = u,
                                        mode = "out")
    all_ditancesdt <- data.table(base_distances)
    all_ditancesdt$origin <- starts$fid

    ## step 5 extracting the distances as a dataframe rather than a matrix

    #first groupby on rows :
    #multiple rows might correspond to the same object (if multiple points like)
    step1 <- all_ditancesdt[, lapply(.SD, min, na.rm = TRUE), by = origin ]
    originid <- step1$origin

    #then transpose and merge
    step2 <- transpose(step1[,2:ncol(step1)])
    step2$vertex <- u
    pts <- data.table(st_drop_geometry(points[c("fid","vertex")]))
    step2.5 <- data.table::merge.data.table(step2,pts,by="vertex",all=TRUE)
    colorder <- c(2:(ncol(step2.5)-1),1,ncol(step2.5))
    step2.5 <- data.table::setcolorder(step2.5,colorder)

    step3 <- step2.5[, lapply(.SD, min, na.rm=TRUE), by=fid ]
    destid <- step3$fid
    #transpose again
    step4 <- transpose(step3)
    dist_matrix <- as.matrix(step4[2:(nrow(step4)-1)])
    dist_matrix <- ifelse(is.na(dist_matrix),Inf,dist_matrix)
    rownames(dist_matrix) <- originid
    colnames(dist_matrix) <- destid

    ##step 6 finaly generate the neighbouring object with the data
    nblist <- lapply(1:nrow(dist_matrix),function(i){
        row <- dist_matrix[i,]
        iid <- originid[i]
        test <- row<=maxdistance & destid!=iid
        if(any(test)){
            nbs <- destid[test]
            okdistances <- row[test]
            okdistances <- ifelse(okdistances==0,mindist,okdistances)
            weights <- dist_func(okdistances)
            if(matrice_type=="B"){
                wtdweights <- rep(1,length(nbs))
            }else if (matrice_type=="W"){
                wtdweights <- weights/sum(weights)
            }
            return(list(nbs,as.numeric(wtdweights)))
        }else{
            return(list(0L, NULL))
        }
    })
    nbs <- lapply(1:length(nblist), function(j) {
        nblist[[j]][[1]]
    })
    names(nbs) <- originid
    weights <- lapply(1:length(nblist), function(j) {
        nblist[[j]][[2]]
    })
    names(weights) <- originid
    return(list(nbs, weights))
}



#' @title Data preparation for network_listw
#'
#' @description Function to prepare selected points and selected lines during
#'   the process.
#'
#' @param is The indices of the quadras to use in the grid
#' @param grid A feature collection of polygons representing the quadras to split
#'   calculus
#' @param snapped_points The start and end points snapped to the lines
#' @param lines The lines representing the network
#' @param maxdistance The maximum distance between two observation to considere
#'   them as neighbours.
#' @return A list of two elements : selected points and selected lines
#' @keywords internal
#' @importFrom sf st_bbox st_buffer
#' @examples
#' #no example provided, this is an internal function
prepare_elements_netlistw <- function(is,grid,snapped_points,lines,maxdistance){

    ## step1 preparing the spatial indices
    lines_tree <- spNetwork::build_quadtree(lines)
    snapped_points_tree <-  spNetwork::build_quadtree(snapped_points)

    snapped_points$oids <- 1:nrow(snapped_points)
    ## step 2 extracting the quadra
    results <- lapply(is,function(i){
        quadra <- grid[i,]
        #selecting the starting points
        start_pts <-  spNetwork::spatial_request(quadra,snapped_points_tree, snapped_points)
        if(nrow(start_pts)==0){
            return(NULL)
        }else{
            #start_pts <- snapped_points[snapped_points$fid %in% start_pts$fid,]
            start_pts$pttype <- "start"
            #selecting the endpoints
            poly <- st_bbox_geom(start_pts)
            buff <- st_buffer(poly, dist = maxdistance)
            if(nrow(start_pts)==nrow(snapped_points)){
                all_pts <- start_pts
            }else{
                end_pts <-  spNetwork::spatial_request(buff, snapped_points_tree, snapped_points)
                end_pts <- subset(end_pts,(end_pts$oids %in% start_pts$oids)== FALSE)
                end_pts$pttype <- "end"
                #combining all the points
                all_pts <- rbind(start_pts,end_pts)
            }
            #selecting the lines
            selected_lines <-  spNetwork::spatial_request(buff, lines_tree, lines)
            #calculating the elements
            return(list(all_pts,selected_lines))
        }
    })
    return(results)
}



#' @title Network distance listw
#'
#' @description Generate listw object (spdep like) based on network distances.
#'
#' @param origins A feature collection of lines, points, or polygons
#' for which the spatial neighbouring list will be built
#' @param lines A feature collection of lines representing the network
#' @param maxdistance The maximum distance between two observations to
#' consider them as neighbours.
#' @param method A string indicating how the starting points will be built.
#' If 'centroid' is used, then the centre of lines or polygons is used. If
#' 'pointsalong' is used, then points will be placed along polygons' borders or
#' along lines as starting and end points. If 'ends' is used (only for lines)
#' the first and last vertices of lines are used as starting and ending points.
#' @param point_dist A float, defining the distance between points when the
#' method 'pointsalong' is selected.
#' @param snap_dist The maximum distance to snap the start and end points on
#' the network.
#' @param line_weight The weighting to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional to the geographical length of the lines.
#' @param mindist The minimum distance between two different observations.
#' It is important for it to be different from 0 when a W style is used.
#' @param direction Indicates a field providing information about authorized
#' travelling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @param dist_func Indicates the function to use to convert the distance
#' between observation in spatial weights. Can be 'identity', 'inverse',
#' 'squared inverse' or a function with one parameter x that will be
#' vectorized internally
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(1,1), so all the calculation is
#' done in one go. It might be necessary to split it if the dataset is large.
#' @param verbose A Boolean indicating if the function should print its
#' progress
#' @param digits The number of digits to retain in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A listw object (spdep like)
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom graphics plot
#' @importFrom sf st_point_on_surface st_length
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg",package = "spNetwork", mustWork = TRUE)
#' mtl_network <- sf::st_read(networkgpkg,layer="mtl_network")
#' listw <- network_listw(mtl_network, mtl_network, maxdistance = 500,
#'         method = "centroid", line_weight = "length",
#'         dist_func = 'squared inverse', matrice_type='B', grid_shape = c(2,2))
#' }
network_listw <- function(origins,lines, maxdistance, method="centroid", point_dist=NULL, snap_dist=Inf, line_weight = "length", mindist=10, direction=NULL, dist_func = "inverse", matrice_type = "B", grid_shape=c(1,1), verbose = FALSE, digits = 3, tol=0.1){

    ## step1 adjusting the weights of the lines
    lines$line_length <- as.numeric(sf::st_length(lines))
    if(line_weight=="length"){
        lines$line_weight <- as.numeric(sf::st_length(lines))
    }else {
        lines$line_weight <- lines[[line_weight]]
    }
    if(min(lines$line_weight)<=0){
        stop("the weights of the lines must be superior to 0")
    }

    ## step2 adjusting the directions of the lines (done now by the graph function)
    #if(is.null(direction) == FALSE){
    #    lines <- lines_direction(lines,direction)
    #}

    ## step3  checking the matrix type
    if (matrice_type %in% c("B", "W") == FALSE) {
        stop("Matrice type must be B (binary) or W (row standardized)")
    }
    ## step 4 creating the vectorized distance function
    vdist_func <- select_dist_function(dist_func)
    lines$tmpid <- 1:nrow(lines)

    ## step5 generating the starting points
    if(verbose){
        print("generating the starting points")
    }
    origins$fid <- 1:nrow(origins)
    if(unique(st_geometry_type(origins))=="POLYGON"){
        if(method=="centroid"){
            centers_geom <- st_point_on_surface(origins$geometry)
            centers <- origins
            centers$geometry <- centers_geom
        }else if(method=="pointsalong"){
            centers <- surrounding_points(origins,point_dist)
        }
    }else if(unique(st_geometry_type(origins))=="POINT"){
        centers <- origins
    }else if(unique(st_geometry_type(origins))=="LINESTRING"){
        if(method=="centroid"){
            if(verbose){
                print("getting the centers of the lines ...")
                centers <- lines_center(origins)
            }else{
                invisible(capture.output(centers <- lines_center(origins)))
            }
        }else if(method=="ends"){
            centers <- lines_extremities(origins)
        }else if(method=="pointsalong"){
            centers <- lines_points_along(origins,point_dist)
        }else{
            stop("if origins class is feature collection of linestrings, method must be on of centroid, ends or pointsalong")
        }
    }else{
        stop("origins must be a feature collection of points, polygons or linestrings")
    }

    if(verbose){
        print("snapping the points to the lines (only once)")
    }
    #snapped_points <- maptools::snapPointsToLines(centers,lines,maxDist = snap_dist, idField="tmpid")
    snapped_points <- snapPointsToLines2(centers,lines, idField="tmpid")
    #snapped_points <- cbind(snapped_points, centers)

    ## step 6 building grid
    mygrid <- build_grid(grid_shape,list(origins,lines))
    if(verbose){
        print("starting the network part")
    }
    ids <- 1:nrow(mygrid)
    list_elements <- prepare_elements_netlistw(ids,mygrid,snapped_points,lines,maxdistance)

    ## step7 iterating over the grid
    listvalues <- lapply(1:nrow(mygrid),function(i){
        if(verbose){
            print(paste("working on quadra : ",i,"/",max(ids),sep=""))
        }
        elements <- list_elements[[i]]

        if(length(elements)==0){
            return()
        }else {
            all_pts <- elements[[1]]
            selected_lines <- elements[[2]]
            #calculating the elements
            values <- network_listw_worker(all_pts, selected_lines,
                                           maxdistance, direction=direction, mindist=mindist,
                                           dist_func = vdist_func, matrice_type = matrice_type,
                                           verbose = verbose, digits = digits, tol=tol)
            return(values)
        }

    })

    if(verbose){
        print("building the final listw object")
    }

    ## step8 combining the results
    okvalues <- listvalues[lengths(listvalues) != 0]
    nblist <- do.call("c", lapply(okvalues, function(value) {
        return(value[[1]])
    }))
    ordered_nblist <- lapply(as.character(origins$fid), function(x) {
        return(nblist[[x]])
    })
    nbweights <- do.call("c", lapply(okvalues, function(value) {
        return(value[[2]])
    }))
    names(nbweights) <- names(nblist)
    ordered_weights <- lapply(as.character(origins$fid), function(x) {
        return(nbweights[[x]])
    })
    # setting the nb attributes
    attr(ordered_nblist, "class") <- "nb"
    attr(ordered_nblist, "region.id") <- as.character(origins$fid)
    if(verbose){
        print("finally generating the listw object ...")
    }
    # setting the final listw attributes
    listw <- spdep::nb2listw(ordered_nblist, glist = ordered_weights,
                             zero.policy = TRUE, style = matrice_type)
    return(listw)

}

#' @title Network distance listw (multicore)
#'
#' @description Generate listw object (spdep like) based on network distances with multicore support.
#'
#' @param origins A feature collection of linestrings, points or polygons for
#' which the spatial neighbouring list will be built.
#' @param lines A feature collection of linestrings representing the network
#' @param maxdistance The maximum distance between two observations to
#' consider them as neighbours.
#' @param method A string indicating how the starting points will be built.
#' If 'centroid' is used, then the centre of lines or polygons is used. If
#' 'pointsalong' is used, then points will be placed along polygons' borders or
#' along lines as starting and end points. If 'ends' is used (only for lines)
#' the first and last vertices of lines are used as starting and ending points.
#' @param point_dist A float, defining the distance between points when the
#' method pointsalong is selected.
#' @param snap_dist the maximum distance to snap the start and end points on
#' the network.
#' @param line_weight The weights to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional with the geographical length of the lines.
#' @param mindist The minimum distance between two different observations.
#' It is important for it to be different from 0 when a W style is used.
#' @param direction Indicates a field giving information about authorized
#' travelling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @param dist_func Indicates the function to use to convert the distance
#' between observation in spatial weights. Can be 'identity', 'inverse',
#' 'squared inverse' or a function with one parameter x that will be
#' vectorized internally
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(1,1), so all the calculation is
#' done in one go. It might be necessary to split it if the dataset is large.
#' @param verbose A Boolean indicating if the function should print its
#' progress
#' @param digits The number of digits to retain in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A listw object (spdep like)
#' @importFrom sf st_length st_point_on_surface
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- sf::st_read(networkgpkg,layer="mtl_network")
#' future::plan(future::multisession(workers=2))
#' listw <- network_listw.mc(mtl_network,mtl_network,maxdistance=500,
#'         method = "centroid", line_weight = "length",
#'         dist_func = 'squared inverse', matrice_type='B', grid_shape = c(2,2))
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'}
network_listw.mc <- function(origins,lines,maxdistance, method="centroid", point_dist=NULL, snap_dist=Inf, line_weight = "length", mindist=10, direction=NULL, dist_func = "inverse", matrice_type = "B", grid_shape=c(1,1), verbose = FALSE, digits = 3, tol=0.1){
    ##adjusting the weights of the lines
    lines$line_length <- as.numeric(st_length(lines))
    if(line_weight=="length"){
        lines$line_weight <- as.numeric(st_length(lines))
    }else {
        lines$line_weight <- lines[[line_weight]]
    }
    if(min(lines$line_weight)<=0){
        stop("the weights of the lines must be superior to 0")
    }

    # adjusting the directions of the lines
    if(is.null(direction) == FALSE){
        lines <- lines_direction(lines,direction)
    }

    ## checking the matrix type
    if (matrice_type %in% c("B", "W") == FALSE) {
        stop("Matrice type must be B (binary) or W (row standardized)")
    }
    ## creating the vectorized distance function
    vdist_func <- select_dist_function(dist_func)
    ## setting the OID : tmpid
    lines$tmpid <- 1:nrow(lines)

    ##generating the starting points
    ## step5 generating the starting points
    if(verbose){
      print("generating the starting points")
    }
    origins$fid <- 1:nrow(origins)
    if(unique(st_geometry_type(origins))=="POLYGON"){
      if(method=="centroid"){
        centers_geom <- st_point_on_surface(origins$geometry)
        centers <- origins
        centers$geometry <- centers_geom
      }else if(method=="pointsalong"){
        centers <- surrounding_points(origins,point_dist)
      }
    }else if(unique(st_geometry_type(origins))=="POINT"){
      centers <- origins
    }else if(unique(st_geometry_type(origins))=="LINESTRING"){
      if(method=="centroid"){
        if(verbose){
          print("getting the centers of the lines ...")
          centers <- lines_center(origins)
        }else{
          invisible(capture.output(centers <- lines_center(origins)))
        }
      }else if(method=="ends"){
        centers <- lines_extremities(origins)
      }else if(method=="pointsalong"){
        centers <- lines_points_along(origins,point_dist)
      }else{
        stop("if origins class is feature collection of linestrings, method must be on of centroid, ends or pointsalong")
      }
    }else{
      stop("origins must be a feature collection of points, polygons or linestrings")
    }
    if(verbose){
        print("snapping the points to the lines (only once)")
    }
    #snapped_points <- maptools::snapPointsToLines(centers,lines,maxDist = snap_dist, idField="tmpid")
    snapped_points <- snapPointsToLines2(centers,lines, idField="tmpid")
    #snapped_points <- cbind(snapped_points, centers)

    ##building grid
    grid <- build_grid(grid_shape,list(origins,lines))

    if(verbose){
        print("starting the network part")
    }

    all_is <- 1:nrow(grid)
    iseq <- list()
    cnt <- 0
    for(i in 1:grid_shape[[1]]){
        start <- cnt*grid_shape[[2]]+1
        iseq[[length(iseq)+1]] <- list(cnt+1,all_is[start:(start+grid_shape[[2]]-1)])
        cnt<-cnt+1
    }
    listelements <- future.apply::future_lapply(iseq,function(ii){
        elements <- prepare_elements_netlistw(ii[[2]],grid,snapped_points,lines,maxdistance)
        return(elements)
    }, future.packages = c("sf", "spNetwork"))
    listelements <- unlist(listelements,recursive = FALSE)
    ##iterating on the grid
    listvalues <- future.apply::future_lapply(listelements,function(elements){
        ##step1 : preparing elements
        if(is.null(elements)){
            return()
        }else {
            all_pts <- elements[[1]]
            selected_lines <- elements[[2]]
            #calculating the elements
            values <- network_listw_worker(all_pts, selected_lines,
                                           maxdistance, direction=direction, mindist=mindist,
                                           dist_func = vdist_func, matrice_type = matrice_type,
                                           verbose = verbose, digits = digits, tol=tol)
            return(values)
        }

    }, future.packages = c("sf", "spNetwork"))

    if(verbose){
        print("building the final listw object")
    }

    okvalues <- listvalues[lengths(listvalues) != 0]
    nblist <- do.call("c", lapply(okvalues, function(value) {
        return(value[[1]])
    }))
    ordered_nblist <- lapply(as.character(origins$fid), function(x) {
        return(nblist[[x]])
    })
    nbweights <- do.call("c", lapply(okvalues, function(value) {
        return(value[[2]])
    }))
    names(nbweights) <- names(nblist)
    ordered_weights <- lapply(as.character(origins$fid), function(x) {
        return(nbweights[[x]])
    })
    ## setting the nb attributes
    attr(ordered_nblist, "class") <- "nb"
    attr(ordered_nblist, "region.id") <- as.character(origins$fid)
    if(verbose){
        print("finally generating the listw object ...")
    }
    ## setting the final listw attributes
    listw <- spdep::nb2listw(ordered_nblist, glist = ordered_weights,
                             zero.policy = TRUE, style = matrice_type)
    return(listw)

}
