#' A helper function to select a function to convert distance to weights
#'if a function is provided, this function will be vectorized.
#'
#'@param dist_func COuld be a name in c('inverse', 'identity',
#''squared inverse') or a function with only one parameter x
#'@return a vectorized function used to convert distance into spatial weights
#'@examples
#'#This is an internal function, no example provided
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


#' The worker function of network_listw_gridded
#'
#' @param points A SpatialPointsDataFrame corresponding to start and end
#' points. It must have a column fid, grouping the points if necessary.
#' @param lines A SpatialLinesDataFrame representing the network
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @param dist_func A vectorized function converting spatial distances into
#' weights.
#' @param mindist The minimum distance between two different observations.
#' It is important for it to be different from 0 when a W style is used.
#' @param mindist The minimum distance between two different observations.
#' It is important for it to be different from 0 when a W style is used.
#' @param direction indicate a field giving informations about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details
#' @param verbose A string indicating how the advance of the process is
#' displayed. Default is all. "text" show only the main steps, "progressbar"
#' show main steps and progress bar for intermediar steps, "all" add a plot
#' of the grid to follow visually the processing, "silent" show nothing
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A list of neihbours as weights.
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom rgeos gLength
#' @importFrom utils capture.output
#' @importFrom data.table data.table .SD transpose
#' @examples
#' #no example provided, this is an internal function
network_listw_worker<-function(points,lines,maxdistance,dist_func, direction=NULL, mindist=10, matrice_type = "B", verbose = "silent", digits = 3, tol=0.1){
    #the points must already be snapped to the lines
    #adding the points to the lines
    if(verbose!="silent"){
        print("adding the points as vertices to nearest lines")
    }
    if(verbose=="progressbar"){
        new_lines <- add_vertices_lines(lines, points, tol = 0.1, check = TRUE)
    }else{
        invisible(capture.output(new_lines<-add_vertices_lines(lines, points, tol = 0.1, check = TRUE)))
    }

    #splitting the lines on vertices and adjusting weights
    if(verbose!="silent"){
        print("splitting the lines for the network")
    }
    graph_lines <- simple_lines(new_lines)
    graph_lines$lx_length <- gLength(graph_lines,byid=T)
    graph_lines$lx_weight <- (graph_lines$lx_length / graph_lines$line_length) * graph_lines$line_weight
    #building the network man !
    if(verbose!="silent"){
        print("generating the network")
    }
    # generating the network
    if (is.null(direction)){
        result_graph <- build_graph(graph_lines, digits = digits,
                                    attrs = T, line_weight = "lx_weight")
    }else{
        dir <- ifelse(graph_lines[[direction]]=="Both",0,1)
        result_graph <- build_graph_directed(graph_lines, digits = digits,
                                    attrs = T, line_weight='line_weight',direction = dir)
    }
    #finding all points in the graph
    points$vertex <- find_vertices(result_graph$spvertices,points,digits = digits, tol = tol)
    if(verbose!="silent"){
        print("calculating the distances on the network")
    }
    starts <- subset(points,points$pttype=="start")
    u <- unique(points$vertex)
    #calculating the distances between them
    base_distances <- igraph::distances(result_graph$graph,v = starts$vertex, to = u, mode = "out")
    all_ditancesdt <- data.table(base_distances)
    #all_ditancesdt$origin <- starts[[oid]]
    all_ditancesdt$origin <- starts$fid
    ##1er groupby sur les rows
    step1 <- all_ditancesdt[, lapply(.SD, min, na.rm=TRUE), by=origin ]
    originid <- step1$origin
    ##transposons et merge
    step2 <- transpose(step1[,2:ncol(step1)])
    step2$vertex <- u
    pts <- data.table(points@data[c("fid","vertex")])
    step2[pts,on="vertex",fid:=fid]
    step3 <- step2[, lapply(.SD, min, na.rm=TRUE), by=fid ]
    destid <- step3$fid
    #derniere retransposition
    step4 <- transpose(step3)
    dist_matrix <- as.matrix(step4[2:(nrow(step4)-1)])
    dist_matrix <- ifelse(is.na(dist_matrix),Inf,dist_matrix)
    #now calculating the neighbouring
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



#' Function to prepare selected points and selected lines during the process
#'
#' @param i The index of the quadra to use in the grid
#' @param grid A SpatialPolygonsDataFrame representing the quadras to split
#' calculus
#' @param snapped_points The start and end points snapped to the lines
#' @param lines The lines representing the network
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @return A list of two elements : selected points and selected lines
#' @examples
#' #no example provided, this is an internal function
prepare_elements_netlistw <- function(i,grid,snapped_points,lines,maxdistance){
    #extracting the quadra
    quadra <- grid[i,]
    #selecting the starting points
    test1 <- as.vector(gIntersects(quadra,snapped_points,byid=T))
    if(any(test1)==F){
        return(NULL)
    }else{
        start_pts <- snapped_points[test1,]
        start_pts <- snapped_points[snapped_points$fid %in% start_pts$fid,]
        start_pts$pttype <- "start"
        #selecting the endpoints
        ext <- raster::extent(start_pts)
        poly <- methods::as(ext, "SpatialPolygons")
        raster::crs(poly) <- raster::crs(start_pts)
        buff <- gBuffer(poly, width = maxdistance)
        test2 <- as.vector(gIntersects(buff,snapped_points,byid=T))
        end_pts <- snapped_points[test2 & test1==F,]
        end_pts$pttype <- "end"
        #combining all the points
        all_pts <- rbind(start_pts,end_pts)
        #selecting the lines
        test3 <- as.vector(gIntersects(buff,lines,byid=T))
        selected_lines <- lines[test3,]
        #calculating the elements
        return(list(all_pts,selected_lines))
    }
}



#' generate listw object (spdep like) based on network distances
#'
#' @param origins A SpatialLinesDataFrame, SpatialPointsDataFrame or
#' SpatialPolygonsDataFrame for which the spatial neighbouring list will be built
#' @param lines A SpatialLinesDataFrame representing the network
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @param method A string indicating how the starting points will be built.
#' If centroid is used, then the center of lines or polygons is used. If
#' pointsalong is used, then points will be placed alon polygons' borders or
#' along lines as starting and end points. If ends is used (only for lines)
#' the first and list vertices of lines are used as startng and ends points.
#' @param point_dist A float, defining the distance between points when the
#' method pointsalong is selected.
#' @param snap_dist the maximum distance to snap the start and end points on
#' the network.
#' @param line_weight The ponderation to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional with the geographical length of the lines.
#' @param mindist The minimum distance between two different observations.
#' It is important for it to be different from 0 when a W style is used.
#' @param direction indicate a field giving informations about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
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
#' @param verbose A string indicating how the advance of the process is
#' displayed. Default is all. "text" show only the main steps, "progressbar"
#' show main steps and progress bar for intermediar steps, "all" add a plot
#' of the grid to follow visually the processing, "silent" show nothing
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A listw object (spdep like)
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects gPointOnSurface
#' @importFrom dplyr left_join group_by summarize summarize_all
#' @importFrom graphics plot
#' @export
#' @examples
#' data(mtl_network)
#' listw <- network_listw(mtl_network,mtl_network,maxdistance=500,
#'         method = "centroid", line_weight = "length",
#'         dist_func = 'squared inverse', matrice_type='B', grid_shape = c(2,2))
network_listw <- function(origins,lines,maxdistance, method="centroid", point_dist=NULL, snap_dist=Inf, line_weight = "length", mindist=10, direction=NULL, dist_func = "inverse", matrice_type = "B", grid_shape=c(1,1), verbose = "all", digits = 3, tol=0.1){
    if(verbose %in% c("silent","progressbar","text","all")==F){
        stop("the verbose argument must be 'silent', 'text' or 'progressbar'")
    }
    ##adjusting the weights of the lines
    lines$line_length <- gLength(lines,byid=T)
    if(line_weight=="length"){
        lines$line_weight <- gLength(lines,byid=T)
    }else {
        lines$line_weight <- lines[[line_weight]]
    }
    if(min(lines$line_weight)<=0){
        stop("the weights of the lines must be superior to 0")
    }

    # adjusting the directions of the lines
    if(is.null(direction)==F){
        lines <- lines_direction(lines,direction)
    }

    ## checking the matrix type
    if (matrice_type %in% c("B", "W") == F) {
        stop("Matrice type must be B (binary) or W (row standardized)")
    }
    ## creating the vectorized distance function
    vdist_func <- select_dist_function(dist_func)
    ## setting the OID : tmpid
    lines$tmpid <- 1:nrow(lines)

    ##generating the starting points
    if(verbose != "silent"){
        print("generating the starting points")
    }
    origins$fid <- 1:nrow(origins)
    if(class(origins)=="SpatialPolygonsDataFrame"){
        if(method=="centroid"){
            centers <- gPointOnSurface(origins,byid = T)
            centers <- SpatialPointsDataFrame(centers,origins@data)
        }else if(method=="pointsalong"){
            centers <- surrounding_points(origins,point_dist)
        }
    }else if(class(origins)=="SpatialPointsDataFrame"){
        centers <- origins
    }else if(class(origins)=="SpatialLinesDataFrame"){
        if(method=="centroid"){
            if(verbose %in% c("progressbar","all")){
                centers <- lines_center(origins)
            }else{
                invisible(capture.output(centers <- lines_center(origins)))
            }
        }else if(method=="ends"){
            centers <- lines_extremities(origins)
        }else if(method=="pointsalong"){
            centers <- lines_points_along(origins,point_dist)
        }else{
            stop("if origins class is SpatialLinesDataFrame, method must be on of centroid, ends or pointsalong")
        }
    }else{
        stop("origins must be a SpatialPointsDataFrame, a SpatialPolygonsDataFrame or a SpatialLinesDataFrame")
    }

    if(verbose != "silent"){
        print("snapping the points to the lines (only once)")
    }
    snapped_points <- maptools::snapPointsToLines(centers,lines,maxDist = snap_dist)
    snapped_points <- cbind(snapped_points, centers)

    ##building grid
    grid <- build_grid(grid_shape,origins)


    if (verbose == "all"){
        sp::plot(lines)
        sp::plot(grid,add=T)
    }

    if(verbose != "silent"){
        print("starting the network part")
    }

    ##iterating on the grid
    listvalues <- lapply(1:length(grid),function(i){
        #extracting the quadra
        quadra <- grid[i,]
        if(verbose=="all"){
            sp::plot(quadra,add=T,col="red")
        }
        if(verbose!="silent"){
            print(paste("working on quadra : ",i,"/",length(grid),sep=""))
        }

        elements <- prepare_elements_netlistw(i,grid,snapped_points,lines,maxdistance)
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

    if(verbose != "silent"){
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
    if(verbose != "silent"){
        print("finally generating the listw object ...")
    }
    ## setting the final listw attributes
    listw <- spdep::nb2listw(ordered_nblist, glist = ordered_weights, zero.policy = T,
                             style = matrice_type)
    return(listw)

}

#' generate listw object (spdep like) based on network distances. This version
#' of the function can use a multiprocess plan defined with the package future
#'
#' @param origins A SpatialLinesDataFrame, SpatialPointsDataFrame or
#' SpatialPolygonsDataFrame for which the spatial neighbouring list will be built
#' @param lines A SpatialLinesDataFrame representing the network
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @param method A string indicating how the starting points will be built.
#' If centroid is used, then the center of lines or polygons is used. If
#' pointsalong is used, then points will be placed alon polygons' borders or
#' along lines as starting and end points. If ends is used (only for lines)
#' the first and list vertices of lines are used as startng and ends points.
#' @param point_dist A float, defining the distance between points when the
#' method pointsalong is selected.
#' @param snap_dist the maximum distance to snap the start and end points on
#' the network.
#' @param line_weight The ponderation to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional with the geographical length of the lines.
#' @param mindist The minimum distance between two different observations.
#' It is important for it to be different from 0 when a W style is used.
#' @param direction indicate a field giving informations about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
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
#' @param verbose A string indicating how the advance of the process is
#' displayed. Default is all. "text" show only the main steps, "progressbar"
#' show main steps and progress bar for intermediar steps,"silent" show nothing
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @param tol A float indicating the spatial tolerance when points are
#' added as vertices to lines.
#' @return A listw object (spdep like)
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects gPointOnSurface
#' @importFrom dplyr left_join group_by summarize summarize_all
#' @export
#' @examples
#' data(mtl_network)
#' future::plan(future::multiprocess(workers=4))
#' listw <- network_listw.mc(mtl_network,mtl_network,maxdistance=500,
#'         method = "centroid", line_weight = "length",
#'         dist_func = 'squared inverse', matrice_type='B', grid_shape = c(2,2))
#' \dontshow{
#'    ## R CMD check: make sure any open connections are closed afterward
#'    if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'}
network_listw.mc <- function(origins,lines,maxdistance, method="centroid", point_dist=NULL, snap_dist=Inf, line_weight = "length", mindist=10, direction=NULL, dist_func = "inverse", matrice_type = "B", grid_shape=c(1,1), verbose = "text", digits = 3, tol=0.1){
    if(verbose %in% c("silent","progressbar","text")==F){
        stop("the verbose argument must be 'silent', 'text' or 'progressbar'")
    }
    ##adjusting the weights of the lines
    lines$line_length <- gLength(lines,byid=T)
    if(line_weight=="length"){
        lines$line_weight <- gLength(lines,byid=T)
    }else {
        lines$line_weight <- lines[[line_weight]]
    }
    if(min(lines$line_weight)<=0){
        stop("the weights of the lines must be superior to 0")
    }

    # adjusting the directions of the lines
    if(is.null(direction)==F){
        lines <- lines_direction(lines,direction)
    }

    ## checking the matrix type
    if (matrice_type %in% c("B", "W") == F) {
        stop("Matrice type must be B (binary) or W (row standardized)")
    }
    ## creating the vectorized distance function
    vdist_func <- select_dist_function(dist_func)
    ## setting the OID : tmpid
    lines$tmpid <- 1:nrow(lines)

    ##generating the starting points
    if(verbose != "silent"){
        print("generating the starting points")
    }
    origins$fid <- 1:nrow(origins)
    if(class(origins)=="SpatialPolygonsDataFrame"){
        if(method=="centroid"){
            centers <- gPointOnSurface(origins,byid = T)
            centers <- SpatialPointsDataFrame(centers,origins@data)
        }else if(method=="pointsalong"){
            centers <- surrounding_points(origins,point_dist)
        }
    }else if(class(origins)=="SpatialPointsDataFrame"){
        centers <- origins
    }else if(class(origins)=="SpatialLinesDataFrame"){
        if(method=="centroid"){
            if(verbose %in% c("progressbar","all")){
                centers <- lines_center(origins)
            }else{
                invisible(capture.output(centers <- lines_center(origins)))
            }
        }else if(method=="ends"){
            centers <- lines_extremities(origins)
        }else if(method=="pointsalong"){
            centers <- lines_points_along(origins,point_dist)
        }else{
            stop("if origins class is SpatialLinesDataFrame, method must be on of centroid, ends or pointsalong")
        }
    }else{
        stop("origins must be a SpatialPointsDataFrame, a SpatialPolygonsDataFrame or a SpatialLinesDataFrame")
    }

    if(verbose != "silent"){
        print("snapping the points to the lines (only once)")
    }
    snapped_points <- maptools::snapPointsToLines(centers,lines,maxDist = snap_dist)
    snapped_points <- cbind(snapped_points, centers)

    ##building grid
    grid <- build_grid(grid_shape,origins)

    if(verbose != "silent"){
        print("starting the network part")
    }

    ##iterating on the grid
    listvalues <- future.apply::future_lapply(1:length(grid),function(i){
        ##step1 : preparing elements
        elements <- prepare_elements_netlistw(i,grid,snapped_points,lines,maxdistance)
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

    })

    if(verbose != "silent"){
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
    if(verbose != "silent"){
        print("finally generating the listw object ...")
    }
    ## setting the final listw attributes
    listw <- spdep::nb2listw(ordered_nblist, glist = ordered_weights, zero.policy = T,
                             style = matrice_type)
    return(listw)

}
