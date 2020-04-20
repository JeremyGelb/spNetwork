# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### worker functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Worker function for the nkde_grided.mc function.
#' @param i The row index of the spatiale grid to use
#' @param grid The SpatialPolygons object reprsenting the grid for splitted
#' calculation
#' @param lines The SpatialLinesDataFrame to use as a network
#' @param lixels The SpatialLinesDataFrame representing the lixels
#' @param points The events to use in the kernel density estimation
#' @param kernel_range The range of the kernel function
#' @param kernel The name of the kernel function to use (must be one of
#' quartic, gaussian or epanechnikov). Default is Quartic
#' @param snap_dist The maximum distance snapping between points and lines
#' @param tol The tolerence for topological operations
#' @param digits The number of digits to keep in the coordinates of the
#' geometries
#' @param direction indicate a field giving informations about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @return A SpatialLinesDataFrame of the lines and their kernel densities
#' @importFrom sp coordinates  SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar capture.output
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @examples
#' #This is an internal function, no example provided
exe_nkde <- function(i, grid, lines, lixels, points, kernel, kernel_range, snap_dist, digits, tol, direction) {
    # extracting the analyzed pixels
    quadra <- grid[i]
    test_lixels <- as.vector(gIntersects(lixels, quadra, byid = T))
    if (any(test_lixels) == FALSE) {
        return()
    } else {
        selected_lixels <- subset(lixels, test_lixels)
        if (nrow(selected_lixels) > 0) {
            # extracting the needed lines
            ext <- raster::extent(selected_lixels)
            poly <- methods::as(ext, "SpatialPolygons")
            raster::crs(poly) <- raster::crs(selected_lixels)
            buff <- gBuffer(poly, width = kernel_range)
            test_lines <- as.vector(gIntersects(lines, buff, byid = T))
            selected_lines <- subset(lines, test_lines)
            # extracting the needed points
            test_points <- as.vector(gIntersects(points, buff, byid = T))
            selected_points <- subset(points, test_points)
            #step3 : extract centroids of lixels
            lixelscenters <- gCentroid(selected_lixels, byid = T)
            lixelscenters <- SpatialPointsDataFrame(lixelscenters, selected_lixels@data)
            # step4 : combining events and lixels centers
            lixelscenters$type <- "lixel"
            lixelscenters$OID <- 1:nrow(lixelscenters)

            if (any(test_points) == FALSE) {
                selected_lixels$density <- 0
                return(selected_lixels)
            }

            selected_points$type <- "event"
            selected_points$OID <- 1:nrow(selected_points)
            allpts <- rbind(lixelscenters[c("type", "OID")],
                            selected_points[c("type","OID")])

            # step 5 : snapping all points on lines
            snapped_points <- maptools::snapPointsToLines(allpts, selected_lines,
                maxDist = snap_dist)
            snapped_points$spOID <- sp_char_index(coordinates(snapped_points),
                digits = digits)
            snapped_points <- cbind(snapped_points,allpts)

            # step6 : adding these points as vertices
            invisible(capture.output(newlines <- add_vertices_lines(selected_lines,
                snapped_points, tol = tol)))


            # step7 : converting lines as simple lines
            lines2 <- simple_lines(newlines)

            #now adjusting the weights of the lixels
            lines2$lx_length <- gLength(lines2,byid=T)
            lines2$lx_weight <- (lines2$lx_length / lines2$line_length) * lines2$line_weight

            # step8 : generating the graph
            if (is.null(direction)){
                ListNet <- build_graph(lines2, digits = digits, line_weight='lx_weight')
            }else{
                dir <- ifelse(lines2[[direction]]=="Both",0,1)
                ListNet <- build_graph_directed(lines2, digits = digits,
                                                line_weight='lx_weight',direction = dir)
            }
            Graph <- ListNet$graph

            # step9 : find starting points
            start_pts <- subset(snapped_points, snapped_points$type == "lixel")
            origins <- find_vertices(ListNet$spvertices, start_pts, tol = tol,
                digits = digits)

            end_pts <- subset(snapped_points, snapped_points$type == "event")
            destinations <- find_vertices(ListNet$spvertices, end_pts, tol = tol,
                digits = digits)

            Values <- calc_NKDE(Graph, origins, destinations, kernel = kernel,
                range = kernel_range, weights = selected_points$tmpweight)
            selected_lixels$density <- Values
            return(selected_lixels)
        }
    }
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### launching functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' nkde with multicore support.
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
#' @param line_weight The ponderation to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional with the geographical length of the lines.
#' @param direction indicate a field giving informations about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
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
#' @param verbose A string indicating how the advance of the process is
#' displayed. Default is progressbar "progressbar" shows main steps and
#' progress bar for intermediar steps, "silent" show nothing
#' @return Vector of kernel densities (one for each origin)
#' @importFrom sp coordinates  SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @export
#' @examples
#' data(mtl_network)
#' data(bike_accidents)
#' future::plan(future::multiprocess(workers=2))
#' lixels_nkde <- nkde.mc(mtl_network, bike_accidents,
#'       snap_dist = 150,
#'       lx_length = 150,
#'       line_weight = "length",
#'       mindist = 50,
#'       kernel_range = 800,
#'       kernel='quartic',
#'       grid_shape=c(5,5)
#' )
#' \dontshow{
#'    ## R CMD check: make sure any open connections are closed afterward
#'    if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'}
nkde.mc <- function(lines, points, snap_dist, lx_length, line_weight="length", direction=NULL, kernel_range, kernel = "quartic", tol = 0.1, digits = 3, mindist = NULL, weights = NULL, grid_shape = c(2, 2), verbose = "progressbar") {
    if(verbose %in% c("silent","progressbar")==F){
        stop("the verbose argument must be 'silent' or 'progressbar'")
    }

    #adjusting the weights of the lines
    lines$line_length <- gLength(lines,byid=T)
    if(line_weight=="length"){
        lines$line_weight <- gLength(lines,byid=T)
    }else {
        lines$line_weight <- lines[[line_weight]]
    }
    if(min(lines$line_weight)<=0){
        stop("the weights of the lines must be superior to 0")
    }

    # adjusting the weights if needed
    if (is.null(weights)) {
        W <- rep(1, nrow(points))
    } else {
        W <- points[[weights]]
    }

    # adjusting the directions of the lines
    if(is.null(direction)==F){
        lines <- lines_direction(lines,field = direction)
    }

    show_progress <- verbose=="progressbar"

    points$tmpweight <- W
    if(verbose != "silent"){
        print("generating the grid...")
    }

    grid <- build_grid(grid_shape, spatial = lines)
    ## step2 : generer les lixels
    if(verbose != "silent"){
        print("generating the lixels...")
    }

    lixels <- lixelize_lines.mc(lines, lx_length = lx_length, mindist = mindist,
        show_progress=show_progress)
    lixels$lxid <- 1:nrow(lixels)
    ## step3 : lancer les iterations
    if(verbose != "silent"){
        print("starting the calculation on the grid")
    }
    iseq <- 1:length(grid)
    if (show_progress) {
        progressr::with_progress({
            p <- progressr::progressor(along = iseq)
            values <- future.apply::future_lapply(iseq, function(i) {
                p(sprintf("i=%g", i))
                return(exe_nkde(i, grid, lines, lixels, points, kernel,
                  kernel_range, snap_dist, digits, tol, direction))
            })
        })
    } else {
        values <- future.apply::future_lapply(iseq, exe_nkde, grid = grid,
            lines = lines, lixels = lixels, points = points,
            kernel_range = kernel_range, kernel = kernel,
            snap_dist = snap_dist, digits = digits, tol = tol,
            direction = direction)
    }
    if(verbose != "silent"){
        print("Combining the results from processes ...")
    }
    okvalues <- values[lengths(values) != 0]
    alllixels <- do.call("rbind", okvalues)
    alllixels <- alllixels[!duplicated(alllixels$lxid), ]
    return(alllixels)
}
