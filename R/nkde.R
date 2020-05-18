# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### available kernels ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' calcualte the density at a point with a quartic function given
#' the distances of the points around, their weights and the kernel range
#'
#' @param d A vector of distances
#' @param w A vector of weights
#' @param r A float representing the range
#' @return The kernel density
#' @examples
#' #This is an internal function, no example provided
quartic_kernel <- function(d, w, r) {
    w <- w[d < r]
    d <- d[d < r]
    u <- (d**2 / r**2)
    K <- (3 / pi) * w * (1 - u)
    tot <-  sum((1 / r) * K)
    return(sum(w) * tot)
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
gaussian_kernel <- function(d, w, r) {
    w <- w[d < r]
    d <- d[d < r]
    K <- (1 / (sqrt(2 * pi))) * exp(-((d**2)/(2*(r**2)))) * w
    tot <-  sum((1 / r) * K)
    return(sum(w) * tot)
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
epanechnikov_kernel <- function(d, w, r) {
    w <- w[d < r]
    d <- d[d < r]
    u <- (d**2 / r**2)
    K <- (3 / 4) * (1 - u) * w
    tot <-  sum((1 / r) * K)
    return(sum(w) * tot)
}

#' select the right kernel function with its name
#'
#' @param name The name of the kernel to use
#' @return A kernel function
#' @examples
#' #This is an internal function, no example provided
select_kernel <- function(name) {
    if (name == "quartic") {
        return(quartic_kernel)
    } else if (name == "gaussian") {
        return(gaussian_kernel)
    } else if (name == "epanechnikov") {
        return(epanechnikov_kernel)
    } else {
        stop("the specified kernel is not implemented. The kernel name must be
             in c('quartic','gaussian','epanechnikov')")
    }
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### NKDE calculation ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
#' @importFrom data.table data.table
#' @examples
#' #This is an internal function, no example provided
calc_NKDE <- function(graph, origins, destinations, range, kernel = "quartic", weights = NULL) {
    # selecting the kernel function
    kernel_fun <- select_kernel(kernel)
    # adjusting the weights if needed
    if (is.null(weights)) {
        weights <- rep(1, length(destinations))
    }
    df <- data.frame(dest = destinations, weights = weights)
    dt <- data.table(df)
    counts <- dt[,list(sweights=sum(weights)),by=dest]
    #counts <- df %>% dplyr::group_by(dest) %>% dplyr::summarise_all(sum)
    destinations <- counts$dest
    weights <- counts$sweights
    i <- 0
    Values <- sapply(origins, function(o) {
        i < i + 1
        alldistances <- igraph::distances(graph, v = o, to = destinations,
            mode = "out")
        d <- alldistances[alldistances < range]
        w <- weights[alldistances < range]
        v <- kernel_fun(d, w, range)
        return(v)
    })
    return(Values)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Perform nkde analysis ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' nkde function.
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
#' @param lixels A SpatialLinesDataFrame created before with the same network.
#' If NULL, then the lixels will be created using the argument lx_lenght and
#' mindist
#' @param line_weight The ponderation to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional with the geographical length of the lines.
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
#' use for splitting the dataset. Default is c(1,1)
#' @param verbose A string indicating how the advance of the process is
#' displayed. Default is all. "text" show only the main steps, "progressbar"
#' show main steps and progress bar for intermediar steps, "all" add a plot
#' of the grid to follow visually the processing, "silent" show nothing
#' @return Vector of kernel densities (one for each origin)
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar capture.output
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @importFrom data.table data.table
#' @export
#' @examples
#' data(mtl_network)
#' data(bike_accidents)
#' lixels_nkde <- nkde(mtl_network, bike_accidents,
#'       snap_dist = 150,
#'       lx_length = 150,
#'       mindist = 50,
#'       kernel_range = 800,
#'       kernel='quartic',
#'       grid_shape=c(5,5)
#' )
nkde<- function(lines, points, snap_dist, lx_length, lixels=NULL, line_weight="length", kernel_range, kernel = "quartic", tol = 0.1, digits = 3, mindist = NULL, weights = NULL, grid_shape = c(1,1), verbose = "all") {
    if(verbose %in% c("silent","progressbar","text","all")==F){
        stop("the verbose argument must be 'all', 'silent', 'progressbar' or 'text'")
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

    # adjusting the weights of the points
    if (is.null(weights)) {
        W <- rep(1, nrow(points))
    } else {
        W <- points[[weights]]
    }

    show_progress <- verbose %in% c("all","progressbar")

    points$tmpweight <- W
    if (verbose != "silent"){
        print("generating the grid...")
    }
    grid <- build_grid(grid_shape, spatial = lines)

    if (verbose != "silent"){
        print("generating the lixels...")
    }

    if (is.null(lixels)){
        if (show_progress){
            lixels <- lixelize_lines(lines, lx_length = lx_length, mindist = mindist)
        }else {
            invisible(capture.output(lixels <- lixelize_lines(lines,
                lx_length = lx_length, mindist = mindist)))
        }
    }

    lixels$lxid <- 1:nrow(lixels)

    if (verbose != "silent"){
        print("starting the calculation on the grid")
    }
    remaininglixels <- lixels

    if (verbose =="all") {
        sp::plot(grid)
        sp::plot(lines, add = T)
    }
    nquadra <- length(grid)
    values <- lapply(1:nquadra, function(i) {
        if (verbose != "silent"){
            print(paste("-------------------Iterating on quadra : ", i, "/",
                        nquadra, "----------------", sep = ""))
        }
        # extracting the analyzed pixels
        quadra <- grid[i]
        if (verbose =="all") {
            sp::plot(quadra, col = "red", add = T)
        }
        test_lixels <- as.vector(gIntersects(remaininglixels, quadra, byid = T))
        if (any(test_lixels) == FALSE) {
            if (verbose != "silent"){
                print("---- passing, empty quadra----")
            }
            return()
        } else {
            selected_lixels <- subset(remaininglixels, test_lixels)
            if (nrow(selected_lixels) > 0) {
                remaininglixels <<- subset(remaininglixels, test_lixels == F)
                # extracting the needed lines
                ext <- raster::extent(selected_lixels)
                poly <- methods::as(ext, "SpatialPolygons")
                raster::crs(poly) <- raster::crs(selected_lixels)
                # selecting the lines close to build the network
                buff <- gBuffer(poly, width = kernel_range)
                test_lines <- as.vector(gIntersects(lines, buff, byid = T))
                selected_lines <- subset(lines, test_lines)
                # extracting the needed points
                test_points <- as.vector(gIntersects(points, buff, byid = T))
                selected_points <- subset(points, test_points)
                ## performing the real job
                if (verbose != "silent"){
                    print("extracting the centroids of the lixels...")
                }
                # step3 : extraire le centroid ce ces lixels
                lixelscenters <- gCentroid(selected_lixels, byid = T)
                lixelscenters <- SpatialPointsDataFrame(lixelscenters, selected_lixels@data)
                if (verbose != "silent"){
                    print("combining lixels and events ...")
                }
                # step4 : combiner les evenements et les centres de lixels
                lixelscenters$type <- "lixel"
                lixelscenters$OID <- 1:nrow(lixelscenters)

                if (any(test_points) == FALSE) {
                  selected_lixels$density <- 0
                  return(selected_lixels)
                }

                selected_points$type <- "event"
                selected_points$OID <- 1:nrow(selected_points)
                allpts <- rbind(lixelscenters[c("type", "OID")], selected_points[c("type",
                  "OID")])

                # step 5 : snapper ces points sur les lignes
                if (verbose != "silent"){
                    print("snapping points on lines...")
                }
                snapped_points <- maptools::snapPointsToLines(allpts, selected_lines,
                  maxDist = snap_dist)
                snapped_points$spOID <- sp_char_index(coordinates(snapped_points),
                  digits = digits)
                snapped_points <- cbind(snapped_points,allpts)

                # step6 : ajouter ces points comme vertices
                if (verbose != "silent"){
                    print("edditing vertices of lines...")
                }
                if (show_progress){
                    newlines <- add_vertices_lines(selected_lines, snapped_points, tol = tol)
                }else {
                    invisible(capture.output(newlines <- add_vertices_lines(selected_lines,
                        snapped_points, tol = tol)))
                }

                # step7 : couversion vers des lignes simple
                if (verbose != "silent"){
                    print("converting polylines to simple lines...")
                }
                lines2 <- simple_lines(newlines)

                #now adjusting the weights of the lixels
                lines2$lx_length <- gLength(lines2,byid=T)
                lines2$lx_weight <- (lines2$lx_length / lines2$line_length) * lines2$line_weight


                # step8 : generation et dessin du graph
                if (verbose != "silent"){
                    print("building the graph...")
                }
                ListNet <- build_graph(lines2, digits = digits,line_weight='lx_weight')
                Graph <- ListNet$graph

                # step9 : retrouver les points de depart
                if (verbose != "silent"){
                    print("getting the starting points...")
                }
                start_pts <- subset(snapped_points, snapped_points$type ==
                  "lixel")
                origins <- find_vertices(ListNet$spvertices, start_pts, tol = tol,
                  digits = digits)

                if (verbose != "silent"){
                    print("getting the destination points...")
                }
                end_pts <- subset(snapped_points, snapped_points$type ==
                  "event")
                destinations <- find_vertices(ListNet$spvertices, end_pts,
                  tol = tol, digits = digits)

                if (verbose != "silent"){
                    print("calulating the density values...")
                }
                Values <- calc_NKDE(Graph, origins, destinations, range = kernel_range,
                  kernel = kernel, weights = selected_points$tmpweight)
                selected_lixels$density <- Values

                return(selected_lixels)
            }
        }
    })

    okvalues <- values[lengths(values) != 0]
    alllixels <- do.call("rbind", okvalues)
    return(alllixels)
}




