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


#' generate listw object (spdep like) based on network distances
#'
#' The distances are calculated from the start and end points of each
#' line to each start and end points of each other lines. The minimum
#' distance is keeped.
#'
#' @param lines A SpatialLinesDataFrame
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @param line_weight The ponderation to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional with the geographical length of the lines.
#' @param dist_func The function to use to convert the distance
#' between observation in spatial weights. Can be 'identity', 'inverse',
#' 'squared inverse' or a function with one parameter x, returning one numeric
#' that will be vectorized internally.
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details.
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(2,2)
#' @param verbose A string indicating how the advance of the process is
#' displayed. Default is all. "text" show only the main steps, "progressbar"
#' show main steps and progress bar for intermediar steps, "all" add a plot
#' of the grid to follow visually the processing, "silent" show nothing
#' @param mindist Indicates the minimum value to replace 0 with. Two lines that
#' share a vertex have a 0 distance between them, wich could create some
#' trouble in the weightings of the neighbouring object.
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @return A listw object (spdep like)
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @importFrom dplyr left_join group_by summarize summarize_all
#' @export
#' @examples
#' data(mtl_network)
#' listw <- line_ext_listw_gridded(mtl_network,maxdistance=800,
#'         line_weight = "length", dist_func = 'squared inverse',
#'         matrice_type='B', grid_shape = c(5,5), mindist = 10)
line_ext_listw_gridded <- function(lines, maxdistance, line_weight = "length", dist_func = "inverse", matrice_type = "B", grid_shape = c(2, 2), verbose = "all", mindist = 10, digits = 3) {
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

    ## checking the matrix type
    if (matrice_type %in% c("B", "W", "C", "U", "minmax", "S") == F) {
        stop("Matrice type must be B, W, C, U, minmax, or S.
             See the documentation of spdep::nb2listw")
    }
    ## creating the vectorized distance function
    vdist_func <- select_dist_function(dist_func)
    ## setting the OID : tmpid
    lines$tmpid <- 1:nrow(lines)
    if(verbose !="silent"){
        print("generating the grid...")
    }
    grid <- build_grid(grid_shape, spatial = lines)

    ## step3 : lancer les iterations
    if(verbose !="silent"){
        print("starting the calculation on the grid")
    }
    remaininglines <- lines
    nquadra <- length(grid)

    if (verbose=="all") {
        sp::plot(grid)
        sp::plot(lines, add = T)
    }

    values <- lapply(1:length(grid), function(i) {
        if(verbose !="silent"){
            print(paste("-------------------Iterating on quadra : ", i, "/",
                        nquadra, "----------------", sep = ""))
        }

        # extracting the analyzed pixels
        quadra <- grid[i]
        if (verbose=="all") {
            sp::plot(quadra, col = "red", add = T)
        }

        # selecting the lines in the quadra
        test_lines <- as.vector(gIntersects(remaininglines, quadra, byid = T))
        if (any(test_lines) == FALSE) {
            return()
        } else {
            if(verbose !="silent"){
                print(paste("---quadra with values : -----", i, sep = ""))
                print("... selecting lines and building graph")
            }
            selected_lines <- subset(remaininglines, test_lines)
            remaininglines <<- subset(remaininglines, test_lines == F)
            ext <- raster::extent(selected_lines)
            poly <- methods::as(ext, "SpatialPolygons")
            raster::crs(poly) <- raster::crs(selected_lines)
            # selecting the lines close to build the network
            buff <- gBuffer(poly, width = maxdistance)
            test_lines2 <- as.vector(gIntersects(lines, buff, byid = T))
            graph_lines <- subset(lines, test_lines2)
            # generating the network
            result_graph <- build_graph(graph_lines, digits = digits,
                attrs = T, line_weight = 'line_weight')
            graph <- result_graph$graph
            graphdf <- igraph::as_long_data_frame(graph)
            # extracting the starting and ending points
            allpts <- lines_extremities(selected_lines)
            allpts_start <- dplyr::left_join(
              subset(allpts@data, allpts@data$pttype == "start"),
              graphdf[c("tmpid", "from")], by = c("tmpid"))
            allpts_start$vertex <- allpts_start$from
            allpts_end <- dplyr::left_join(
              subset(allpts@data, allpts@data$pttype == "end"),
              graphdf[c("tmpid", "to")], by = c("tmpid"))
            allpts_end$vertex <- allpts_end$to
            # now taking the first set of vertices and the second set
            if(verbose !="silent"){
                print("... calculating distances between vertices")
            }
            origin_vertices1 <- allpts_start$vertex
            origin_vertices2 <- allpts_end$vertex
            dest_vertices <- as.numeric(igraph::V(graph))
            # calculating the distances
            alldistances1 <- igraph::distances(graph, v = origin_vertices1,
                to = dest_vertices, mode = "out")
            alldistances2 <- igraph::distances(graph, v = origin_vertices2,
                to = dest_vertices, mode = "out")
            matdist <- ifelse(alldistances1 < alldistances2, alldistances1, alldistances2)
            colnames(matdist) <- dest_vertices
            rownames(matdist) <- selected_lines$tmpid
            if(verbose !="silent"){
                print("... building the sub-distance matrix")
            }

            show_progress <- verbose %in% c("all","progressba")
            if (show_progress){
                pb <- txtProgressBar(min = 0, max = nrow(selected_lines),style = 3)
            }

            # We now have a nice distance matrix
            nblist <- lapply(1:nrow(selected_lines), function(j) {
                # lets get all the distances from the begining and the end of each line
                if (show_progress){
                    setTxtProgressBar(pb, j)
                }
                line <- selected_lines[j, ]
                dfdistances <- data.frame(dest = colnames(matdist), distance = matdist[j,])
                # now filtering the two long distances and NA
                dfdistances <- subset(dfdistances, dfdistances$distance <
                  maxdistance & is.na(dfdistances$distance) == FALSE)
                dfdistances$dest <- as.integer(as.character(dfdistances$dest))
                # finding the corresponding roads
                dfdistances_p1 <- left_join(dfdistances, graphdf, by = c(dest = "from"))
                dfdistances_p2 <- left_join(dfdistances, graphdf, by = c(dest = "to"))
                combined <- rbind(dfdistances_p1[c("tmpid", "distance")],
                  dfdistances_p2[c("tmpid", "distance")])
                combined <- combined %>% group_by(tmpid) %>% summarise_all(min)
                combined <- subset(combined, combined$tmpid != line$tmpid &
                  is.na(combined$tmpid) == F)
                if (nrow(combined) > 0) {
                  if (matrice_type == "B") {
                    combined$weights <- vdist_func(combined$distance)
                    combined$stdweights <- rep(1, length(combined$weights))
                  } else if (matrice_type == "W") {
                    combined$distance <- ifelse(combined$distance == 0, mindist,
                      combined$distance)
                    combined$weights <- vdist_func(combined$distance)
                    combined$stdweights <- combined$weights / sum(combined$weights)
                  } else {
                    combined$distance <- ifelse(combined$distance == 0, mindist,
                      combined$distance)
                    combined$stdweights <- vdist_func(combined$distance)
                  }
                  return(list(combined$tmpid, combined$stdweights))
                } else {
                  return(list(0L, NULL))
                }

            })

            nbs <- lapply(1:length(nblist), function(j) {
                nblist[[j]][[1]]
            })
            names(nbs) <- selected_lines$tmpid
            weights <- lapply(1:length(nblist), function(j) {
                nblist[[j]][[2]]
            })
            names(weights) <- selected_lines$weights
            return(list(nbs, weights))
        }
    })
    # now, combining everything
    if(verbose != "silent"){
        print("combining all the matrices ...")
    }
    okvalues <- values[lengths(values) != 0]
    nblist <- do.call("c", lapply(okvalues, function(value) {
        return(value[[1]])
    }))
    ordered_nblist <- lapply(as.character(lines$tmpid), function(x) {
        return(nblist[[x]])
    })
    nbweights <- do.call("c", lapply(okvalues, function(value) {
        return(value[[2]])
    }))
    names(nbweights) <- names(nblist)
    ordered_weights <- lapply(as.character(lines$tmpid), function(x) {
        return(nbweights[[x]])
    })
    ## setting the nb attributes
    attr(ordered_nblist, "class") <- "nb"
    attr(ordered_nblist, "region.id") <- as.character(lines$tmpid)
    if(verbose != "silent"){
        print("finally generating the listw object ...")
    }
    ## setting the final listw attributes
    listw <- spdep::nb2listw(ordered_nblist, glist = ordered_weights, zero.policy = T,
        style = matrice_type)
    return(listw)
}




#' generate listw object (spdep like) based on network distances
#'
#' The distances are calculated from the centroid of each line to the centroid
#' of each other lines.
#'
#' @param lines A SpatialLinesDataFrame
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @param line_weight The ponderation to use for lines. Default is "length"
#' (the geographical length), but can be the name of a column. The value is
#' considered proportional with the geographical length of the lines.
#' @param dist_func Indicates the function to use to convert the distance
#' between observation in spatial weights. Can be 'identity', 'inverse',
#' 'squared inverse' or a function with one parameter x that will be
#' vectorized internally
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(2,2)
#' @param verbose A string indicating how the advance of the process is
#' displayed. Default is all. "text" show only the main steps, "progressbar"
#' show main steps and progress bar for intermediar steps, "all" add a plot
#' of the grid to follow visually the processing, "silent" show nothing
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @return A listw object (spdep like)
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @importFrom dplyr left_join group_by summarize summarize_all
#' @export
#' @examples
#' data(mtl_network)
#' listw <- line_center_listw_gridded(mtl_network,maxdistance=500,
#'         line_weight = "length",dist_func = 'squared inverse',
#'         matrice_type='B', grid_shape = c(5,5))
line_center_listw_gridded <- function(lines, maxdistance, line_weight = "length", dist_func = "inverse", matrice_type = "B", grid_shape = c(2, 2), verbose = "all", digits = 3) {
    if(verbose %in% c("silent","progressbar","text","all")==F){
        stop("the verbose argument must be 'all', 'silent', 'progressbar' or 'text'")
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

    ## checking the matrix type
    if (matrice_type %in% c("B", "W") == F) {
        stop("Matrice type must be B (binary) or W (row standardized)")
    }
    ## creating the vectorized distance function
    vdist_func <- select_dist_function(dist_func)
    ## setting the OID : tmpid
    lines$tmpid <- 1:nrow(lines)
    ## generating the centers of the lines
    if(verbose != "silent") {
        print("generating the centers of the lines...")
    }
    centers <- lines_center(lines)
    ## adding the centers to lines
    if(verbose != "silent") {
        print("adding these vertices to lines (only once do not worry)...")
    }
    lines <- add_center_lines(lines)
    if(verbose != "silent") {
        print("generating the grid...")
    }
    grid <- build_grid(grid_shape, spatial = lines)
    ## step3 : lancer les iterations
    if(verbose != "silent") {
        print("starting the calculation on the grid")
    }
    remaininglines <- lines
    remainingcenters <- centers
    nquadra <- length(grid)

    if (verbose =="all") {
        sp::plot(grid)
        sp::plot(lines, add = T)
    }

    values <- lapply(1:length(grid), function(i) {
        if (verbose != "silent"){
            print(paste("-------------------Iterating on quadra : ", i, "/",
                        nquadra, "----------------", sep = ""))
        }

        # extracting the analyzed pixels
        quadra <- grid[i]
        if (verbose == "all") {
            sp::plot(quadra, col = "red", add = T)
        }

        # selecting the lines in the quadra
        test_lines <- as.vector(gIntersects(remaininglines, quadra, byid = T))
        if (any(test_lines) == FALSE) {
            return()
        } else {
            if (verbose != "silent"){
                print(paste("---quadra with values : -----", i, sep = ""))
                print("... selecting lines and building graph")
            }
            selected_lines <- subset(remaininglines, test_lines)
            selected_centers <- subset(remainingcenters, test_lines)
            remaininglines <<- subset(remaininglines, test_lines == F)
            remainingcenters <<- subset(remainingcenters, test_lines == F)
            ext <- raster::extent(selected_lines)
            poly <- methods::as(ext, "SpatialPolygons")
            raster::crs(poly) <- raster::crs(selected_lines)
            # selecting the lines close to build the network
            buff <- gBuffer(poly, width = maxdistance)
            test_lines2 <- as.vector(gIntersects(lines, buff, byid = T))
            graph_lines <- subset(lines, test_lines2)
            graph_centers <- subset(centers, test_lines2)
            # generating the simpe lines
            graph_lines2 <- simple_lines(graph_lines)
            graph_lines2$lx_length <- gLength(graph_lines2,byid=T)
            graph_lines2$lx_weight <- (graph_lines2$lx_length / graph_lines2$line_length) * graph_lines2$line_weight

            # generating the network
            result_graph <- build_graph(graph_lines2, digits = digits,
                attrs = T, line_weight = "lx_weight")
            graph <- result_graph$graph
            graphdf <- igraph::as_long_data_frame(graph)
            # finding the interesting vertices
            vertices <- find_vertices(spvertices = result_graph$spvertices,
                points = selected_centers, digits = digits)
            selected_centers$vertex <- vertices
            u <- find_vertices(spvertices = result_graph$spvertices, points = graph_centers,
                digits = digits)
            # calculating the distances
            matdist <- igraph::distances(graph, v = vertices, to = u, mode = "out")

            colnames(matdist) <- u
            rownames(matdist) <- selected_lines$tmpid

            if (verbose != "silent"){
                print("... building the sub-distance matrix")
            }

            show_progress <- verbose %in% c("all","progressbar")

            if (show_progress) {
                pb <- txtProgressBar(min = 0, max = nrow(selected_lines),
                  style = 3)
            }

            # we now have a nice distance matrix
            nblist <- lapply(1:nrow(selected_lines), function(j) {
                # lets get all the distances from the begining and the end of each line
                if (show_progress){
                    setTxtProgressBar(pb, j)
                }
                line <- selected_lines[j, ]
                dfdistances <- data.frame(dest = colnames(matdist), distance = matdist[j,])
                # now filtering the two long distances and NA
                dfdistances <- subset(dfdistances, dfdistances$distance <
                  maxdistance & is.na(dfdistances$distance) == FALSE)
                dfdistances$dest <- as.integer(as.character(dfdistances$dest))
                # finding the corresponding roads
                dfdistances_p1 <- left_join(dfdistances, graphdf, by = c(dest = "from"))
                dfdistances_p2 <- left_join(dfdistances, graphdf, by = c(dest = "to"))
                combined <- rbind(dfdistances_p1[c("tmpid", "distance")],
                  dfdistances_p2[c("tmpid", "distance")])
                combined <- combined %>% group_by(tmpid) %>% summarise_all(min)
                combined <- subset(combined, combined$tmpid != line$tmpid &
                  is.na(combined$tmpid) == F)
                if (nrow(combined) > 0) {
                  if (matrice_type == "B") {
                    combined$weights <- vdist_func(combined$distance)
                    combined$stdweights <- rep(1, length(combined$weights))
                  } else if (matrice_type == "W") {
                    combined$weights <- vdist_func(combined$distance)
                    combined$stdweights <- combined$weights / sum(combined$weights)
                  }
                  return(list(combined$tmpid, combined$stdweights))
                } else {
                  return(list(0L, NULL))
                }

            })

            nbs <- lapply(1:length(nblist), function(j) {
                nblist[[j]][[1]]
            })
            names(nbs) <- selected_lines$tmpid
            weights <- lapply(1:length(nblist), function(j) {
                nblist[[j]][[2]]
            })
            names(weights) <- selected_lines$weights
            return(list(nbs, weights))
        }
    })
    # now, combining everything
    if(verbose != "silent"){
        print("combining all the matrices ...")
    }
    okvalues <- values[lengths(values) != 0]
    nblist <- do.call("c", lapply(okvalues, function(value) {
        return(value[[1]])
    }))
    ordered_nblist <- lapply(as.character(lines$tmpid), function(x) {
        return(nblist[[x]])
    })
    nbweights <- do.call("c", lapply(okvalues, function(value) {
        return(value[[2]])
    }))
    names(nbweights) <- names(nblist)
    ordered_weights <- lapply(as.character(lines$tmpid), function(x) {
        return(nbweights[[x]])
    })
    ## setting the nb attributes
    attr(ordered_nblist, "class") <- "nb"
    attr(ordered_nblist, "region.id") <- as.character(lines$tmpid)
    if(verbose != "silent"){
        print("finally generating the listw object ...")
    }
    ## setting the final listw attributes
    listw <- spdep::nb2listw(ordered_nblist, glist = ordered_weights, zero.policy = T,
        style = matrice_type)
    return(listw)
}
