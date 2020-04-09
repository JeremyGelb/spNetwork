# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### worker functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Worker function for the line_ext_listw_gridded.mc function.
#' @param i The row index of the spatiale grid to use
#' @param lines The SpatialLinesDataFrame to use as a network
#' @param grid The SpatialPolygons object reprsenting the grid for splitted
#' calculation
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @param dist_func Indicates the function to use to convert the distance
#' between observation in spatial weights. Can be 'identity', 'inverse',
#' 'squared inverse' or a function with one parameter x that will be
#' vectorized internally
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(2,2)
#' @param show_progress A boolean indicating if a plot should be displayed to
#' track the progressing of the calculation.
#' @param mindist Indicates the minimum value to replace 0 with. Two lines that
#' share a vertex have a 0 distance between them, wich could create some
#' trouble in the weightings of the neighbouring object.
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @return A list, containing a neighbours list and a weights list
#' @importFrom rgeos gIntersects gBuffer
#' @importFrom dplyr group_by %>% summarise summarise_all
#' @examples
#' #This is an internal function, no example provided
exe_line_ext_listw <- function(i, lines, grid, maxdistance, digits, matrice_type, vdist_func, mindist) {
    # extracting the analyzed pixels
    quadra <- grid[i]
    # selecting the lines in the quadra
    test_lines <- as.vector(gIntersects(lines, quadra, byid = T))
    if (any(test_lines) == FALSE) {
        return()
    } else {
        selected_lines <- subset(lines, test_lines)
        ext <- raster::extent(selected_lines)
        poly <- as(ext, "SpatialPolygons")
        raster::crs(poly) <- raster::crs(selected_lines)
        # selecting the lines close to build the network
        buff <- gBuffer(poly, width = maxdistance)
        test_lines2 <- as.vector(gIntersects(lines, buff, byid = T))
        graph_lines <- subset(lines, test_lines2)
        # generating the network
        result_graph <- build_graph(graph_lines, digits = digits, attrs = T)
        graph <- result_graph$graph
        graphdf <- igraph::as_long_data_frame(graph)
        # extracting the starting and ending points
        allpts <- lines_extremities(selected_lines)
        allpts_start <- dplyr::left_join(
          subset(allpts@data, allpts@data$pttype == "start"),
          graphdf[c("tmpid", "from")], by = c("tmpid"))
        allpts_start$vertex <- allpts_start$from
        allpts_end <- dplyr::left_join(subset(
          allpts@data, allpts@data$pttype == "end"),
          graphdf[c("tmpid", "to")], by = c("tmpid"))
        allpts_end$vertex <- allpts_end$to
        # now taking the first set of vertices and the second set
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
        # on a maintenant une belle matrice de distance entre les origines et les
        # destinations
        nblist <- lapply(1:nrow(selected_lines), function(j) {
            # lets get all the distances from the begining and the end of each line
            line <- selected_lines[j, ]
            dfdistances <- data.frame(dest = colnames(matdist), distance = matdist[j,
                ])
            # now filtering the two long distances and NA
            dfdistances <- subset(dfdistances, dfdistances$distance < maxdistance &
                is.na(dfdistances$distance) == FALSE)
            dfdistances$dest <- as.integer(as.character(dfdistances$dest))
            # finding the corresponding roads
            dfdistances_p1 <- left_join(dfdistances, graphdf, by = c(dest = "from"))
            dfdistances_p2 <- left_join(dfdistances, graphdf, by = c(dest = "to"))
            combined <- rbind(dfdistances_p1[c("tmpid", "distance")],
                dfdistances_p2[c("tmpid","distance")])
            combined <- combined %>% group_by(tmpid) %>% summarise_all(min)
            combined <- subset(combined,
                combined$tmpid != line$tmpid & is.na(combined$tmpid) == F)
            if (nrow(combined) > 0) {
                if (matrice_type == "B") {
                  combined$weights <- vdist_func(combined$distance)
                  combined$stdweights <- rep(1, length(combined$weights))
                } else if (matrice_type == "W") {
                  combined$distance <- ifelse(combined$distance == 0, mindist,
                    combined$distance)
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
}


#' Worker function for the line_center_listw_gridded.mc function.
#' @param i The row index of the spatiale grid to use
#' @param lines The SpatialLinesDataFrame to use as a network
#' @param grid The SpatialPolygons object reprsenting the grid for splitted
#' calculation
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @param dist_func Indicates the function to use to convert the distance
#' between observation in spatial weights. Can be 'identity', 'inverse',
#' 'squared inverse' or a function with one parameter x that will be
#' vectorized internally
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(2,2)
#' @param show_progress A boolean indicating if a plot should be displayed to
#' track the progressing of the calculation.
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @return A list, containing a neighbours list and a weights list
#' @importFrom rgeos gIntersects gBuffer
#' @importFrom dplyr group_by %>% summarise summarise_all
#' @examples
#' #This is an internal function, no example provided
exe_line_center_listw <- function(i, lines, centers, grid, maxdistance, digits, matrice_type, vdist_func) {
    # extracting the analyzed pixels
    quadra <- grid[i]
    # selecting the lines in the quadra
    test_lines <- as.vector(gIntersects(lines, quadra, byid = T))
    if (any(test_lines) == FALSE) {
        return()
    } else {
        selected_lines <- subset(lines, test_lines)
        selected_centers <- subset(centers, test_lines)
        # selecting the lines close to build the network
        ext <- raster::extent(selected_lines)
        poly <- as(ext, "SpatialPolygons")
        raster::crs(poly) <- raster::crs(selected_lines)
        # selecting the lines close to build the network
        buff <- gBuffer(poly, width = maxdistance)
        test_lines2 <- as.vector(gIntersects(lines, buff, byid = T))
        graph_lines <- subset(lines, test_lines2)
        graph_centers <- subset(centers, test_lines2)
        # generating the simpe lines
        graph_lines2 <- simple_lines(graph_lines)
        # generating the network
        result_graph <- build_graph(graph_lines2, digits = digits, attrs = T)
        graph <- result_graph$graph
        graphdf <- igraph::as_long_data_frame(graph)
        # finding the interesting vertices
        vertices <- find_vertices(spvertices = result_graph$spvertices,
            points = selected_centers, digits = digits)
        selected_centers$vertex <- vertices
        # u <- igraph::V(graph)
        u <- find_vertices(spvertices = result_graph$spvertices,
            points = graph_centers, digits = digits)
        # calculating the distances
        matdist <- igraph::distances(graph, v = vertices, to = u, mode = "out")

        colnames(matdist) <- u
        rownames(matdist) <- selected_lines$tmpid

        # on a maintenant une belle matrice de distance entre les origines et les
        # destinations
        nblist <- lapply(1:nrow(selected_lines), function(j) {
            # lets get all the distances from the begining and the end of each line
            line <- selected_lines[j, ]
            dfdistances <- data.frame(dest = colnames(matdist), distance = matdist[j,])
            # now filtering the two long distances and NA
            dfdistances <- subset(dfdistances, dfdistances$distance < maxdistance &
                is.na(dfdistances$distance) == FALSE)
            dfdistances$dest <- as.integer(as.character(dfdistances$dest))
            # finding the corresponding roads
            dfdistances_p1 <- left_join(dfdistances, graphdf, by = c(dest = "from"))
            dfdistances_p2 <- left_join(dfdistances, graphdf, by = c(dest = "to"))
            combined <- rbind(dfdistances_p1[c("tmpid", "distance")],
                dfdistances_p2[c("tmpid","distance")])
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
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### launching functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' generate listw object (spdep like) based on network distances using
#' multiple processes
#'
#' The distances are calculated from the start and end point of each
#' line to each start and end point of each other lines. The minimum
#' distance is keeped.
#'
#' @param lines A SpatialLinesDataFrame
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @param dist_func Indicates the function to use to convert the distance
#' between observation in spatial weights. Can be 'identity', 'inverse',
#' 'squared inverse' or a function with one parameter x that will be
#' vectorized internally
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(2,2)
#' @param show_progress A boolean indicating if a plot should be displayed to
#' track the progressing of the calculation.
#' @param mindist Indicates the minimum value to replace 0 with. Two lines that
#' share a vertex have a 0 distance between them, wich could create some
#' trouble in the weightings of the neighbouring object.
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @return A listw object (spdep like)
#' @importFrom sp coordinates  SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @importFrom dplyr left_join group_by summarize summarize_all
#' @export
#' @examples
#' data(mtl_network)
#' future::plan(future::multiprocess(workers=2))
#' listw <- line_ext_listw_gridded.mc(mtl_network,maxdistance=800,
#'         dist_func = 'squared inverse', matrice_type='B',
#'         grid_shape = c(5,5),
#'         mindist = 10)
line_ext_listw_gridded.mc <- function(lines, maxdistance, dist_func = "inverse", matrice_type = "B", grid_shape = c(2, 2), show_progress = TRUE, mindist = 10, digits = 3) {
    ## checking the matrix type
    if (matrice_type %in% c("B", "W") == F) {
        stop("Matrice type must be B (binary) or W (row standardized)")
    }
    # starting the clusters creating the vectorized distance function
    vdist_func <- select_dist_function(dist_func)
    ## setting the OID : tmpid
    lines$tmpid <- 1:nrow(lines)
    print("generating the grid...")
    grid <- build_grid(grid_shape, spatial = lines)
    ## step3 : startint the iterations
    print("starting the calculation on the grid")
    iseq <- 1:length(grid)
    if (show_progress) {
        progressr::with_progress({
            p <- progressr::progressor(along = iseq)
            values <- future.apply::future_lapply(iseq, function(i) {
                p(sprintf("i=%g", i))
                return(exe_line_ext_listw(i = i, lines = lines, grid = grid,
                  maxdistance = maxdistance, digits = digits, matrice_type = matrice_type,
                  mindist = mindist, vdist_func = vdist_func))
            })
        })
    } else {
        values <- future.apply::future_lapply(iseq, exe_line_ext_listw,
            lines = lines, grid = grid, maxdistance = maxdistance,
            digits = digits, matrice_type = matrice_type,
            mindist = mindist, vdist_func = vdist_func)
    }
    # now, combining everything
    print("combining all the matrices ...")
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
    print("finally generating the listw object ...")
    # setting the final listw attributes
    listw <- spdep::nb2listw(ordered_nblist, glist = ordered_weights,
                    zero.policy = T, style = matrice_type)
    return(listw)
}




#' generate listw object (spdep like) based on network distances using
#' multiple processes
#'
#' The distances are calculated from the centroid of each line to the centroid
#' of each other lines.
#'
#' @param lines A SpatialLinesDataFrame
#' @param maxdistance The maximum distance between two observation to
#' considere them as neighbours.
#' @param dist_func Indicates the function to use to convert the distance
#' between observation in spatial weights. Can be 'identity', 'inverse',
#' 'squared inverse' or a function with one parameter x that will be
#' vectorized internally
#' @param matrice_type The type of the weighting scheme. Can be 'B' for Binary,
#' 'W' for row weighted, see the documentation of spdep::nb2listw for details
#' @param grid_shape A vector of length 2 indicating the shape of the grid to
#' use for splitting the dataset. Default is c(2,2)
#' @param show_progress A boolean indicating if a plot should be displayed to
#' track the progressing of the calculation.
#' @param mindist Indicates the minimum value to replace 0 with. Two lines that
#' share a vertex have a 0 distance between them, wich could create some
#' trouble in the weightings of the neighbouring object.
#' @param digits the number of digits to keep in the spatial coordinates (
#' simplification used to reduce risk of topological error)
#' @return A listw object (spdep like)
#' @importFrom sp coordinates  SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gCentroid gLength gBuffer gIntersects
#' @importFrom dplyr left_join group_by summarize summarize_all
#' @export
#' @examples
#' data(mtl_network)
#' future::plan(future::multiprocess(workers=2))
#' listw <- line_center_listw_gridded.mc(mtl_network,maxdistance=800,
#'         dist_func = 'squared inverse', matrice_type='B',
#'         grid_shape = c(5,5))
line_center_listw_gridded.mc <- function(lines, maxdistance, dist_func = "inverse", matrice_type = "B", grid_shape = c(2, 2), digits = 3, show_progress = TRUE) {
    ## checking the matrix type
    if (matrice_type %in% c("B", "W") == F) {
        stop("Matrice type must be B (binary) or W (row standardized)")
    }
    ## creating the vectorized distance function
    vdist_func <- select_dist_function(dist_func)
    ## setting the OID : tmpid
    lines$tmpid <- 1:nrow(lines)
    ## generating the centers of the lines
    print("generating the centers of the lines...")
    centers <- lines_center(lines)
    ## add vertices to lines
    print("adding these vertices to lines (only once do not worry)...")
    lines <- add_center_lines.mc(lines, show_progress)
    print("generating the grid...")
    grid <- build_grid(grid_shape, spatial = lines)
    ## step3 : start the iteration process
    print("starting the calculation on the grid")
    iseq <- 1:length(grid)
    if (show_progress) {
        progressr::with_progress({
            p <- progressr::progressor(along = iseq)
            values <- future.apply::future_lapply(iseq, function(i) {
                p(sprintf("i=%g", i))
                return(exe_line_center_listw(i = i, lines = lines,
                    centers = centers, grid = grid, maxdistance = maxdistance,
                    digits = digits, matrice_type = matrice_type,
                    vdist_func = vdist_func))
            })
        })
    } else {
        values <- future.apply::future_lapply(iseq, exe_line_center_listw,
            lines = lines, centers = centers, grid = grid,
            maxdistance = maxdistance, digits = digits,
            matrice_type = matrice_type, mindist = mindist,
            vdist_func = vdist_func)
    }
    # now, combining everything
    print("combining all the matrices ...")
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
    print("finally generating the listw object ...")
    ## setting the final listw attributes
    listw <- spdep::nb2listw(ordered_nblist, glist = ordered_weights,
                zero.policy = T, style = matrice_type)
    return(listw)
}
