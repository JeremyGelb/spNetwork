
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### helper functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate a character vector based on a coodinates matrix and
#' the maximum number of digits to keep
#'
#' @param coords A n * 2 matrix representing the coordinates
#' @param digits The number of digits to keep from the coordinates
#' @return A vector character vector of length n
#' @importFrom data.table data.table tstrsplit :=
#' @examples
#' #This is an internal function, no example provided
sp_char_index <- function(coords, digits) {
    tempdf <- data.frame(Xs = as.character(coords[, 1]),
                         Ys = as.character(coords[, 2]))
    tempdt <- data.table(tempdf)

    tempdt[,  c("xint", "xdec") := tstrsplit(Xs, ".", fixed=TRUE)]
    tempdt[,  c("yint", "ydec") := tstrsplit(Ys, ".", fixed=TRUE)]

    X <- paste(tempdt$xint, substr(tempdt$xdec, start = 1, stop = digits),
        sep = ".")
    Y <- paste(tempdt$yint, substr(tempdt$ydec, start = 1, stop = digits),
        sep = ".")
    return(paste(X, Y, sep = "_"))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### lines manipulations ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate a SpatialPointsDataFrame with the first and last vertex of each
#' line in a SpatialLinesDataFrame.
#'
#' @param lines A SpatialLinesDataFrame
#' @return A SpatialPointsDataFrame
#' @importFrom sp coordinates
#' @examples
#' #This is an internal function, no example provided
lines_extremities <- function(lines) {
    coords <- coordinates(lines)
    ptcoords <- lapply(coords, function(line) {
        linecoords <- line[[1]]
        ptscoords <- linecoords[c(1, nrow(linecoords)), ]
        return(ptscoords)
    })
    newpts <- do.call(rbind, ptcoords)
    ids <- lapply(1:nrow(lines), function(i) {
        return(c(i, i))
    })
    types <- lapply(1:nrow(lines), function(i) {
        return(c("start", "end"))
    })
    ids <- do.call("c", ids)
    types <- do.call("c", types)
    data <- lines@data[ids, ]
    data$X <- newpts[, 1]
    data$Y <- newpts[, 2]
    data$pttype <- types
    sp::coordinates(data) <- cbind(data$X, data$Y)
    raster::crs(data) <- raster::crs(lines)
    return(data)
}

#' A function to deal with the directions of lines
#'
#' @param lines A SpatialLinesDataFrame
#' @param field Indicate a field giving informations about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @return A SpatialLinesDataFrame
#' @importFrom sp coordinates Line Lines SpatialLines SpatialLinesDataFrame
#' @examples
#' #This is an internal function, no example provided
lines_direction <- function(lines,field){
    listlines <- lapply(1:nrow(lines),function(i){
        line <- lines[i,]
        if(line[[field]]=="TF"){
            coords <- coordinates(line)[[1]][[1]]
            new_line <- Lines(list(Line(coords[nrow(coords):1,])), ID = i)
        }else{
            coords <- coordinates(line)[[1]][[1]]
            new_line <- Lines(list(Line(coords)), ID = i)
        }
        return(new_line)
    })
    spLines <- SpatialLines(listlines)
    df <- SpatialLinesDataFrame(spLines,lines@data,match.ID = F)
    raster::crs(df) <- raster::crs(lines)
    return(df)
}


#' Add vertices (SpatialPoints) to a single line (SpatialLines), may fail
#' if the lines geometries are self intersecting
#'
#' @param line The SpatialLine to modify
#' @param points The SpatialPoints to add to as vertex to the lines
#' @param i The index of the line (for recombining after all lines)
#' @return A matrix of coordinates to build the new line
#' @importFrom rgeos gProject
#' @importFrom sp coordinates SpatialPoints Line Lines
#' @examples
#' #This is an internal function, no example provided
add_vertices <- function(line, points, i) {
    # extract coordinates
    line_coords <- coordinates(line)[[1]][[1]]
    original_distances <- sapply(1:nrow(line_coords),function(i){
        if (i==0){
            return(0)
        }else{
            return(sqrt(sum((line_coords[i,]-line_coords[i-1,])**2)))
        }
    })
    line_coords <- cbind(line_coords, original_distances)
    pt_coords <- coordinates(points)
    # calculate lengths
    lengths <- gProject(line, points)
    pt_coords <- cbind(pt_coords, lengths)
    all_coords <- rbind(line_coords,pt_coords)
    # reorder the coordinate matrix
    ord_coords <- all_coords[order(all_coords[,3]), ]
    return(ord_coords[,1:2])
}



#' Add vertices (SpatialPoints) to many lines (SpatialLines), may fail
#' if the lines geometries are self intersecting
#'
#' @param lines The SpatialLinesDataframe to modify
#' @param points The SpatialPoints to add to as vertex to the lines
#' @param tol The max distance to between lines and points so that they are
#'   added
#' @param check A boolean indicating if the function must check if all points
#' were added
#' @return An object of class SpatialLinesDataFrame (package sp)
#' @importFrom rgeos gIntersects gBuffer gDistance
#' @importFrom sp coordinates SpatialPoints SpatialLinesDataFrame Line
#'   SpatialLines
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' #This is an internal function, no example provided
add_vertices_lines <- function(lines, points, tol = 0.1, check = TRUE) {
    pb <- txtProgressBar(min = 0, max = nrow(lines), style = 3)
    alldistances <- gDistance(lines,points,byid=T)<=tol
    ptscheck <- rep(FALSE,nrow(points))
    new_lines_list <- lapply(1:nrow(lines), function(i) {
        setTxtProgressBar(pb, i)
        line <- lines[i, ]
        testpts <- alldistances[,i]
        if (any(testpts)) {
            okpts <- subset(points,testpts)
            if(check){
                ptscheck <<- ifelse(testpts,TRUE,ptscheck)
            }
            newline <- add_vertices(line, okpts, i)
            return(newline)
        } else {
            sline <- coordinates(line)[[1]][[1]]
            return(sline)
        }

    })
    if (check) {
        # performing a check to ensure that all points were added to the lines
        if(sum(ptscheck)<nrow(points)){
            print(paste("remaining points : ", sum(ptscheck) ,"/" ,nrow(points)), sep = "")
            stop("some points were not added as vertices bro... try to increase
                  the tolerance or to snap these points")
        }
    }

    final_lines <- do.call(raster::spLines,new_lines_list)
    final_lines <- SpatialLinesDataFrame(final_lines, lines@data,match.ID = F)
    raster::crs(final_lines) <- raster::crs(lines)
    return(final_lines)
}



#' Cut a SpatialLines object into lixels with a specified minimal distance
#' may fail if the lines geometries are self intersecting.
#'
#' @param lines The SpatialLinesDataframe to modify
#' @param lx_length The length of a lixel
#' @param mindist The minimum length of a lixel. After cut, if the length of
#' the final lixel is shorter than the minimum distance, then it is added to
#' the previous lixel. if NULL, then mindist = maxdist/10. Note that the
#' segments that are already shorter than the minimum distance are not
#' modified.
#' @return An object of class SpatialLinesDataFrame (package sp)
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom  sp coordinates Line Lines SpatialLines SpatialLinesDataFrame SpatialPoints
#' @importFrom rgeos gLength gInterpolate
#' @export
#' @examples
#' data('mtl_network')
#' lixels <- lixelize_lines(mtl_network,150,50)
lixelize_lines <- function(lines, lx_length, mindist = NULL) {
    if (is.null(mindist)) {
        mindist <- lx_length/10
    }
    pb <- txtProgressBar(min = 0, max = nrow(lines), style = 3)
    cnt <- 1
    newlixels <- lapply(1:nrow(lines), function(i) {
        setTxtProgressBar(pb, i)
        line <- lines[i, ]
        tot_length <- gLength(line)
        if (tot_length < lx_length+mindist) {
            coords <- coordinates(line)
            lixel <- list(Lines(list(Line(coords[[1]][[1]])), ID = cnt))
            cnt <<- cnt + 1
            return(lixel)
        } else {
            # producing the points to snapp on
            distances <- seq(lx_length, tot_length, lx_length)
            if ((tot_length - distances[[length(distances)]]) < mindist) {
                distances <- distances[1:(length(distances) - 1)]
            }
            points <- t(sapply(distances, function(d) {
                return(coordinates(gInterpolate(line, d)))
            }))
            points <- data.frame(x = points[, 1], y = points[, 2], distance = distances,
                type = "cut")
            # extracting the original coordinates
            coords <- SpatialPoints(coordinates(line))
            xy <- coordinates(coords)
            points2 <- data.frame(x = xy[, 1], y = xy[, 2], distance = gProject(line,
                coords), type = "base")
            # merging both and sorting
            allpts <- rbind(points, points2)
            allpts <- allpts[order(allpts$distance), ]
            # and now splitting this motherfucker
            indices <- c(0, which(allpts$type == "cut"), nrow(allpts))
            lixels <- lapply(1:(length(indices) - 1), function(j) {
                pts <- allpts[indices[[j]]:indices[[j + 1]], ]
                lixel <- Lines(list(Line(pts[, 1:2])), ID = cnt)
                cnt <<- cnt + 1
                return(lixel)
            })
            return(lixels)
        }
    })
    oids <- sapply(1:nrow(lines), function(i) {
        return(rep(i, length(newlixels[[i]])))
    })
    oids <- do.call("c", oids)
    new_lines <- SpatialLines(unlist(newlixels))
    new_splines <- SpatialLinesDataFrame(new_lines, lines@data[oids, ], match.ID = F)
    raster::crs(new_splines) <- raster::crs(lines)
    return(new_splines)
}


#' Cut a SpatialLines object into lixels with a specified minimal distance
#' may fail if the lines geometries are self intersecting. This version can
#' use a plan defined with the package future
#'
#' @param lines The SpatialLinesDataframe to modify
#' @param lx_length The length of a lixel
#' @param mindist The minimum length of a lixel. After cut, if the length of
#' the final lixel is shorter than the minimum distance, then it is added to
#' the previous lixel. if NULL, then mindist = maxdist/10
#' @param show_progress A boolean indicating if a progress bar must be displayed
#' @param chunk_size The size of a chunk used for multiprocessing. Default is 100.
#' @return An object of class SpatialLinesDataFrame (package sp)
#' @export
#' @importFrom utils capture.output
#' @examples
#' data('mtl_network')
#' future::plan(future::multiprocess(workers=2))
#' lixels <- lixelize_lines.mc(mtl_network,150,50)
#' \dontshow{
#'    ## R CMD check: make sure any open connections are closed afterward
#'    if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'}
lixelize_lines.mc <- function(lines, lx_length, mindist = NULL, show_progress = T, chunk_size = 100) {
    chunks <- split(1:nrow(lines), rep(1:ceiling(nrow(lines) / chunk_size),
                each = chunk_size, length.out = nrow(lines)))
    chunks <- lapply(chunks,function(x){return(lines[x,])})
    # step2 : starting the function
    iseq <- 1:length(chunks)
    if (show_progress) {
        progressr::with_progress({
            p <- progressr::progressor(along = iseq)
            values <- future.apply::future_lapply(iseq, function(i) {
                p(sprintf("i=%g", i))
                chunk_lines <- chunks[[i]]
                invisible(capture.output(new_lines <- lixelize_lines(chunk_lines,
                    lx_length, mindist)))
                return(new_lines)
            })
        })
    } else {
        values <- future.apply::future_lapply(iseq, function(i) {
            chunk_lines <- chunks[[i]]
            invisible(capture.output(new_lines <- lixelize_lines(chunk_lines,
                lx_length, mindist)))
            return(new_lines)
        })
    }

    final_lines <- do.call("rbind", values)
    return(final_lines)
}


#' Split the polylines of a SpatialLinesDataFrame object in simple lines
#'
#' @param lines The SpatialLinesDataframe to modify
#' @return An object of class SpatialLinesDataFrame (package sp)
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom sp coordinates SpatialLinesDataFrame SpatialLines Lines Line
#' @examples
#' #This is an internal function, no example provided
simple_lines <- function(lines) {
    ## extracting the coordinates of the lines
    allcoords <- coordinates(lines)
    counts <- sapply(1:length(allcoords), function(i) {
        return(nrow(allcoords[[i]][[1]]) - 1)
    })
    oids <- sapply(1:nrow(lines), function(i) {
        return(rep(i, counts[[i]]))
    })
    oids <- do.call("c", oids)

    ## using the coordinates to create newlines
    new_lines <- lapply(1:length(allcoords), function(i) {
        coords <- allcoords[[i]][[1]]
        segment_lines <- lapply(1:(nrow(coords) - 1), function(i) {
            mat <- coords[i:(i + 1), ]
            return(mat)
        })
        return(segment_lines)
    })
    data <- lines@data[oids, ]

    final_lines <- do.call(raster::spLines,unlist(new_lines,recursive = F))
    final_lines <- SpatialLinesDataFrame(final_lines,data, match.ID = F)
    raster::crs(final_lines) <- raster::crs(lines)
    return(final_lines)
}

#' Generate a SpatialPointsDataFrame with line center points. The points are
#' located at middle of the line based on the length of the line
#'
#' @param lines The SpatialLinesDataframe to use
#' @return An object of class SpatialPointsDataFrame (package sp)
#' @importFrom sp coordinates SpatialPointsDataFrame SpatialPoints
#' @importFrom rgeos gInterpolate
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @examples
#' data('mtl_network')
#' centers <- lines_center(mtl_network)
lines_center <- function(lines) {
    pb <- txtProgressBar(min = 0, max = nrow(lines), style = 3)
    listpt <- sapply(1:nrow(lines), function(i) {
        setTxtProgressBar(pb, i)
        line <- lines[i, ]
        pt <- gInterpolate(line, gLength(line)/2)
        return(coordinates(pt))
    })
    ptcoords <- t(listpt)
    centerpoints <- SpatialPointsDataFrame(SpatialPoints(ptcoords), lines@data)
    raster::crs(centerpoints) <- raster::crs(lines)
    return(centerpoints)
}

#' Add to each line of a SpatialLinesDataFrame an additionnal vertex at its
#' center
#'
#' @param lines The SpatialLinesDataframe to use
#' @return An object of class SpatialLinesDataframe (package sp)
#' @importFrom sp coordinates SpatialPointsDataFrame SpatialPoints
#' @importFrom rgeos gInterpolate
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' #This is an internal function, no example provided
add_center_lines <- function(lines) {
    pb <- txtProgressBar(min = 0, max = nrow(lines), style = 3)
    listline <- sapply(1:nrow(lines), function(i) {
        setTxtProgressBar(pb, i)
        line <- lines[i, ]
        pt <- gInterpolate(line, gLength(line) / 2)
        newline <- add_vertices(line, pt, i)
        return(newline)
    })
    final_lines <- SpatialLinesDataFrame(SpatialLines(listline), lines@data,
        match.ID = F)
    raster::crs(final_lines) <- raster::crs(lines)
    return(final_lines)
}


#' Add to each line of a SpatialLinesDataFrame an additionnal vertex at its
#' center (multicore version).
#'
#' @param lines The SpatialLinesDataframe to use
#' @param show_progress A boolean indicating if a progress bar must be displayed
#' @param chunk_size The size of a chunk used for multiprocessing. Default is 100.
#' @return An object of class SpatialLinesDataframe (package sp)
#' @importFrom utils capture.output
#' @examples
#' #This is an internal function, no example provided
#' \dontshow{
#'    ## R CMD check: make sure any open connections are closed afterward
#'    if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'}
add_center_lines.mc <- function(lines, show_progress = T, chunk_size = 100) {
    # step1 : splitting the data into chunks
    chunks <- split(1:nrow(lines), rep(1:ceiling(nrow(lines) / chunk_size), each = chunk_size,
        length.out = nrow(lines)))
    chunks <- lapply(chunks, function(x){return(lines[x,])})
    # step2 : starting the function
    iseq <- 1:length(chunks)
    if (show_progress) {
        progressr::with_progress({
            p <- progressr::progressor(along = iseq)
            values <- future.apply::future_lapply(iseq, function(i) {
                p(sprintf("i=%g", i))
                chunk_lines <- chunks[[i]]
                invisible(capture.output(new_lines <- add_center_lines(chunk_lines)))
                return(new_lines)
            })
        })
    } else {
        values <- future.apply::future_lapply(iseq, function(i) {
            chunk_lines <- chunks[[i]]
            invisible(capture.output(new_lines <- add_center_lines(chunk_lines)))
            return(new_lines)
        })
    }

    final_lines <- do.call("rbind", values)
    return(final_lines)
}



lines_points_along <- function(lines,dist){
    lenghts <- gLength(lines, byid = T)
    list_pts <- lapply(1:nrow(lines),function(i){
        line <- lines[i,]
        line_lenght <- lenghts[i]
        distances <- seq(0,line_lenght,dist)
        pts <- gInterpolate(line,distances)
        return(pts)
    })
    oids <- lapply(1:length(list_pts),function(i){rep(i,length(list_pts[[i]]))})
    oids <- do.call("c",oids)
    all_pts <- do.call(rbind,list_pts)
    data <- lines@data[oids,]
    all_pts <- sp::SpatialPointsDataFrame(all_pts,data)
    return(all_pts)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### functions on polygons ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Generate a SpatialPointsDataFrame by placing points along the border of
#' polygons of a SpatialPolygonsDataFrame
#'
#' @param polygons A SpatialPolygonsDataFrame
#' @param dist The distance between the points
#' @importFrom rgeos gBoundary
#' @return a SpatialPolygonsDataFrame representing the grid
#' @examples
#' #This is an internal function, no example provided
surrounding_points <- function(polygons,dist){
    #extracting the boundaries and their lengths
    boundaries <- gBoundary(polygons, byid=T)
    df <- sp::SpatialLinesDataFrame(boundaries,polygons@data)
    all_pts <- lines_points_along(df,dist)
    return(df)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### gridding function ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate a grid of a specified shape in the bbox of a Spatial object
#'
#' @param grid_shape A numeric vector of length 2 indicating the number
#' of rows and the numbers of columns of the grid
#' @param spatial An object of class SpatialLinesDataFrame (package sp)
#' @return a SpatialPolygonsDataFrame representing the grid
#' @examples
#' #This is an internal function, no example provided
build_grid <- function(grid_shape, spatial) {
    if(prod(grid_shape)==1){
        ext <- raster::extent(spatial)
        poly <- methods::as(ext, "SpatialPolygons")
        raster::crs(poly) <- raster::crs(spatial)
        return(poly)
    }else{
        ## step1 : creating the grid
        box <- sp::bbox(spatial)
        x <- seq(from = box[1, 1], to = box[1, 2], length.out = grid_shape[[1]])
        y <- seq(from = box[2, 1], to = box[2, 2], length.out = grid_shape[[2]])
        xy <- expand.grid(x = x, y = y)
        grid.pts <- sp::SpatialPointsDataFrame(coords = xy, data = xy)
        raster::crs(grid.pts) <- raster::crs(spatial)
        sp::gridded(grid.pts) <- TRUE
        grid <- methods::as(grid.pts, "SpatialPolygons")
        return(grid)
    }

}
