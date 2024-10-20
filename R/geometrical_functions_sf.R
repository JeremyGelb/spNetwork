
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### helper functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#defining some global variables to prevent check error (weird flex but ok)
utils::globalVariables(c("Xs", "Ys", "L1", "lineID"))


#' @title Coordinates to unique character vector
#'
#' @description Generate a character vector based on a coordinates matrix and
#' the maximum number of digits to keep.
#'
#' @param coords A n * 2 matrix representing the coordinates
#' @param digits The number of digits to keep from the coordinates
#' @return A vector character vector of length n
#' @importFrom data.table data.table tstrsplit :=
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
sp_char_index <- function(coords, digits) {
    tempdf <- data.frame(Xs = as.character(coords[, 1]),
                         Ys = as.character(coords[, 2]))
    tempdf$Xs <- ifelse(grepl(".",tempdf$Xs,fixed = TRUE), tempdf$Xs, paste(tempdf$Xs,".0",sep=""))
    tempdf$Ys <- ifelse(grepl(".",tempdf$Ys,fixed = TRUE), tempdf$Ys, paste(tempdf$Ys,".0",sep=""))
    tempdt <- data.table(tempdf)

    tempdt[,  c("xint", "xdec") := tstrsplit(Xs, ".", fixed = TRUE)]
    tempdt[,  c("yint", "ydec") := tstrsplit(Ys, ".", fixed = TRUE)]

    X <- paste(tempdt$xint, substr(tempdt$xdec, start = 1, stop = digits),
        sep = ".")
    X <- gsub("NA","0",X,fixed = TRUE)
    Y <- paste(tempdt$yint, substr(tempdt$ydec, start = 1, stop = digits),
        sep = ".")
    Y <- gsub("NA","0",Y,fixed = TRUE)
    return(paste(X, Y, sep = "_"))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### lines manipulations ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Get lines extremities
#'
#' @description Generate a feature collection of points with the first and last vertex
#'   of each line in a feature collection of linestrings.
#'
#' @param lines A feature collection of linestrings (simple Linestrings)
#' @return A feature collection of points
#' @importFrom sf st_coordinates st_drop_geometry st_crs st_as_sf
#' @importFrom data.table as.data.table .SD .N
#' @export
#' @examples
#' data(mtl_network)
#' points <- lines_extremities(mtl_network)
lines_extremities <- function(lines) {
  coords <- st_coordinates(lines)
  data <- st_drop_geometry(lines)
  remove_it <- FALSE

  if(ncol(data) == 1){
    data$jgtmpoid <- 1:nrow(data)
    remove_it <- TRUE
  }
  df <- as.data.table(coords)
  pts <- as.data.frame(df[, .SD[c(1,.N)], by=L1])
  pts$pttype <- rep(c("start","end"), nrow(pts)/2)

  data <- cbind(data[pts$L1,], pts[c("pttype","X","Y")])

  sf_pts <- st_as_sf(
    x = data,
    coords = c("X","Y"),
    crs = st_crs(lines)
  )
  sf_pts$X <- data$X
  sf_pts$Y <- data$Y
  if(remove_it){
    sf_pts$jgtmpoid <- NULL
  }
  return(sf_pts)
}


#' @title Remove loops
#'
#' @description Remove from a sf object with linestring type geometries the lines that have the
#' same starting and ending point.
#'
#' @param lines A sf object with linestring type geometries
#' @param digits An integer indicating the number of digits to keep for the
#' spatial coordinates
#' @return A sf object with linestring type geometries
#' @importFrom sf st_coordinates
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
remove_loop_lines <- function(lines, digits){
    lines_ext <- lines_extremities(lines)
    starts <- subset(lines_ext,lines_ext$pttype=="start")
    ends <- subset(lines_ext,lines_ext$pttype=="end")
    starts$spoid <- sp_char_index(st_coordinates(starts),digits)
    ends$spoid <- sp_char_index(st_coordinates(ends),digits)
    test <- starts$spoid != ends$spoid
    return(subset(lines,test))
}

#' @title Unify lines direction
#'
#' @description A function to deal with the directions of lines. It ensures
#' that only From-To situation are present by reverting To-From lines. For
#' the lines labelled as To-From, the order of their vertices is reverted.
#'
#' @param lines A sf object with linestring type geometries
#' @param field Indicate a field giving information about authorized
#' travelling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @return A sf object with linestring type geometries
#' @importFrom sf st_reverse
#' @export
#' @examples
#' data(mtl_network)
#' mtl_network$length <- as.numeric(sf::st_length(mtl_network))
#' mtl_network$direction <- "Both"
#' mtl_network[6, "direction"] <- "TF"
#' mtl_network_directed <- lines_direction(mtl_network, "direction")
lines_direction <- function(lines,field){

  lines$tmpjgoid <- 1:nrow(lines)
  test <- lines[[field]] == "TF"
  rev_sub <- lines[test,]
  sub <- lines[!test, ]

  reversed <- st_reverse(rev_sub)
  combined <- rbind(sub,reversed)
  combined <- combined[order(combined$tmpjgoid),]

  return(combined)
}


#' @title Reverse lines
#'
#' @description A function to reverse the order of the vertices of lines
#'
#' @param lines A sf object with linestring type geometries
#' @return A sf object with linestring type geometries
#' @importFrom sf st_coordinates
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
reverse_lines <- function(lines){
    return(st_reverse(lines))
}

#' @title Lines coordinates as list
#'
#' @description A function to get the coordinates of some lines as a list of matrices
#'
#' @param lines A sf object with linestring type geometries
#' @return A list of matrices
#' @importFrom sf st_coordinates
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
lines_coordinates_as_list <- function(lines){
  line_coords <- st_coordinates(lines)
  lineid <- line_coords[,3]
  line_coords <- line_coords[,1:2]
  line_coords <- split(data.frame(line_coords), f = lineid)
  line_coords <- lapply(line_coords,as.matrix)
  return(line_coords)
}


# list_coordinates_as_lines_OLD <- function(coord_list, crs){
#
#   lids <- do.call(c,lapply(1:length(coord_list), function(e){rep(e, nrow(coord_list[[e]]))}))
#   new_coords <- data.frame(do.call(rbind, coord_list))
#   new_coords$lineID <- lids
#   pts <- st_as_sf(new_coords, coords = c(1,2), crs = crs)
#
#   final_lines <- pts %>%
#     group_by(lineID) %>%
#     summarise(do_union = FALSE) %>%
#     st_cast("LINESTRING")
#
#   return(final_lines)
# }


#' @title List of coordinates as lines
#'
#' @description A function to convert a list of matrices to as sf object with linestring geometry type
#'
#' @param coord_list A list of matrices
#' @param crs The CRS to use to create the lines
#' @return A sf object with linestring type geometries
#' @importFrom sf st_coordinates st_cast st_as_sf
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
list_coordinates_as_lines <- function(coord_list, crs){

  lids <- do.call(c,lapply(1:length(coord_list), function(e){rep(e, nrow(coord_list[[e]]))}))
  new_coords <- data.frame(do.call(rbind, coord_list))
  new_coords$lineID <- lids
  colnames(new_coords) <- c('X', 'Y', 'lineID')

  final_lines <- sfheaders::sf_linestring(new_coords, x = "X", y = "Y", linestring_id = "lineID", keep = TRUE)
  final_lines <- st_as_sf(final_lines)
  st_crs(final_lines) <- crs

  return(final_lines)
}



#' @title Add vertices to a feature collection of linestrings
#'
#' @description Add vertices (feature collection of points) to their nearest lines
#'   (feature collection of linestrings), may fail if the line geometries are self intersecting.
#'
#' @param lines The feature collection of linestrings to modify
#' @param points The feature collection of points to add to as vertex to the lines
#' @param nearest_lines_idx For each point, the index of the nearest line
#' @param mindist The minimum distance between one point and the extremity of
#'   the line to add the point as a vertex.
#' @return A feature collection of linestrings
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
add_vertices_lines <- function(lines, points, nearest_lines_idx, mindist) {

    line_coords <- lines_coordinates_as_list(lines)
    pt_coords <- st_coordinates(points)

    new_lines_list <- add_vertices_lines_cpp(as.matrix(pt_coords),
                                             line_coords,
                                             nearest_lines_idx,
                                             mindist)

    final_lines <- list_coordinates_as_lines(new_lines_list, st_crs(lines))
    final_lines <- cbind(final_lines, st_drop_geometry(lines))
    return(final_lines)
}


# the old version is kept for debug
# lixelize_lines <- function(lines, lx_length, mindist = NULL, verbose = FALSE) {
#     if (is.null(mindist)) {
#         mindist <- lx_length/10
#     }
#     if(verbose){
#         pb <- txtProgressBar(min = 0, max = nrow(lines), style = 3)
#     }
#     newlixels <- lapply(1:nrow(lines), function(i) {
#         if(verbose){
#             setTxtProgressBar(pb, i)
#         }
#         line <- lines[i, ]
#         tot_length <- gLength(line)
#         if (tot_length < lx_length+mindist) {
#             coords <- coordinates(line)
#             return(list(coords[[1]][[1]]))
#         } else {
#             # producing the points to snapp on
#             distances <- seq(lx_length, tot_length, lx_length)
#             if ((tot_length - distances[[length(distances)]]) < mindist) {
#                 distances <- distances[1:(length(distances) - 1)]
#             }
#             points <- t(sapply(distances, function(d) {
#                 return(coordinates(gInterpolate(line, d)))
#             }))
#             points <- data.frame(x = points[, 1], y = points[, 2],
#                                  distance = distances,
#                                  type = "cut")
#             # extracting the original coordinates
#             coords <- SpatialPoints(coordinates(line))
#             xy <- coordinates(coords)
#             points2 <- data.frame(x = xy[, 1],
#                                   y = xy[, 2],
#                                   distance = gProject(line,coords),
#                                   type = "base")
#             if(points2$distance[nrow(points2)]==0){
#                 points2$distance[nrow(points2)] <- tot_length
#             }
#             # merging both and sorting
#             allpts <- rbind(points, points2)
#             allpts <- allpts[order(allpts$distance, allpts$type), ]
#             # I should remove duplicates here
#             allpts <- allpts[!duplicated(allpts[1:3]),]
#             # and now splitting this motherfucker
#             indices <- c(0, which(allpts$type == "cut"), nrow(allpts))
#             lixels <- lapply(1:(length(indices) - 1), function(j) {
#                 pts <- allpts[indices[[j]]:indices[[j + 1]], ]
#                 return(as.matrix(pts[, 1:2]))
#             })
#             return(lixels)
#         }
#     })
#     oids <- lapply(1:nrow(lines), function(i) {
#         return(rep(i, length(newlixels[[i]])))
#     })
#     oids <- do.call("c", oids)
#
#     new_lines <- do.call(raster::spLines,unlist(newlixels,recursive = FALSE))
#     df <- lines@data[oids, ]
#     if(class(df) != "dataframe"){
#         df <- data.frame("oid" = df)
#     }
#     new_splines <- SpatialLinesDataFrame(new_lines, df, match.ID = FALSE)
#     raster::crs(new_splines) <- raster::crs(lines)
#     return(new_splines)
#
# }

#' @title Cut lines into lixels
#'
#' @description Cut the lines of a feature collection of linestrings into lixels with a specified minimal
#'   distance may fail if the line geometries are self intersecting.
#'
#' @param lines The sf object with linestring geometry type to modify
#' @param lx_length The length of a lixel
#' @param mindist The minimum length of a lixel. After cut, if the length of the
#'   final lixel is shorter than the minimum distance, then it is added to the
#'   previous lixel. if NULL, then mindist = maxdist/10. Note that the segments
#'   that are already shorter than the minimum distance are not modified.
#' @return An sf object with linestring geometry type
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' lixels <- lixelize_lines(mtl_network,150,50)
#' }
lixelize_lines <- function(lines, lx_length, mindist = NULL) {

    if (is.null(mindist)) {
        mindist <- lx_length/10
    }

    if(lx_length <= mindist){
        stop("The lixel length (lx_length) must be greater than minimum length (mindist)")
    }

    coords <- lines_coordinates_as_list(lines)
    result <- lixelize_lines_cpp(coords, lx_length, mindist)

    final_lines <- list_coordinates_as_lines(result[[1]], st_crs(lines))

    df <- st_drop_geometry(lines)[result[[2]]+1, ]


    if(inherits(df, "data.frame") == FALSE){
        df <- as.data.frame(df)
    }

    final_lines <- cbind(final_lines, df)

    return(final_lines)
}



#'@title Cut lines into lixels (multicore)
#'
#' @description Cut the lines of a feature collection of linestrings into lixels with a specified minimal
#'  distance may fail if the line geometries are self intersecting with
#'  multicore support.
#'
#'@param lines A feature collection of linestrings to convert to lixels
#'@param lx_length The length of a lixel
#'@param mindist The minimum length of a lixel. After cut, if the length of the
#'  final lixel is shorter than the minimum distance, then it is added to the
#'  previous lixel. If NULL, then mindist = maxdist/10
#'@param verbose A Boolean indicating if a progress bar must be displayed
#'@param chunk_size The size of a chunk used for multiprocessing. Default is
#'  100.
#'@return A feature collection of linestrings
#'@export
#'@importFrom utils capture.output
#' @examples
#' \donttest{
#' data(mtl_network)
#' future::plan(future::multisession(workers=1))
#' lixels <- lixelize_lines.mc(mtl_network,150,50)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")){
#' future::plan(future::sequential)
#' }
#'}
lixelize_lines.mc <- function(lines, lx_length, mindist = NULL, verbose = TRUE, chunk_size = 100) {
    chunks <- split(1:nrow(lines), rep(1:ceiling(nrow(lines) / chunk_size),
                each = chunk_size, length.out = nrow(lines)))
    chunks <- lapply(chunks,function(x){return(lines[x,])})
    # step2 : starting the function
    iseq <- 1:length(chunks)
    if (verbose) {
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


#' @title LineString to simple Line
#'
#' @description Split the polylines of a feature collection of linestrings in simple
#' segments at each vertex. The values of the columns are duplicated for each segment.
#'
#' @param lines The featue collection of linestrings to modify
#' @return An featue collection of linestrings
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom sf st_crs st_drop_geometry
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' new_lines <- simple_lines(mtl_network)
#' }
simple_lines <- function(lines) {
    ## extracting the coordinates of the lines
    allcoords <- lines_coordinates_as_list(lines)

    counts <- sapply(1:length(allcoords), function(i) {
        return(nrow(allcoords[[i]]) - 1)
    })
    oids <- sapply(1:nrow(lines), function(i) {
        return(rep(i, counts[[i]]))
    })
    dim(oids) <- NULL # just to be sure !
    if(nrow(lines)>1 & is.list(oids)){
        oids <- do.call("c", oids)
    }

    ## using the coordinates to create newlines
    new_lines <- lapply(1:length(allcoords), function(i) {
        coords <- allcoords[[i]]
        segment_lines <- lapply(1:(nrow(coords) - 1), function(i) {
            mat <- coords[i:(i + 1), ]
            return(mat)
        })
        return(segment_lines)
    })
    new_lines <- unlist(new_lines, recursive = FALSE)
    final_lines <- list_coordinates_as_lines(new_lines, st_crs(lines))

    df2 <- data.frame(st_drop_geometry(lines)[oids,])
    names(df2) <- names(st_drop_geometry(lines))
    final_lines <- cbind(final_lines, df2)

    return(final_lines)
}

#' @title Centre points of lines
#'
#' @description Generate a feature collection of points at the centre of the
#'   lines of a feature collection of linestrings. The length of the lines is
#'   used to determine their centres.
#'
#' @param lines A feature collection of linestrings to use
#' @return A feature collection of points
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' centers <- lines_center(mtl_network)
#' }
lines_center <- function(lines) {

    lines$baseorder <- 1:nrow(lines)

    all_lengths <- as.numeric(st_length(lines))
    no_length <- subset(lines, all_lengths == 0)
    with_length <- subset(lines, all_lengths>0)

    ## preparing the new points for lines with a length > 0
    if(nrow(with_length) > 0){
      lines_coords <- lines_coordinates_as_list(with_length)
      pts_coords1 <- points_at_lines_centers_cpp(lines_coords)
      data1 <- data.frame(st_drop_geometry(with_length))
      names(data1) <- names(st_drop_geometry(with_length))
      data1$X <- pts_coords1[,1]
      data1$Y <- pts_coords1[,2]
      pts1 <- st_as_sf(data1, coords = c("X","Y"), crs = st_crs(lines))
    }else{
      pts1 <- NULL
    }

    ## preparing the new points for lines with a length == 0
    if(nrow(no_length) > 0){
      lines_coords <- lines_coordinates_as_list(no_length)
      pts_coords2 <- lapply(lines_coords, function(i){
        return(i[1,])
      })
      pts_coords2 <- do.call(rbind, pts_coords2)
      data2 <- data.frame(st_drop_geometry(no_length))
      names(data2) <-names(st_drop_geometry(no_length))
      data2$X <- pts_coords2[,1]
      data2$Y <- pts_coords2[,2]
      pts2 <- st_as_sf(data2, coords = c("X","Y"), crs = st_crs(lines))
    }else{
      pts2 <- NULL
    }

    if(is.null(pts1) == FALSE & is.null(pts2) == FALSE){
      all_pts <- rbind(pts1,pts2)
    }else if (is.null(pts1) & is.null(pts2) == FALSE){
      all_pts <- pts2
    } else if (is.null(pts1) == FALSE & is.null(pts2)){
      all_pts <- pts1
    }else{
      stop("impossible to create center points for the provided lines... Check your data")
    }

    all_pts <- all_pts[order(all_pts$baseorder),]
    all_pts$baseorder <- NULL

    return(all_pts)
}

#' @title Add center vertex to lines
#'
#' @description Add to each feature of a feature collection of lines an additional vertex at its
#' center.
#'
#' @param lines The feature collection of linestrings to use
#' @return A feature collection of points
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
add_center_lines <- function(lines) {
    coords <- lines_coordinates_as_list(lines)
    new_coords <- add_center_lines_cpp(coords)

    final_lines <- list_coordinates_as_lines(new_coords, crs = st_crs(lines))
    df <- data.frame(st_drop_geometry(lines))
    names(df) <- names(st_drop_geometry(lines))
    final_lines <- cbind(final_lines, df)

    return(final_lines)
}


#' @title Points along lines
#'
#' @description Generate a feature collection of points along the lines of
#' feature collection of Linestrings.
#'
#' @param lines A feature collection of linestrings to use
#' @param dist The distance between the points along the lines
#' @return A feature collection of points
#' @importFrom utils capture.output
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' new_pts <- lines_points_along(mtl_network,50)
#' }
lines_points_along <- function(lines,dist){

    coords <- lines_coordinates_as_list(lines)
    new_coords <- points_along_lines_cpp(coords, dist)

    df <- data.frame(st_drop_geometry(lines))
    deleteit <- FALSE
    if(ncol(df) == 1){
      df$tmpjgid <- 0
      deleteit <- TRUE
    }
    df <- df[new_coords[,3]+1,]
    if(deleteit){
      df$tmpjgid <- NULL
    }
    df$X <- new_coords[,1]
    df$Y <- new_coords[,2]
    all_pts <- st_as_sf(df, coords = c("X","Y"), crs = st_crs(lines))

    return(all_pts)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### functions on polygons ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Points along polygon boundary
#'
#' @description Generate a feature collection of points by placing points along the border of
#' polygons of a feature collection.
#'
#' @param polygons A feature collection of polygons
#' @param dist The distance between the points
#' @return A feature collection of points representing the points arrond the polygond
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
surrounding_points <- function(polygons,dist){
  boundaries <- st_boundary(polygons)
  pts <- lines_points_along(boundaries,dist)
  return(pts)
}

# this is the previous version of the function, kept for debuging
# lines_points_along <- function(lines,dist){
#     lenghts <- gLength(lines, byid = TRUE)
#     list_pts <- lapply(1:nrow(lines),function(i){
#         line <- lines[i,]
#         line_lenght <- lenghts[i]
#         distances <- seq(0,line_lenght,dist)
#         pts <- gInterpolate(line,distances)
#         return(pts)
#     })
#     oids <- lapply(1:length(list_pts),function(i){rep(i,length(list_pts[[i]]))})
#     oids <- do.call("c",oids)
#     all_pts <- do.call(rbind,list_pts)
#     # adding a useless column to avoid a bug when lines has only one column
#     lines$tmpOID <- 1:nrow(lines)
#     data <- lines@data[oids,]
#     all_pts <- sp::SpatialPointsDataFrame(all_pts,data)
#     all_pts$tmpOID <- NULL
#     return(all_pts)
# }

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### gridding function ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title sf geometry bbox
#'
#' @description Generate polygon as the bounding box of a feature collection
#'
#' @param x A feature collection
#' @return A feature collection of polygons
#' @importFrom sf st_as_sfc st_bbox
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
st_bbox_geom <- function(x){
  vec <- st_bbox(x)
  poly <- st_as_sf(st_as_sfc(vec),crs = st_crs(x))
  poly$oid <- 1
  return(poly)
}


#' @title Spatial grid
#'
#' @description Generate a grid of a specified shape in the bbox of a Spatial object.
#'
#' @param grid_shape A numeric vector of length 2 indicating the number
#' of rows and the numbers of columns of the grid
#' @param spatial A list of spatial feature collections objects (package sf)
#' @return A feature collection of polygons representing the grid
#' @importFrom sf st_make_grid
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
build_grid <- function(grid_shape, spatial) {
    boxes <- lapply(spatial,function(i){matrix(st_bbox(i),nrow = 2, byrow = FALSE)})
    boxes <- do.call(cbind,boxes)
    v1 <- as.numeric(c(min(boxes[1,]),max(boxes[1,])))
    v2 <- as.numeric(c(min(boxes[2,]),max(boxes[2,])))
    box <- rbind(v1,v2)
    vec <- c(box[1,1], box[2,1], box[1,2], box[2,2])
    names(vec) <- c("xmin", "ymin", "xmax", "ymax")
    class(vec) <- "bbox"
    if(prod(grid_shape)==1){
      poly <- st_as_sf(st_as_sfc(vec),crs = st_crs(spatial[[1]]))
      poly$oid <- 1
      return(poly)
    }else{
        ## step1 : creating the grid
        grid <- st_as_sf(st_make_grid(vec, n = grid_shape,
                             crs = st_crs(spatial[[1]]),
                             what = "polygons"))
        grid$oid <- 1:nrow(grid)
        return(grid)
    }
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### snapping function ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Nearest point on segment
#'
#' @description Find the nearest projected point on a segment (from maptools)
#'
#' @param s The coordinates of the segment
#' @param p The coordinates of the point
#'
#' @return A numeric vector with the coordinates of the projected point
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
nearestPointOnSegment <- function(s, p){
    # Adapted from http://pastebin.com/n9rUuGRh
    ap <- c(p[1] - s[1,1], p[2] - s[1,2])
    ab <- c(s[2,1] - s[1,1], s[2,2] - s[1,2])
    t <- sum(ap*ab) / sum(ab*ab)
    t <- ifelse(t<0,0,ifelse(t>1,1,t))
    # when start and end of segment are identical t is NA
    t <- ifelse(is.na(t), 0, t)
    x <- s[1,1] + ab[1] * t
    y <- s[1,2] + ab[2] * t
    # Return nearest point and distance
    result <- c(x, y, sqrt((x-p[1])^2 + (y-p[2])^2))
    names(result) <- c("X","Y","distance")
    return(result)
}


#' @title Nearest point on Line
#'
#' @description Find the nearest projected point on a LineString (from maptools)
#'
#' @param coordsLine The coordinates of the line (matrix)
#' @param coordsPoint The coordinates of the point
#'
#' @return A numeric vector with the coordinates of the projected point
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
nearestPointOnLine <- function(coordsLine, coordsPoint){
    nearest_points <- vapply(2:nrow(coordsLine), function(x){
        nearestPointOnSegment(coordsLine[(x-1):x,], coordsPoint)
    }, FUN.VALUE=c(0,0,0))

    # Return coordinates of the nearest point on this line
    return(nearest_points[1:2, which.min(nearest_points[3,])])

}

# NOTE: this code is kept if bugs are encountered in the new version using
# BH instead of sf.
# @title Nearest feature
#
# @description Find the nearest feature from set X in set Y
#
# @param x A feature collection of points
# @param y A feature collection of points
#
# @return A numeric vector with the index of the nearest features
# @keywords internal
# @examples
# #This is an internal function, no example provided
# nearest <- function(x,y){
#     sf1 <- sf::st_as_sf(x)
#     sf2 <- sf::st_as_sf(y)
#     sf::st_crs(sf1) <- sf::st_crs(sf2)
#     best <- sf::st_nearest_feature(sf1,sf2)
#     return(best)
# }

#' @title Nearest line for points
#' @description Find for each point its nearest LineString
#' @param points A feature collection of points
#' @param lines A feature collection of linestrings
#' @param snap_dist A distance (float) given to find for each point its
#' nearest line in a spatial index. A too big value will produce
#' unnecessary distance calculations and a too short value will lead to
#' more iterations to find neighbours. In extrem cases, a too short value
#' could lead to points not associated with lines (index = -1).
#' @param max_iter An integer indicating how many iteration the search
#' algorithm must perform in the spatial index to find lines close to a
#' point. At each iteration, the snap_dist is doubled to find candidates.
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
nearest_lines <- function(points, lines, snap_dist = 300, max_iter = 10){

    # getting the coordinates of the lines
    list_lines <- lines_coordinates_as_list(lines)

    # getting the coordinates of the points
    coords <- st_coordinates(points)

    # getting the indexes
    idx <- find_nearest_object_in_line_rtree(coords, list_lines, snap_dist, max_iter)

    # adding 1 to match with c++ indexing
    return(idx+1)
}


# snapPointsToLines2_OLD <- function(points, lines ,idField = NA, snap_dist = 300, max_iter = 10) {
#
#     #nearest_line_index <- nearest(points,lines)
#     if(is.na(idField)){
#       lines$tmpjgid <- 1:nrow(lines)
#       idField <- "tmpjgid"
#     }
#
#     nearest_line_index <- nearest_lines(points, lines, snap_dist, max_iter)
#     coordsLines <- lines_coordinates_as_list(lines)
#     coordsPoints <- st_coordinates(points)
#
#     # Get coordinates of nearest points lying on nearest lines
#     mNewCoords <- vapply(1:nrow(points),function(x){
#             return(nearestPointOnLine(coordsLines[[nearest_line_index[x]]],
#                                coordsPoints[x,]))}, FUN.VALUE=c(0,0))
#
#     # Recover lines' Ids (If no id field has been specified, take the sp-lines id)
#     nearest_line_id <- lines[[idField]][nearest_line_index]
#
#     # Create data frame and sp points
#     df <- data.frame(st_drop_geometry(points))
#
#     df$nearest_line_id <- nearest_line_id
#     mNewCoords <- t(mNewCoords)
#     df$Xx <- mNewCoords[,1]
#     df$Yy <- mNewCoords[,2]
#     final_points <- st_as_sf(df, coords = c("Xx","Yy"), crs = st_crs(lines))
#     final_points$Xx <- NULL
#     final_points$Yy <- NULL
#
#     return(final_points)
#
# }



#' @title Snap points to lines
#'
#' @description Snap points to their nearest lines (edited from maptools)
#'
#' @param points A feature collection of points
#' @param lines A feature collection of linestrings
#' @param idField The name of the column to use as index for the lines
#' @param snap_dist A distance (float) given to find for each point its
#' nearest line in a spatial index. A too big value will produce
#' unnecessary distance calculations and a too short value will lead to
#' more iterations to find neighbours. In extrem cases, a too short value
#' could lead to points not associated with lines (index = -1).
#' @param max_iter An integer indicating how many iteration the search
#' algorithm must perform in the spatial index to find lines close to a
#' point. At each iteration, the snap_dist is doubled to find candidates.
#'
#' @return A feature collection of points with the projected geometries
#' @keywords internal
#' @importFrom methods slot
#' @export
#' @examples
#' # reading the data
#' data(mtl_network)
#' data(bike_accidents)
#' mtl_network$LineID <- 1:nrow(mtl_network)
#' # snapping point to lines
#' snapped_points <- snapPointsToLines2(bike_accidents,
#'     mtl_network,
#'     "LineID"
#' )
snapPointsToLines2 <- function(points, lines ,idField = NA, snap_dist = 300, max_iter = 100) {

  #nearest_line_index <- nearest(points,lines)
  if(is.na(idField)){
    lines$tmpjgid <- 1:nrow(lines)
    idField <- "tmpjgid"
  }

  #nearest_line_index <- nearest_lines(points, lines, snap_dist, max_iter)

  # we start by finding the nearest line for each point
  nearest_line_index <- sf::st_nearest_feature(points, lines)
  lines2 <- lines[nearest_line_index,]

  # we can then project the points on theses lines
  proj_lines <- sf::st_nearest_points(points, lines2, pairwise = TRUE, tolerance = 1e-9)
  proj_pts <- st_cast(proj_lines, "POINT")
  final_points <- st_as_sf(proj_pts[seq(2,length(proj_pts),2)])

  # we can then add the ids of the lines
  final_points$nearest_line_id <- lines2[[idField]]

  # and we can create the final object
  final_points <- cbind(final_points, st_drop_geometry(points))
  st_geometry(final_points) <- 'geometry'

  return(final_points)

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Network Cleaning Functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Heal edges
#'
#' @description Merge Lines if they form a longer linestring without external intersections (experimental)
#'
#' @param lines A feature collection of linestrings
#' @param digits An integer indicating the number of digits to keep in coordinates
#' @param verbose A boolean indicating if a progress bar should be displayed
#'
#' @return A feature collection of linestrings with the eventually merged geometries. Note
#' that if lines are merged, only the attributes of the first line are preserved
#' @keywords internal
#' @importFrom sf st_line_merge st_sfc st_cast st_geometry st_geometry<-
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' #This is an internal function, no example provided
heal_edges <- function(lines,digits = 3, verbose = TRUE){
    ## healing step
    geoms <- st_geometry(lines)
    lines <- st_drop_geometry(lines)
    sf::st_geometry(lines) <- geoms

    #first getting the coordinates of each point
    lines$tmpOID <- 1:nrow(lines)
    coords <- lines_extremities(lines)
    coords$spindex <- sp_char_index(st_coordinates(coords),digits)

    #then counting for each point how many times it appears
    counts <- table(coords$spindex)
    countdf <- data.frame(
        count = as.vector(counts),
        spindex = names(counts)
    )

    #if a point appears exactly twice, it might require a healing
    cases <- subset(countdf, countdf$count==2)
    consid_points <- subset(coords, coords$spindex %in% cases$spindex)

    #adding the start and end sp index to each line
    startpts <- data.frame(st_drop_geometry(subset(coords, coords$pttype=="start")))
    endpts <- data.frame(st_drop_geometry(subset(coords, coords$pttype=="end")))
    names(startpts) <- names(st_drop_geometry(coords))
    names(endpts) <- names(startpts)

    tempDT <- as.data.table(lines)
    lines$startidx <- setDT(tempDT)[startpts, on = "tmpOID", "startidx" := startpts$spindex][]$startidx
    lines$endidx <- setDT(tempDT)[endpts, on = "tmpOID", "endidx" := endpts$spindex][]$endidx


    test <- (lines$startidx %in% consid_points$spindex | lines$endidx %in% consid_points$spindex)
    consid_lines <- subset(lines, test)
    keeped_lines <- subset(lines, test == FALSE)

    consid_lines$tmpOID <- as.character(consid_lines$tmpOID)
    # generating a neighbouring list from the extremies

    neighbouring <- lapply(1:nrow(consid_lines), function(i){
        line <- consid_lines[i,]
        val1 <- line$startidx
        val2 <- line$endidx
        if(val1 %in% consid_points$spindex){
            test1 <- (consid_lines$startidx == val1 | consid_lines$endidx == val1) &
                               ((consid_lines$startidx == val1 & consid_lines$endidx == val2)==FALSE
            )
        }else{
            test1 <- rep(FALSE, nrow(consid_lines))
        }
        if(val2 %in% consid_points$spindex){
            test2 <- (consid_lines$startidx == val2 | consid_lines$endidx == val2) &
                ((consid_lines$startidx == val1 & consid_lines$endidx == val2)==FALSE
                )
        }else{
            test2 <- rep(FALSE, nrow(consid_lines))
        }
        neighbours <- subset(consid_lines, test1 | test2)

        codes <- neighbours$tmpOID
        return(codes)
    })
    names(neighbouring) <- consid_lines$tmpOID

    ## ok now we must find the routes that we should merge
    merged <- c()
    if(verbose){
        pb <- txtProgressBar(min = 0, max = nrow(consid_lines), style = 3)
    }
    merge_with <- lapply(1:nrow(consid_lines), function(i){
        line <- consid_lines[i,]
        if(verbose){
            setTxtProgressBar(pb, i)
        }
        all_neighbours <- c()
        next_check <- line$tmpOID
        while (length(next_check) > 0){
            neighbours <- do.call(c,lapply(next_check, function(x){
                neighbouring[[x]]
            }))
            neighbours <- neighbours[(neighbours %in% merged) == FALSE]
            next_check <- neighbours
            merged <<- c(merged,neighbours)
            all_neighbours <- c(all_neighbours,neighbours)
        }
        return(all_neighbours)
    })
    merge_with <- Filter(function(x) length(x) > 0, merge_with)

    merged_lines <- st_sfc(do.call(rbind,lapply(merge_with, function(x){
        sub <- subset(consid_lines, consid_lines$tmpOID %in% x)
        return(st_line_merge(st_cast(st_geometry(sub),"MULTILINESTRING", ids = rep(1, nrow(sub)))))
    })))

    df <- do.call(rbind,lapply(merge_with, function(x){
        sub <- subset(consid_lines, consid_lines$tmpOID %in% x)
        return(sub[1,])
    }))

    df$geometry <- merged_lines
    df$startidx <- NULL
    df$endidx <- NULL
    keeped_lines$startidx <- NULL
    keeped_lines$endidx <- NULL
    st_crs(df) <- st_crs(keeped_lines)
    new_lines <- rbind(keeped_lines, df)
    new_lines$tmpOID <- NULL
    return(new_lines)
}

#' @title Remove mirror edges
#'
#' @description Keep unique edges based on start and end point
#'
#' @param lines A feature collection of linestrings
#' @param keep_shortest A boolean, if TRUE, then the shortest line is keeped if
#' several lines have the same starting point and ending point. if FALSE, then the
#' longest line is keeped.
#' @param digits An integer indicating the number of digits to keep in coordinates
#'
#' @return A feature collection of linestrings with the mirror edges removed
#' @keywords internal
#' @importFrom sf st_length
#' @importFrom data.table as.data.table setDT :=
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' #This is an internal function, no example provided
remove_mirror_edges <- function(lines, keep_shortest = TRUE, digits = 3, verbose = TRUE){

    # step1 : get start and end coordinates of each line
    lines$tmpOID <- 1:nrow(lines)
    coords <- lines_extremities(lines)
    coords$spindex <- sp_char_index(st_coordinates(coords),digits)

    starts <- subset(coords, coords$pttype == "start")
    starts$geometry <- NULL
    ends <- subset(coords, coords$pttype == "end")
    ends$geometry <- NULL
    # merging with data.table
    dfline <- lines
    dfline$geometry <- NULL

    tempDT <- data.table::as.data.table(dfline)
    lines$startidx <- setDT(tempDT)[starts, on = "tmpOID", "startidx" := starts$spindex][]$startidx
    lines$endidx <- setDT(tempDT)[ends, on = "tmpOID", "endidx" := ends$spindex][]$endidx

    # create two sp_index (one start->end and end->start)
    lines$idx1 <- paste(lines$startidx,lines$endidx,sep=" - ")
    lines$idx2 <- paste(lines$endidx,lines$startidx,sep=" - ")

    # step2 : for each line, find if another has the same paire of points (or reversed)
    # if it is the case, flag the longest or shortest for a removal after
    if(verbose){
        pb <- txtProgressBar(min = 0, max = nrow(lines), style = 3)
    }
    to_remove <- sapply(1:nrow(lines), function(i){
        line <- lines[i,]
        idx1 <- line$idx1
        if(verbose){
            setTxtProgressBar(pb, i)
        }
        test <- (lines$idx1 == idx1 | lines$idx2 == idx1)
        if (sum(test)>1){
            sub <- subset(lines, test)
            lengths <- as.numeric(st_length(sub))
            # if all the lengths are equal
            if(length(unique(lengths))==1){
                if(sum(line$tmpOID > sub$tmpOID)==0){
                    return(TRUE)
                }else{
                    return(FALSE)
                }
            }
            # if some lengths are longer
            l <- as.numeric(st_length(line))
            longest <- max(lengths)
            if(keep_shortest){
                if(longest > l){
                    return(FALSE)
                }else{
                    return(TRUE)
                }
            }else{
                if(longest < l){
                    return(FALSE)
                }else{
                    return(TRUE)
                }
            }
        }else{
            return(FALSE)
        }
    })
    keeped <- subset(lines, to_remove == FALSE)
    return(keeped)
}


#' @title Simplify a network
#'
#' @description Simplify a network by applying two corrections: Healing edges and
#' Removing mirror edges (experimental).
#'
#' @details Healing is the operation to merge two connected linestring if the are
#' intersecting at one extremity and do not intersect any other linestring. It helps
#' to reduce the complexity of the network and thus can reduce calculation time.
#' Removing mirror edges is the operation to remove edges that have the same
#' extremities. If two edges start at the same point and end at the same point,
#' they do not add information in the network and one can be removed to simplify
#' the network. One can decide to keep the longest of the two edges or the shortest.
#' NOTE: the edge healing does not consider lines directions currently!
#'
#' @param lines A feature collection of linestrings
#' @param digits An integer indicating the number of digits to keep in coordinates
#' @param heal A boolean indicating if the healing operation must be performed
#' @param mirror A boolean indicating if the mirror edges must be removed
#' @param keep_shortest A boolean, if TRUE, then the shortest line is kept from
#' mirror edges. if FALSE, then the longest line is kept.
#' @param verbose A boolean indicating if messages and a progress bar should be displayed
#' @return A feature collection of linestrings
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' edited_lines <- simplify_network(mtl_network, digits = 3, verbose = FALSE)
#' }
simplify_network <- function(lines, digits = 3, heal = TRUE, mirror = TRUE, keep_shortest = TRUE, verbose = TRUE){ # nocov start

    if(heal){
        if(verbose){
            print("healing the connected edges...")
        }
        lines <- heal_edges(lines, digits, verbose)
    }

    if(mirror){
        if(verbose){
            print("removing mirror edges")
        }
        lines <- remove_mirror_edges(lines, keep_shortest, digits, verbose)
    }

    return(lines)
} # nocov end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Development ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# # # this is the old function kept if bugs are found
# split_at_vertices <- function(line, points, i, mindist) {
#     # extract coordinates
#     line_coords <- coordinates(line)[[1]][[1]]
#     original_distances <- sapply(1:nrow(line_coords),function(i){
#         if (i==0){
#             return(0)
#         }else{
#             return(sqrt(sum((line_coords[i,]-line_coords[i-1,])**2)))
#         }
#     })
#     line_coords <- cbind(line_coords, cumsum(original_distances),0)
#     tot_lengths <- max(line_coords[,3])
#     # calculate lengths
#     pt_coords <- coordinates(points)
#     lengths <- gProject(line, points)
#     pt_coords <- cbind(pt_coords, lengths)
#     pt_coords <- pt_coords[(lengths>mindist & lengths<(tot_lengths-mindist)),]
#     if(is.null(nrow(pt_coords))){
#         pt_coords <- c(pt_coords,1)
#     }else{
#         pt_coords <- cbind(pt_coords,1)
#     }
#
#     all_coords <- rbind(line_coords,pt_coords)
#     # reorder the coordinate matrix
#     ord_coords <- all_coords[order(all_coords[,3]), ]
#
#     #split on new coordinates
#     ruptidx <- unique(c(1,(1:nrow(ord_coords))[ord_coords[,4] == 1],nrow(ord_coords)))
#     final_coords <- lapply(1:(length(ruptidx)-1), function(j){
#         els <- ord_coords[ruptidx[[j]]:ruptidx[[j+1]],1:2]
#     })
#
#     return(final_coords)
# }

# # this is the previous function, kept for debug
# split_lines_at_vertex2 <- function(lines, points, nearest_lines_idx, mindist) {
#     new_lines_list <- lapply(1:nrow(lines), function(i) {
#         line <- lines[i, ]
#         testpts <- nearest_lines_idx == i
#         if (any(testpts)) {
#             okpts <- subset(points,testpts)
#             newline <- split_at_vertices(line, okpts, i, mindist)
#             return(newline)
#         } else {
#             sline <- list(coordinates(line)[[1]][[1]])
#             return(sline)
#         }
#     })
#     final_lines <- do.call(raster::spLines,unlist(new_lines_list, recursive = FALSE))
#     idxs <- do.call(c,lapply(1:length(new_lines_list), function(j){
#         el <- new_lines_list[[j]]
#         if (class(el) == "list"){
#             return(rep(j,times = length(el)))
#         }else{
#             return(j)
#         }
#     }))
#     final_lines <- SpatialLinesDataFrame(final_lines,
#                                          lines@data[idxs,],match.ID = FALSE)
#     raster::crs(final_lines) <- raster::crs(lines)
#     return(final_lines)
# }


#' @title Split lines at vertices in a feature collection of linestrings
#'
#' @description Split lines (feature collection of linestrings) at their nearest vertices
#' (feature collection of points), may fail if the line geometries are self intersecting.
#'
#' @param lines The feature collection of linestrings to split
#' @param points The feature collection of points to add to as vertex to the lines
#' @param nearest_lines_idx For each point, the index of the nearest line
#' @param mindist The minimum distance between one point and the extremity of
#'   the line to add the point as a vertex.
#' @return A feature collection of linestrings
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @examples
#' \donttest{
#' # reading the data
#' data(mtl_network)
#' data(bike_accidents)
#' # aggregating points within a 5 metres radius
#' bike_accidents$weight <- 1
#' agg_points <- aggregate_points(bike_accidents, 5)
#' mtl_network$LineID <- 1:nrow(mtl_network)
#' # snapping point to lines
#' snapped_points <- snapPointsToLines2(agg_points,
#'     mtl_network,
#'     "LineID"
#' )
#' # splitting lines
#' new_lines <- split_lines_at_vertex(mtl_network, snapped_points,
#'     snapped_points$nearest_line_id, 1)
#' }
split_lines_at_vertex <- function(lines, points, nearest_lines_idx, mindist) {

    coords <- st_coordinates(points)
    # step1 : remove points that are not far enough from the extremities
    exts <- lines_extremities(lines)

    coords2 <- st_coordinates(exts)

    #min_dists <- FNN::knnx.dist(data = coords2, query = coords,k = 1)
	min_dists <- dbscan::kNN(x = coords2, query = coords,k = 1)$dist
    coords <- subset(coords, min_dists > mindist)
    nearest_lines_idx <- nearest_lines_idx[min_dists > mindist]
    # no need to split here, great !
    if(nrow(coords) == 0){
        return(lines)
    }else{
        #lines_coords <- unlist(sp::coordinates(lines), recursive = FALSE)
        lines_coords <- lines_coordinates_as_list(lines)
        elements <- split_lines_at_points_cpp(coords, lines_coords, nearest_lines_idx, mindist)
        new_lines_list <- elements[[1]]
        idxs <- elements[[2]]

        #final_lines <- do.call(raster::spLines,new_lines_list)
        final_lines <- list_coordinates_as_lines(new_lines_list, crs = st_crs(lines))
        final_lines <- cbind(final_lines, st_drop_geometry(lines[idxs,]))
        return(final_lines)
    }
}




#' @title Cut lines at a specified distance
#'
#' @description Cut lines in a feature collection of linestrings at a specified distance from the
#' begining of the lines.
#'
#' @param lines The feature collection of linestrings to cut
#' @param dists A vector of distances, if only one value is given,
#' each line will be cut at that distance.
#' @return A feature collection of linestrings
#' @keywords internal
#' @examples
#' # This is an interal function, no example provided
cut_lines_at_distance <- function(lines, dists){

    # step 1 : create a list of coordinates
    coord_lists <- lines_coordinates_as_list(lines)

    new_coords <- cut_lines_at_distances_cpp(coord_lists, dists)

    final_lines <- list_coordinates_as_lines(new_coords, crs = st_crs(lines))
    final_lines <- cbind(final_lines, st_drop_geometry(lines))

    return(final_lines)
}
