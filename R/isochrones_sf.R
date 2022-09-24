
#' @title Isochrones calculation
#'
#' @description Calculate isochrones on a network
#'
#' @details An isochrone is the set of reachable lines around a node in a network within
#' a specified distance (or time). This function perform dynamic segmentation to return the
#' part of the edges reached and not only the fully covered edges. Several start points and
#' several distances can be given. The network can also be directed. The lines returned
#' by the function are the most accurate representation of the isochrones. However, if
#' polygons are required for mapping, the vignette "Calculating isochrones" shows
#' how to create smooth polygons from the returned sets of lines.
#'
#' @param lines A feature collection of lines representing the edges of the network
#' @param dists A vector of the size of the desired isochrones
#' @param start_points A feature collection of points representing the starting
#' points if the isochrones
#' @param donught A boolean indicating if the returned lines must overlap for
#' each distance (FALSE, default) or if the lines must be cut between each
#' distance step (TRUE).
#' @param mindist The minimum distance between two points. When two points are
#' too close, they might end up snapped at the same location on a line.
#' Default is 1.
#' @param weight The name of the column in lines to use an edge weight. If NULL,
#' the geographical length is used. Note that if lines are split during the
#' network creation, the weight column is recalculated proportionally to the new lines
#' length.
#' @param direction The name of the column indicating authorized
#' travelling direction on lines. if NULL, then all lines can be used in both
#' directions (undirected). The values of the column must be "FT" (From - To),
#' "TF" (To - From) or "Both".
#' @return A feature collection of lines representing the isochrones with the
#' following columns
#' \itemize{
#'         \item point_id: the index of the point at the centre of the isochrone;
#'         \item distance: the size of the isochrone
#' }
#' @importFrom utils strcapture
#' @export
#' @examples
#' library(sf)
#' # creating a simple network
#' wkt_lines <- c(
#'   "LINESTRING (0.0 0.0, 5.0 0.0)",
#'   "LINESTRING (0.0 -5.0, 5.0 -5.0)",
#'   "LINESTRING (5.0 0.0, 5.0 5.0)",
#'   "LINESTRING (5.0 -5.0, 5.0 -10.0)",
#'   "LINESTRING (5.0 0.0, 5.0 -5.0)",
#'   "LINESTRING (5.0 0.0, 10.0 0.0)",
#'   "LINESTRING (5.0 -5.0, 10.0 -5.0)",
#'   "LINESTRING (10.0 0, 10.0 -5.0)",
#'   "LINESTRING (10.0 -10.0, 10.0 -5.0)",
#'   "LINESTRING (15.0 -5.0, 10.0 -5.0)",
#'   "LINESTRING (10.0 0.0, 15.0 0.0)",
#'   "LINESTRING (10.0 0.0, 10.0 5.0)")
#'
#' linesdf <- data.frame(wkt = wkt_lines,
#'                       id = paste("l",1:length(wkt_lines),sep=""))
#'
#' lines <- st_as_sf(linesdf, wkt = "wkt")
#'
#' # and the definition of the starting point
#' start_points <- data.frame(x=c(5),
#'                            y=c(-2.5))
#' start_points <- st_as_sf(start_points, coords = c("x","y"))
#'
#' # setting the directions
#'
#' lines$direction <- "Both"
#' lines[6,"direction"] <- "TF"
#'
#' isochrones <- calc_isochrones(lines,dists = c(10,12),
#'                               donught = TRUE,
#'                               start_points = start_points,
#'                               direction = "direction")
calc_isochrones <- function(lines, dists, start_points, donught = FALSE, mindist = 1, weight = NULL, direction = NULL){

  # step1 : Check that some points are not too close to each other
  # before snapping
  if(nrow(start_points) >1){
    xy <- st_coordinates(start_points)
    #min_dists <- FNN::get.knn(xy, k = 1)
	min_dists <- dbscan::kNN(xy, k = 1)
    if(any(min_dists$dist < mindist)){
      stop("some points are closer to each other than the mindist argument value before snapping")
    }
  }

  # step2 : snapping hte points on the lines
  lines$OID <- 1:nrow(lines)
  snapped_points <- snapPointsToLines2(start_points, lines, "OID")

  # step3 : check the points closeness again
  if(nrow(start_points) >1){
    xy <- st_coordinates(snapped_points)
    #min_dists <- FNN::get.knn(xy, k = 1)
	min_dists <- dbscan::kNN(xy, k = 1)
    if(any(min_dists$dist < mindist)){
      stop("some points are closer to each other than the mindist argument value after snapping them on lines")
    }
  }


  # step4 : splitting the lines with the points
  lines$tot_length <- as.numeric(st_length(lines))
  new_lines <- split_lines_at_vertex(lines, snapped_points, snapped_points$nearest_line_id, mindist = 0)
  new_lines$length <-  as.numeric(st_length(new_lines))

  if(is.null(weight) == FALSE){
    new_lines[[weight]] <- new_lines[[weight]] * (new_lines$length / new_lines$tot_length)
  }else{
    weight <- "length"
  }

  # step5 : building the graph
  if (is.null(direction)){
    graph_result <- build_graph(lines = new_lines, line_weight = weight, digits = 2, attrs = TRUE)
  } else{
    vals <- unique(new_lines[[direction]])
    u <- union(vals, c("FT","TF","Both"))
    if(length(u) > 3){
      stop("when indicating line direction, accepted values are TF, FT and Both")
    }
    graph_result <- build_graph_directed(new_lines, digits = 2, line_weight = weight,
                                         direction = direction, attrs = TRUE)
  }


  # finding for each start points its corresponding node in the graph
  maxdist <- max(dists)
  xynodes <- st_coordinates(graph_result$spvertices)
  xy_points <- st_coordinates(snapped_points)

  start_nodes <- dbscan::kNN(xynodes, query = xy_points, k=1)
  df_start <- data.frame(
    "ptOID" = 1:nrow(snapped_points),
    "node_idx" = as.vector(start_nodes$id),
    "node_id" = graph_result$spvertices$id[start_nodes$id],
    "node_name" = graph_result$spvertices$name[start_nodes$id]
  )


  # tm_shape(lines) + tm_lines("black") + tm_shape(start_points) + tm_dots("red", size = 0.5)

  # creating the isochrones: sets of linestrings
  all_multi_lignes <- lapply(1:nrow(df_start), function(i){
    row <- df_start[i,]
    # pour chaque point de depart calculer les distances a tous les autres points
    all_dists <- t(igraph::distances(graph_result$graph,
                                   to = row$node_idx,
                                   v = graph_result$spvertices$id,
                                   mode = "out"))
    dist_df <- data.frame(
      "node_id" = graph_result$spvertices$id,
      "dist" = as.vector(all_dists)
    )

    # on va maintenant iterer sur toutes les distances demandees
    dists2 <- c(0,dists[1:(length(dists)-1)])
    all_lignes <- lapply(1:length(dists), function(j){
      d <- dists[[j]]
      dd <- dists2[[j]]

      # on retrouve tous les noeuds a utiliser
      ok_nodes <- subset(dist_df, dist_df$dist <= d)
      # et toutes les edges
      ok_edges <- subset(graph_result$spedges,
                         graph_result$spedges$start_oid %in% ok_nodes$node_id |
                           graph_result$spedges$end_oid %in% ok_nodes$node_id
                           )

      # on trouve maintenant les edges sur lesquelles on doit retravailler un peu
      df1 <- ok_edges
      df2 <- ok_nodes

      df1 <- merge(df1,df2, by.x = "start_oid",
                       by.y = "node_id", all.x = TRUE,
                       all.y = FALSE)

      df1$start_dist <- df1$dist
      df1$dist <- NULL


      df1 <- merge(df1,df2, by.x = "end_oid",
                              by.y = "node_id", all.x = TRUE,
                              all.y = FALSE)

      df1$end_dist <- df1$dist
      df1$dist <- NULL

      if(donught){
        df1 <- subset(df1,
                      df1$start_dist >= dd & df1$start_dist <= d |
                      df1$end_dist >= dd & df1$end_dist <= d |
                      is.na(df1$start_dist) |   is.na(df1$end_dist)
                        )
      }
      # tm_shape(lines) + tm_lines("black") + tm_shape(start_points) + tm_dots("red", size = 0.5) + tm_shape(df1) + tm_lines("blue")

      # If we are in a directed graph, some edges need to be removed here
      if(is.null(direction) == FALSE){
        df1 <- subset(df1, (is.na(df1$end_dist) & df1$direction == "TF") == FALSE)
        df1 <- subset(df1, (is.na(df1$start_dist) & df1$direction == "FT") == FALSE)
      }
      ok_lines <- trim_lines_at(df1, graph_result, d, dd, i, donught)
      # # we now have to cut the remaining edges
      # test <- is.na(df1$start_dist)==FALSE & is.na(df1$end_dist)==FALSE
      # no_cut <- subset(df1, test)
      #
      # # cutting the lines by the start
      # to_cut <- subset(df1, !test)
      # to_cut$node_okid <- ifelse(is.na(to_cut$start_dist), to_cut$end_oid,
      #                            to_cut$start_oid
      #                            )
      # to_cut$ok_dist <- ifelse(is.na(to_cut$start_dist), to_cut$end_dist,
      #                          to_cut$start_dist
      # )
      #
      # # reordering line if required
      # ext <- lines_extremities(to_cut)
      # ok_nodes <- graph_result$spvertices[to_cut$node_okid,]
      # test1 <- subset(ext, ext$pttype == "start")
      # test2 <- subset(ext, ext$pttype == "end")
      #
      # d1 <- sqrt((test1$X - ok_nodes$x)**2 + (test1$Y - ok_nodes$y)**2)
      # d2 <- sqrt((test2$X - ok_nodes$x)**2 + (test2$Y - ok_nodes$y)**2)
      # to_keep <- subset(to_cut, d1 < d2)
      # to_reverse <- subset(to_cut, d1 >= d2)
      #
      # # and final cut
      # all_dists <- d - c(to_keep$ok_dist, to_reverse$ok_dist)
      # if(nrow(to_reverse) > 0){
      #   all_cuts <- rbind(to_keep, reverse_lines(to_reverse))
      # }else{
      #   all_cuts <- to_keep
      # }
      #
      # cut_lines <- cut_lines_at_distance(all_cuts, all_dists)
      #
      # # saving
      # ok_lines <- rbind(no_cut[c("end_oid","start_oid","edge_id","weight" )],
      #                   cut_lines[c("end_oid","start_oid","edge_id","weight" )])
      # ok_lines$distance <- d
      # ok_lines$point_id <- i
      # ok_lines <- ok_lines[c("point_id", "distance")]
      return(ok_lines)

    })
    return(all_lignes)

  })

  all_multi_lignes <- do.call(rbind, unlist(all_multi_lignes, recursive = FALSE))
  return(all_multi_lignes)

}

#' @title Helper for isochrones lines cutting
#'
#' @description last operation for isochrone calculation, cutting the lines at
#' their begining and ending. This is a worker function for calc_isochrones.
#'
#' @param df1 A features collection of linestrings with some specific fields.
#' @param graph_result A list produced by the functions build_graph_directed or build_graph.
#' @param d the end distance of this isochrones.
#' @param dd the start distance of this isochrones.
#' @param i the actual iteration.
#' @param donught A boolean indicating if the returned isochrone will be plained or a donught.
#'
#' @return A feature collection of lines
#'
#' @keywords internal
trim_lines_at <- function(df1, graph_result, d, dd, i, donught){

  # first, we need to determine which edges will not need to be cut
  if(donught){
    test <- is.na(df1$start_dist) == FALSE &
      is.na(df1$end_dist) == FALSE &
      df1$start_dist >= dd &
      df1$end_dist >= dd
  }else{
    test <- is.na(df1$start_dist)==FALSE & is.na(df1$end_dist)==FALSE
  }

  if(sum(!test) > 0){
    no_cut <- subset(df1, test)
    to_cut <- subset(df1, !test)
    # here I need to know if the start point and end point of each line match with
    # the start node and end nodes

    # first, I find the coordinates of the nodes in the graph
    to_cut$start_nodeX <- graph_result$spvertices[to_cut$start_oid,]$x
    to_cut$end_nodeX <- graph_result$spvertices[to_cut$end_oid,]$x
    to_cut$start_nodeY <- graph_result$spvertices[to_cut$start_oid,]$y
    to_cut$end_nodeY <- graph_result$spvertices[to_cut$end_oid,]$y

    # second, I find the coordinates of the extremities of the lines
    ext <- lines_extremities(to_cut)
    coords <- st_coordinates(ext)
    start_coords <- subset(coords, ext$pttype == "start")

    # third, I calculate the distances between the nodes and the start point
    dist1 <- sqrt(((start_coords[,1] - to_cut$start_nodeX) ** 2) + ((start_coords[,2] - to_cut$start_nodeY) ** 2))
    dist2 <- sqrt(((start_coords[,1] - to_cut$end_nodeX) ** 2) + ((start_coords[,2] - to_cut$end_nodeY) ** 2))

    # If the distance to start_node is not the smallest, then I must reverse the geometry
    # it will insure that for each line, start_node is also the first vertex of the line
    test <- dist1 < dist2
    to_cut2 <- rbind(
      subset(to_cut, test),
      reverse_lines(subset(to_cut, !test))
    )

    # on va passer le tout a une fonction en c++ qui va iterer et regler son cas a chaque ligne
    list_of_lines <- lines_coordinates_as_list(to_cut2)
    start_dists <- ifelse(is.na(to_cut2$start_dist),-1,to_cut2$start_dist)
    end_dists <- ifelse(is.na(to_cut2$end_dist),-1,to_cut2$end_dist)

    cut_coords <- trim_lines_for_isos_cpp(list_of_lines, start_dists, end_dists, donught, d, dd)
    new_lines <- list_coordinates_as_lines(cut_coords, st_crs(to_cut2))
    to_cut2$geometry <- new_lines$geometry


    # saving
    ok_lines <- rbind(no_cut[c("end_oid","start_oid","edge_id","weight" )],
                      to_cut2[c("end_oid","start_oid","edge_id","weight" )])
  }else{
    ok_lines <- no_cut[c("end_oid","start_oid","edge_id","weight" )]
  }


  ok_lines$distance <- d
  ok_lines$point_id <- i
  ok_lines <- ok_lines[c("point_id", "distance")]
  return(ok_lines)
}



# trim_lines_at <- function(df1, graph_result, d, dd, i, donught){
#
#   # we now have to cut the remaining edges
#   if(donught){
#     test <- is.na(df1$start_dist)==FALSE &
#       is.na(df1$end_dist)==FALSE &
#       df1$start_dist >= dd &
#       df1$end_dist >= dd
#   }else{
#     test <- is.na(df1$start_dist)==FALSE & is.na(df1$end_dist)==FALSE
#   }
#
#   no_cut <- subset(df1, test)
#   to_cut <- subset(df1, !test)
#
#   if(donught){
#     ok_first_cut <- is.na(to_cut$start_dist) | is.na(to_cut$end_dist)
#   }else{
#
#   }
#
#   # cutting the lines by the start
#   to_cut$node_okid <- ifelse(is.na(to_cut$start_dist), to_cut$end_oid,
#                              to_cut$start_oid
#   )
#   to_cut$ok_dist <- ifelse(is.na(to_cut$start_dist), to_cut$end_dist,
#                            to_cut$start_dist
#   )
#
#
#   # reordering line if required
#   ext <- lines_extremities(to_cut)
#   ok_nodes <- graph_result$spvertices[to_cut$node_okid,]
#   test1 <- subset(ext, ext$pttype == "start")
#   test2 <- subset(ext, ext$pttype == "end")
#
#   d1 <- sqrt((test1$X - ok_nodes$x)**2 + (test1$Y - ok_nodes$y)**2)
#   d2 <- sqrt((test2$X - ok_nodes$x)**2 + (test2$Y - ok_nodes$y)**2)
#   to_keep <- subset(to_cut, d1 < d2)
#   to_reverse <- subset(to_cut, d1 >= d2)
#
#   # and final cut
#   all_dists <- d - c(to_keep$ok_dist, to_reverse$ok_dist)
#   if(nrow(to_reverse) > 0){
#     all_cuts <- rbind(to_keep, reverse_lines(to_reverse))
#   }else{
#     all_cuts <- to_keep
#   }
#
#   cut_lines <- cut_lines_at_distance(all_cuts, all_dists)
#
#   # HERE : TODO THINK ABOUT HOW TO ALSO CUT THE LINES AT THEIR BEGINING IF REQUIRED !
#   if(donught){
#     # HERE : TODO THINK ABOUT HOW TO ALSO CUT THE LINES AT THEIR BEGINING IF REQUIRED !
#     start_cuts <- dd - cut_lines$ok_dist
#     need_2nd_cut <- start_cuts > 0
#     if(sum(need_2nd_cut) > 0){
#       part1 <- subset(cut_lines, need_2nd_cut)
#       part2 <- subset(cut_lines, !need_2nd_cut)
#
#       cut_lines2 <-  cut_lines_at_distance(reverse_lines(part1),
#                                            as.vector(st_length(part1)) - start_cuts[need_2nd_cut]
#       )
#       cut_lines2$lineID.1 <- NULL
#       cut_lines <- rbind(cut_lines2,part2)
#     }
#   }
#
#   # saving
#   ok_lines <- rbind(no_cut[c("end_oid","start_oid","edge_id","weight" )],
#                     cut_lines[c("end_oid","start_oid","edge_id","weight" )])
#   ok_lines$distance <- d
#   ok_lines$point_id <- i
#   ok_lines <- ok_lines[c("point_id", "distance")]
#   return(ok_lines)
# }
