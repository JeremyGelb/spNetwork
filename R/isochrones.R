
#' @title Isochrones calculation
#'
#' @description Calculate isochrones on a network
#'
#' @param lines A SpatialLinesDataFrame representing the edges of the network
#' @param dists A vector of distances to use to calculate for each starting
#' point to create the isochrones
#' @param start_points A SpatialPointsDataFrame representing the starting
#' points if the isochrones
#' @param mindist The minimum distance between two points. When two points are
#' too close, they might end up snapped at the same location on a line.
#' Default is 1.
#' @param weight The name of the column in lines to use a edge weight. If NULL,
#' the geographical length is used. Note that when lines are split create the
#' network, the weight column is recalculated proportionally to the new lines
#' length.
#' @param direction The name of the column giving informations about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions (undirected). The values of the column must be "FT" (From - To),
#' "TF" (To - From) or "Both".
#' @return A SpatialLinesDataFrame representing the isochrones with the
#' following columns
#' \itemize{
#'         \item point_id: the index of the point at the centre of the isochrone
#'         \item distance: the size of the isochrone
#' }
#' @importFrom sp coordinates SpatialPoints
#' @importFrom utils strcapture
#' @export
#' @examples
#' library(sp)
#' # creating a simple network
#'wkt_lines <- c(
#'  "LINESTRING (0.0 0.0, 5.0 0.0)",
#'  "LINESTRING (0.0 -5.0, 5.0 -5.0)",
#'  "LINESTRING (5.0 0.0, 5.0 5.0)",
#'  "LINESTRING (5.0 -5.0, 5.0 -10.0)",
#'  "LINESTRING (5.0 0.0, 5.0 -5.0)",
#'  "LINESTRING (5.0 0.0, 10.0 0.0)",
#'  "LINESTRING (5.0 -5.0, 10.0 -5.0)",
#'  "LINESTRING (10.0 0, 10.0 -5.0)",
#'  "LINESTRING (10.0 -10.0, 10.0 -5.0)",
#'  "LINESTRING (15.0 -5.0, 10.0 -5.0)",
#'  "LINESTRING (10.0 0.0, 15.0 0.0)",
#'  "LINESTRING (10.0 0.0, 10.0 5.0)")
#'
#'linesdf <- data.frame(wkt = wkt_lines,
#'                      id = paste("l",1:length(wkt_lines),sep=""))
#'
#'geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
#'  txt <- as.character(linesdf[i,]$wkt)
#'  geom <- rgeos::readWKT(txt,id=i)
#'  return(geom)
#'}))
#'
#'lines <- SpatialLinesDataFrame(geoms, linesdf, match.ID = FALSE)
#'
#'# and the definition of the starting point
#'start_points <- data.frame(x=c(5),
#'                           y=c(-2.5))
#'coordinates(start_points) <- cbind(start_points$x,start_points$y)
#'
#'# setting the directions
#'
#'lines$direction <- "Both"
#'lines[6,"direction"] <- "TF"
#'
#'isochrones <- calc_isochrones(lines, dists = c(10,12),
#'                              start_points = start_points,
#'                              direction = "direction")
calc_isochrones <- function(lines, dists, start_points, mindist = 1, weight = NULL, direction = NULL){

  # step1 : verifier que certains points ne sont pas a une distance
  # inferieure a mindist de d'autres points
  if(nrow(start_points) >1){
    xy <- sp::coordinates(start_points)
    min_dists <- FNN::get.knn(xy, k = 1)
    if(any(min_dists < mindist)){
      stop("some points are closer to each other than the mindist argument value before snapping")
    }
  }


  # step2 : snapper les points sur les lignes
  lines$OID <- 1:nrow(lines)
  snapped_points <- snapPointsToLines2(start_points, lines, "OID")

  # step3 : reverifier la proximite
  if(nrow(start_points) >1){
    xy <- sp::coordinates(snapped_points)
    min_dists <- FNN::get.knn(xy, k = 1)
    if(any(min_dists < mindist)){
      stop("some points are closer to each other than the mindist argument value after snapping them on lines")
    }
  }


  # step4 : splitting the lines with the points
  lines$tot_length <- rgeos::gLength(lines, byid = TRUE)
  new_lines <- split_lines_at_vertex(lines, snapped_points, snapped_points$nearest_line_id, mindist = 0)
  new_lines$length <- rgeos::gLength(new_lines, byid = TRUE)

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
    #new_lines <- lines_direction(new_lines, direction)
    #new_lines$b_direction <- ifelse(new_lines[[direction]] == "Both", 0,1)
    graph_result <- build_graph_directed(new_lines, digits = 2, line_weight = weight,
                                         direction = direction, attrs = TRUE)
  }


  # trouver pour chaque point de depart le noeud correspondant
  maxdist <- max(dists)
  xynodes <- sp::coordinates(graph_result$spvertices)
  xy_points <- sp::coordinates(snapped_points)

  start_nodes <- FNN::get.knnx(xynodes, xy_points, k=1)
  df_start <- data.frame(
    "ptOID" = 1:nrow(snapped_points),
    "node_idx" = start_nodes[[1]],
    "node_id" = graph_result$spvertices$id[start_nodes[[1]]],
    "node_name" = graph_result$spvertices$name[start_nodes[[1]]]
  )

  # creer les isochrones, soit des ensembles de multilignes
  all_multi_lignes <- lapply(1:nrow(df_start), function(i){
    row <- df_start[i,]

    # pour chaque point de depart calculer les distances a tous les autres points
    all_dists <- t(igraph::distances(graph_result$graph,
                                   v = row$node_idx,
                                   to = graph_result$spvertices$id,
                                   mode = "out"))
    dist_df <- data.frame(
      "node_id" = 1:nrow(all_dists),
      "dist" = as.vector(all_dists)
    )

    # on va maintenant iterer sur toutes les distances demandees
    all_lignes <- lapply(dists, function(d){
      # on retrouve les noeuds a utiliser
      ok_nodes <- subset(dist_df, dist_df$dist <= d)
      # et toutes les edges
      ok_edges <- subset(graph_result$spedges,
                         graph_result$spedges$start_oid %in% ok_nodes$node_id |
                           graph_result$spedges$end_oid %in% ok_nodes$node_id
                           )

      # on trouve maintenant les edges sur lesquelles on doit retravailler un peu
      df1 <- ok_edges
      df2 <- ok_nodes

      df1<- sp::merge(df1,df2, by.x = "start_oid",
                       by.y = "node_id", all.x = TRUE,
                       all.y = FALSE)

      df1$start_dist <- df1$dist
      df1$dist <- NULL


      df1 <- sp::merge(df1,df2, by.x = "end_oid",
                              by.y = "node_id", all.x = TRUE,
                              all.y = FALSE)

      df1$end_dist <- df1$dist
      df1$dist <- NULL

      # si on est en mode directed, il faut degager certaines edges
      if(is.null(direction) == FALSE){
        df1 <- subset(df1, (is.na(df1$start_dist) & df1$direction == "TF") == FALSE)
        df1 <- subset(df1, (is.na(df1$end_dist) & df1$direction == "FT") == FALSE)
      }

      # il reste encore a decouper les edges en question si necessaire
      test <- is.na(df1$start_dist)==FALSE & is.na(df1$end_dist)==FALSE
      no_cut <- subset(df1, test)

      # decoupage des lignes par le debut
      to_cut <- subset(df1, !test)
      to_cut$node_okid <- ifelse(is.na(to_cut$start_dist), to_cut$end_oid,
                                 to_cut$start_oid
                                 )
      to_cut$ok_dist <- ifelse(is.na(to_cut$start_dist), to_cut$end_dist,
                               to_cut$start_dist
      )

      # remise dans l'ordre des lignes si necessaire
      ext <- lines_extremities(to_cut)
      ok_nodes <- graph_result$spvertices[to_cut$node_okid,]
      test1 <- subset(ext, ext$pttype == "start")
      test2 <- subset(ext, ext$pttype == "end")

      d1 <- sqrt((test1$X - ok_nodes$x)**2 + (test1$Y - ok_nodes$y)**2)
      d2 <- sqrt((test2$X - ok_nodes$x)**2 + (test2$Y - ok_nodes$y)**2)
      to_keep <- subset(to_cut, d1 < d2)
      to_reverse <- subset(to_cut, d1 >= d2)

      # et decoupage final
      all_dists <- d - c(to_keep$ok_dist, to_reverse$ok_dist)
      all_cuts <- rbind(to_keep, reverse_lines(to_reverse))
      cut_lines <- cut_lines_at_distance(all_cuts, all_dists)

      # enregistrement
      ok_lines <- rbind(no_cut[c("end_oid","start_oid","edge_id","weight" )],
                        cut_lines[c("end_oid","start_oid","edge_id","weight" )])
      ok_lines$distance <- d
      ok_lines$point_id <- i
      ok_lines <- ok_lines[c("point_id", "distance")]
      return(ok_lines)

    })
    return(all_lignes)

  })

  all_multi_lignes <- do.call(rbind, unlist(all_multi_lignes))
  return(all_multi_lignes)

}
