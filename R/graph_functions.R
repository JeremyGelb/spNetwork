#' @title Network generation
#'
#' @description Generate an igraph object from a SpatialLinesDataFrame.
#'
#' @param lines A SpatialLinesDataFrame
#' @param digits The number of digits to keep from the coordinates
#' @param line_weight The name of a field that represent the cost to use a line
#' @param attrs A Boolean indicating if the original lines attributes
#' must be added to the graph lines
#' @return A list containing the following elements:
#' \itemize{
#'         \item graph: an igraph object that preserves the original lines
#'         characteristics
#'         \item linelist: the dataframe used to build the graph
#'         \item lines: the original SpatialLinesDataFrame
#'         \item spvertices: a SpatialPointsDataFrame representing the vertices
#'         of the graph
#'         \item digits : the number of digits kept for the coordinates
#' }
#' @importFrom sp coordinates SpatialPoints
#' @importFrom utils strcapture
#' @export
#' @examples
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
#' mtl_network$length <- rgeos::gLength(mtl_network, byid = TRUE)
#' graph_result <- build_graph(mtl_network, 2, "length", attrs = TRUE)
build_graph <- function(lines, digits, line_weight, attrs = FALSE) {
    # extracting lines coordinates lines_coords
    extremites <- lines_extremities(lines)
    start_coords <- extremites@data[extremites$pttype == "start", c("X","Y")]
    end_coords <- extremites@data[extremites$pttype == "end", c("X", "Y")]
    # extracting the coordinates of the starting and end points start_coords
    start <- sp_char_index(start_coords, digits)
    end <- sp_char_index(end_coords, digits)

    weights <- lines[[line_weight]]

    # building the line list
    linelist <- data.frame(start = start, end = end, weight = weights,
        graph_id = 1:nrow(lines), wkt= rgeos::writeWKT(lines,byid = TRUE))
    if (attrs) {
        linelist <- cbind(linelist, lines@data)
    }
    ## generating the graph
    graph <- igraph::graph_from_data_frame(linelist, directed = FALSE,
                                           vertices = NULL)


    ## building a spatial object for the vertices
    vertices <- igraph::V(graph)
    dfvertices <- data.frame(name = names(vertices), id = as.vector(vertices))
    dfvertices$name <- as.character(dfvertices$name)
    cols <- strcapture("(.*)_(.*)",dfvertices$name,data.frame(x = "", y = ""))
    dfvertices$x <- as.numeric(cols$x)
    dfvertices$y <- as.numeric(cols$y)
    points <- SpatialPoints(dfvertices[c("x", "y")])
    points <- SpatialPointsDataFrame(points, dfvertices)
    raster::crs(points) <- raster::crs(lines)

    ##building a spatial object for the lines
    edge_attrs <- igraph::get.edge.attribute(graph)
    edge_df <- data.frame(
      "edge_id" = as.numeric(igraph::E(graph)),
      "weight" = edge_attrs$weight
    )
    edge_df$wkt <- edge_attrs$wkt

    geoms <- do.call(rbind,lapply(1:nrow(edge_df),function(i){
      wkt <- edge_df[i,"wkt"]
      geom <- rgeos::readWKT(wkt,id=i)
      return(geom)
    }))

    spedges <- SpatialLinesDataFrame(geoms, edge_df,match.ID = FALSE)
    raster::crs(spedges) <- raster::crs(lines)
    vertex_df <- igraph::ends(graph,spedges$edge_id,names = FALSE)
    spedges$start_oid <- vertex_df[,1]
    spedges$end_oid <- vertex_df[,2]

    vertex_df <- igraph::ends(graph,linelist$graph_id,names = FALSE)
    linelist$start_oid <- vertex_df[,1]
    linelist$end_oid <- vertex_df[,2]

    return(list(graph = graph, linelist = linelist, lines = lines,
                spvertices = points, digits = digits, spedges = spedges))
}

#' @title Directed network generation
#'
#' @description Generate a directed igraph object from a SpatialLinesDataFrame.
#'
#' @param lines A SpatialLinesDataFrame
#' @param digits The number of digits to keep from the coordinates
#' @param line_weight The name of a field that represent the cost to use a line
#' @param attrs A boolean indicating if the original lines attributes
#' @param direction Indicate a field giving informations about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @return A list containing the following elements:
#' \itemize{
#'         \item graph: an igraph object that preserves the original lines
#'         characteristics
#'         \item linelist: the dataframe used to build the graph
#'         \item lines: the original SpatialLinesDataFrame
#'         \item spvertices: a SpatialPointsDataFrame representing the vertices
#'         of the graph
#'         \item digits : the number of digits kept for the coordinates
#' }
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils strcapture
#' @export
#' @examples
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
#' mtl_network$length <- rgeos::gLength(mtl_network, byid = TRUE)
#' mtl_network$direction <- "Both"
#' mtl_network[6, "direction"] <- "TF"
#' mtl_network_directed <- lines_direction(mtl_network, "direction")
#' graph_result <- build_graph_directed(lines = mtl_network_directed,
#'         digits = 2,
#'         line_weight = "length",
#'         direction = "direction",
#'         attrs = TRUE)
build_graph_directed <- function(lines, digits, line_weight, direction, attrs = FALSE) {

  # doubling the lines if needed
  all_lines <- lines_direction(lines, direction)
  dir <- ifelse(all_lines[[direction]] =="Both", 0,1)
  all_lines <- direct_lines(lines, dir)

  # extracting lines coordinates
  extremites <- lines_extremities(all_lines)
  start_coords <- extremites@data[extremites$pttype == "start", c("X","Y")]
  end_coords <- extremites@data[extremites$pttype == "end", c("X", "Y")]
  # extracting the coordinates of the starting and end points start_coords
  start <- sp_char_index(start_coords, digits)
  end <- sp_char_index(end_coords, digits)

  weights <- all_lines[[line_weight]]

  # building the line list
  linelist <- data.frame(start = start,
                         end = end,
                         weight = weights,
                         wkt= rgeos::writeWKT(all_lines,byid = TRUE),
                         graph_id = 1:nrow(all_lines))
  if (attrs) {
    linelist <- cbind(linelist, all_lines@data)
  }
  graph <- igraph::graph_from_data_frame(linelist, directed = TRUE, vertices = NULL)
  vertices <- igraph::V(graph)
  dfvertices <- data.frame(name = names(vertices), id = as.vector(vertices))
  dfvertices$name <- as.character(dfvertices$name)
  cols <- strcapture("(.*)_(.*)",dfvertices$name,data.frame(x = "", y = ""))
  dfvertices$x <- as.numeric(cols$x)
  dfvertices$y <- as.numeric(cols$y)
  points <- SpatialPoints(dfvertices[c("x", "y")])
  points <- SpatialPointsDataFrame(points, dfvertices)
  raster::crs(points) <- raster::crs(lines)

  ##building a spatial object for the lines
  edge_attrs <- igraph::get.edge.attribute(graph)
  edge_df <- data.frame(
    "edge_id" = as.numeric(igraph::E(graph)),
    "weight" = edge_attrs[[line_weight]],
    "direction" = edge_attrs[[direction]]
  )
  edge_df$wkt <- edge_attrs$wkt

  geoms <- do.call(rbind,lapply(1:nrow(edge_df),function(i){
    wkt <- edge_df[i,"wkt"]
    geom <- rgeos::readWKT(wkt,id=i)
    return(geom)
  }))

  spedges <- SpatialLinesDataFrame(geoms, edge_df,match.ID = FALSE)
  raster::crs(spedges) <- raster::crs(lines)
  vertex_df <- igraph::ends(graph,spedges$edge_id,names = FALSE)
  spedges$start_oid <- vertex_df[,1]
  spedges$end_oid <- vertex_df[,2]

  vertex_df <- igraph::ends(graph,linelist$graph_id,names = FALSE)
  linelist$start_oid <- vertex_df[,1]
  linelist$end_oid <- vertex_df[,2]


  return(list(graph = graph, linelist = linelist, lines = all_lines,
              spvertices = points, digits = digits, spedges = spedges))
}


#'@title Match nodes and points
#'
#'@description Function to match some points (SpatialPointsDataFrame) to the
#'  vertices of a graph.
#'
#'@param spvertices The spatial vertices of a graph (produced with build_graph)
#'@param points the SpatialPointsDataFrame to match
#'@param digits The number of digits to keep from the coordinates
#'@param tol the maximum distance between a point and a vertex
#'@return A vector of the corresponding vertices id, multiple points may belong
#'  to the same vertex
#'@importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#'@importFrom utils txtProgressBar setTxtProgressBar
#'@importFrom rgeos gIntersects gBuffer
#'@keywords internal
#' @examples
#' #This is an internal function, no example provided
find_vertices <- function(spvertices, points, digits, tol = 0.1) {
    # step1 : calculate the spatialnameindex
    coords <- coordinates(points)
    points$tempoid <- 1:nrow(points)
    points$spIndex <- sp_char_index(coords, digits)
    # step2 : check which points are already well matched
    matching <- data.table(points@data)
    B <- data.table(spvertices@data)
    matching[B, on = c("spIndex"="name"), names(B) := mget(paste0("i.", names(B)))]
    matching <- matching[, .SD[1], "tempoid"]

    if (any(is.na(matching$id)) == FALSE) {
        return(matching$id)
    } else {
        # so we have some points that need to be matched
        pb <- txtProgressBar(min = 0, max = nrow(matching), style = 3)
        matching <- SpatialPointsDataFrame(coords, matching)
        raster::crs(matching) <- raster::crs(points)
        # si on a des cas manquant, allons les chercher !
        values <- sapply(1:nrow(matching), function(i) {
            setTxtProgressBar(pb, i)
            pt <- matching[i, ]
            if (is.na(pt$id)) {
                test <- as.vector(gIntersects(spvertices, gBuffer(pt, width = tol),
                  prepared = TRUE, byid = TRUE))
                if (any(test) == FALSE) {
                  stop(paste("no node find at the demanded distance here ",
                    pt$spIndex, sep = ""))
                } else {
                  return(subset(spvertices, test)[["id"]][[1]])
                }
            } else {
                return(pt$id)
            }
        })
    }
}


#' @title Make a network directed
#'
#' @description Function to create complementary lines for a directed network.
#'
#' @param lines The original SpatialLinesDataFrame
#' @param direction A vector of integers. 0 indicates a bidirectional line and 1
#' an unidirectional line
#' @return A SpatialLinesDataFrame with some lines duplicated according to
#' direction
#' @importFrom sp coordinates Line Lines SpatialLinesDataFrame SpatialLines
#' @importFrom raster crs
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
direct_lines<-function(lines,direction){
  ##producing all the lines
  cnt <- 1
  allcoordinates <- coordinates(lines)
  listlines <- lapply(1:nrow(lines),function(i){
    coords <- allcoordinates[[i]][[1]]
    if(direction[[i]]==0){
      c1 <- coords
      c2 <- coords[nrow(coords):1,]
      L1 <- Lines(list(Line(c1)), ID = cnt)
      cnt <<- cnt+1
      L2 <- Lines(list(Line(c2)), ID = cnt)
      cnt <<- cnt+1
      return (list(L1,L2))
    }else{
      c1 <- coords
      L1 <- Lines(list(Line(c1)), ID = cnt)
      cnt <<- cnt+1
      return (list(L1))
    }
  })
  alllines <- unlist(listlines)
  splines <- SpatialLines(alllines)
  ##producing all the data
  oids <- lapply(1:nrow(lines),function(i){
    if(direction[[i]]==0){
      return(c(i,i))
    }else{
      return(c(i))
    }
  })
  oids <- do.call("c",oids)
  data <- lines@data[oids,]
  #combining the lines and the data
  df <- SpatialLinesDataFrame(splines,data,match.ID = FALSE)
  raster::crs(df) <- raster::crs(lines)
  return(df)
}



#' @title Plot graph
#'
#' @description Function to plot a graph (useful to check connectivity).
#'
#' @param graph A graph object (produced with build_graph)
#' @keywords internal
#' @importFrom utils strcapture
#' @examples
#' #This is an internal function, no example provided
plot_graph <- function(graph) {
    N <- data.frame(name = names(igraph::V(graph)), id = as.vector(igraph::V(graph)))
    cols <- strcapture("(.*)_(.*)",N$name,data.frame(x = "", y = ""))
    N$x <- as.numeric(cols$x)
    N$y <- as.numeric(cols$y)

    graphics::plot(graph, vertex.size = 0.01,
                   layout = as.matrix(N[c("x", "y")]),
                   vertex.label.cex = 0.1)

}

#' @title Topological error
#'
#' @description A utility function to find topological errors in a network.
#'
#' @param lines A SpatialLinesDataFrame representing the network
#' @param digits An integer indicating the number of digits to retain for
#'   coordinates
#' @param tol A float indicating the tolerance distance to identify a dangle
#'   node
#' @return A list with two elements. The first is a SpatialPointsDataFrame
#'   indicating for each node of the network to which component it belongs. The
#'   second is a SpatialPointsDataFrame with the dangle nodes of the network.
#' @importFrom rgeos gLength
#' @importFrom sp coordinates
#' @export
#' @examples
#' \donttest{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
#' topo_errors <- graph_checking(mtl_network, 2, 2)
#' }
graph_checking <- function(lines,digits, tol){

  ##step1 : adjusting the lines
  lines$length <- gLength(lines,byid = TRUE)
  lines <- subset(lines, lines$length>0)
  lines$oid <- 1:nrow(lines)

  lines <- simple_lines(lines)
  lines$length <- gLength(lines)

  ##step2 : building the graph
  graph_results <- build_graph(lines, digits, "length", attrs = FALSE)

  ##step3 : identify components
  parts <- igraph::components(graph_results$graph)
  graph_results$spvertices$component <- parts$membership

  ##step4 : identify dangle nodes
  graph_results$spvertices$degree <- igraph::degree(graph_results$graph)
  potential_error <- subset(graph_results$spvertices,
                            graph_results$spvertices$degree==1)

  node_tree <- build_quadtree(graph_results$spvertices)
  close_nodes_idx <- sapply(1:nrow(potential_error),function(i){
    row <- potential_error[i,]
    knn2 <- SearchTrees::knnLookup(node_tree,newdat = coordinates(row))
    return(knn2[[2]])
  })
  close_nodes <- graph_results$spvertices[close_nodes_idx,]
  XY1 <- coordinates(potential_error)
  XY2 <- coordinates(close_nodes)
  dist <- sqrt(rowSums((XY1 - XY2)**2))
  dangle_nodes <- subset(potential_error,dist<tol)
  return(list("dangle_nodes" = dangle_nodes,
              "vertex_components" = graph_results$spvertices))

}


#' @title Distance matrix with dupicated
#'
#' @description Function to Create a distance matrix when some vertices are duplicated.
#'
#' @param graph The Graph to use
#' @param start The vertices to use as starting points
#' @param end The vertices to use as ending points
#' @return A matrix with the distances between the vertices
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
dist_mat_dupl <- function(graph, start, end ){
  start <- as.numeric(start)
  end <- as.numeric(end)
  final_cols <- lapply(start, function(i){
    rows <- data.frame(start = rep(i,length(end)),
                       end = end)
    return(rows)
  })
  final_cols <- data.frame(do.call(rbind, final_cols))
  final_cols$oid <- paste(final_cols$start,final_cols$end, sep="_")
  ustart <- unique(start)
  udend <- unique(end)
  distmat <- igraph::distances(graph,ustart,udend, mode = "out")
  start_names <- row.names(distmat)
  start_codes <- as.numeric(igraph::V(graph)[start_names])
  end_names <- colnames(distmat)
  end_codes <- as.numeric(igraph::V(graph)[end_names])
  all_distances <- lapply(1:nrow(distmat),function(i){
    row <- as.numeric(distmat[i,])
    start_name <- start_codes[[i]]
    rows <- data.frame(start = rep(start_name,length(row)),
                       end = end_codes,
                       dist = row)
    return(rows)
  })
  all_distances <- do.call(rbind, all_distances)
  all_distances$oid <- paste(all_distances$start,all_distances$end, sep="_")
  ok_distances <- merge(final_cols, all_distances[c("oid","dist")], by = "oid", all.x = TRUE, all.y = FALSE)
  sorted_distances <- ok_distances[order(ok_distances$start, ok_distances$end),]
  # unmelting now
  jump <- length(end)
  rows <- lapply(1:length(start), function(i){
    if (i==1){
      i1 <- (i-1) * jump
      i2 <- i1 + 150
    }else {
      i1 <- ((i-1) * jump)
      i2 <- i1 + 150
      i1 <- i1 + 1
    }
    dists <- sorted_distances[i1:i2,"dist"]
    return(dists)
  })
  mat <- do.call(rbind,rows)
  row.names(mat) <- start
  colnames(mat) <- end
  return(mat)
}


