#' @title Network generation with igraph
#'
#' @description Generate an igraph object from a feature collection of linestrings
#'
#' @details This function can be used to generate an undirected graph object (igraph
#'   object). It uses the coordinates of the linestrings extremities to create
#'   the nodes of the graph. This is why the number of digits in the coordinates
#'   is important. Too high precision (high number of digits) might break some
#'   connections.
#'
#' @param lines A feature collection of lines
#' @param digits The number of digits to keep from the coordinates
#' @param line_weight The name of the column giving the weight of the lines
#' @param attrs A boolean indicating if the original lines' attributes should be
#'   stored in the final object
#' @return A list containing the following elements:
#' \itemize{
#'         \item graph: an igraph object;
#'         \item linelist: the dataframe used to build the graph;
#'         \item lines: the original feature collection of linestrings;
#'         \item spvertices: a feature collection of points representing the vertices
#'         of the graph;
#'         \item digits : the number of digits kept for the coordinates.
#' }
#' @importFrom sf st_as_text st_length st_geometry st_geometry<-
#' @importFrom utils strcapture
#' @export
#' @examples
#' data(mtl_network)
#' mtl_network$length <- as.numeric(sf::st_length(mtl_network))
#' graph_result <- build_graph(mtl_network, 2, "length", attrs = TRUE)
build_graph <- function(lines, digits, line_weight, attrs = FALSE) {
    # extracting lines coordinates lines_coords
    extremites <- lines_extremities(lines)
    start_coords <- st_drop_geometry(extremites[extremites$pttype == "start", c("X","Y")])
    end_coords <- st_drop_geometry(extremites[extremites$pttype == "end", c("X", "Y")])
    # extracting the coordinates of the starting and end points start_coords
    start <- sp_char_index(start_coords, digits)
    end <- sp_char_index(end_coords, digits)

    weights <- lines[[line_weight]]

    # building the line list
    linelist <- data.frame(start = start, end = end, weight = weights,
        graph_id = 1:nrow(lines), wkt= st_as_text(st_geometry(lines)))

    if (attrs) {
        linelist <- cbind(linelist, st_drop_geometry(lines))
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

    points <- sf_pts <- st_as_sf(
      x = dfvertices,
      coords = c("x","y"),
      crs = st_crs(lines)
    )
    points$x <- dfvertices$x
    points$y <- dfvertices$y

    ##building a spatial object for the lines
    edge_attrs <- igraph::get.edge.attribute(graph)
    edge_df <- data.frame(
      "edge_id" = as.numeric(igraph::E(graph)),
      "weight" = edge_attrs$weight
    )
    edge_df$wkt <- edge_attrs$wkt

    spedges <- st_as_sf(
      edge_df,
      wkt = "wkt",
      crs = st_crs(lines)
    )

    geoms <- spedges$wkt
    spedges <- st_drop_geometry(spedges)
    sf::st_geometry(spedges) <- geoms
    spedges$wkt <- edge_df$wkt

    vertex_df <- igraph::ends(graph,spedges$edge_id,names = FALSE)
    spedges$start_oid <- vertex_df[,1]
    spedges$end_oid <- vertex_df[,2]

    vertex_df <- igraph::ends(graph,linelist$graph_id,names = FALSE)
    linelist$start_oid <- vertex_df[,1]
    linelist$end_oid <- vertex_df[,2]

    return(list(graph = graph, linelist = linelist, lines = lines,
                spvertices = points, digits = digits, spedges = spedges))
}


#' @title Network generation with cppRouting
#'
#' @description Generate an cppRouting object from a feature collection of linestrings
#'
#' @details This function can be used to generate an undirected graph object (cppRouting
#'   object). It uses the coordinates of the linestrings extremities to create
#'   the nodes of the graph. This is why the number of digits in the coordinates
#'   is important. Too high precision (high number of digits) might break some
#'   connections.
#'
#' @param lines A feature collection of lines
#' @param digits The number of digits to keep from the coordinates
#' @param line_weight The name of the column giving the weight of the lines
#' @param attrs A boolean indicating if the original lines' attributes should be
#'   stored in the final object
#' @return A list containing the following elements:
#' \itemize{
#'         \item graph: a cppRouting object;
#'         \item linelist: the dataframe used to build the graph;
#'         \item lines: the original feature collection of linestrings;
#'         \item spvertices: a feature collection of points representing the vertices
#'         of the graph;
#'         \item digits : the number of digits kept for the coordinates.
#' }
#' @importFrom sf st_as_text st_length st_geometry st_geometry<-
#' @importFrom utils strcapture
#' @export
#' @keywords internal
#' @examples
#' \donttest{
#' data(mtl_network)
#' mtl_network$length <- as.numeric(sf::st_length(mtl_network))
#' graph_result <- build_graph_cppr(mtl_network, 2, "length", attrs = TRUE)
#' }
build_graph_cppr <- function(lines, digits, line_weight, attrs = FALSE, direction = NULL) {

  if(is.null(direction)){

    all_lines <- direct_lines(lines, rep(0, nrow(lines)))

  }else{
    # doubling the lines if needed
    all_lines <- lines_direction(lines, direction)
    dir <- ifelse(all_lines[[direction]] =="Both", 0,1)
    all_lines <- direct_lines(lines, dir)

  }

  lines <- all_lines

  # extracting lines coordinates lines_coords
  extremites <- lines_extremities(lines)
  extremites$name <- sp_char_index(st_coordinates(extremites), digits)

  # extracting the coordinates of the starting and end points start_coords
  start <- extremites[extremites$pttype == "start",]$name
  end <- extremites[extremites$pttype == "end",]$name

  weights <- lines[[line_weight]]

  # building the line list
  linelist <- data.frame(start = start, end = end, weight = weights,
                         graph_id = 1:nrow(lines), wkt= st_as_text(st_geometry(lines)))

  if (attrs) {
    linelist <- cbind(linelist, st_drop_geometry(lines))
  }

  # establishing the coordinates
  all_nodes <- data.frame(
    name = unique(c(linelist$start, linelist$end)))
  all_nodes$id <- 1:nrow(all_nodes)
  ids <- match(all_nodes$name,extremites$name)
  all_nodes$X <- extremites$X[ids]
  all_nodes$Y <- extremites$Y[ids]

  linelist$from <- all_nodes$id[match(linelist$start, all_nodes$name)]
  linelist$to <- all_nodes$id[match(linelist$end, all_nodes$name)]


  ## generating the graph
  graph <- cppRouting::makegraph(linelist[c('from','to','weight')], coords = all_nodes[c('id','X','Y')])

  # NOTE : in the dict, ref is the original id given, id is the code used in data.
  # WARNING, in graph$coords, the column id is not the id, but the ref


  ## building a spatial object for the vertices
  vertices <- graph$coords
  vertices$ref <- vertices$id
  vertices$id <- graph$dict$id[match(vertices$id, graph$dict$ref)]

  points <- sf_pts <- st_as_sf(
    x = vertices,
    coords = c("X","Y"),
    crs = st_crs(lines)
  )
  points$x <- vertices$X
  points$y <- vertices$Y


  ##building a spatial object for the lines
  ok_cols <- names(linelist)[!(names(linelist) %in% c('from', 'to', 'dist'))]
  edge_df <- cbind(graph$data, linelist[ok_cols])

  spedges <- st_as_sf(
    edge_df,
    wkt = "wkt",
    crs = st_crs(lines)
  )

  return(list(graph = graph,
              linelist = linelist,
              lines = lines,
              spvertices = points,
              digits = digits,
              spedges = spedges))
}



#' @title Directed network generation
#'
#' @description Generate a directed igraph object from a feature collection of linestrings
#'
#' @details This function can be used to generate a directed graph object (igraph
#'   object). It uses the coordinates of the linestrings extremities to create
#'   the nodes of the graph. This is why the number of digits in the coordinates
#'   is important. Too high precision (high number of digits) might break some
#'   connections. The column used to indicate directions can only have the
#'   following values: "FT" (From-To), "TF" (To-From) and "Both".
#'
#' @param lines A feature collection of linestrings
#' @param digits The number of digits to keep from the coordinates
#' @param line_weight The name of the column giving the weight of the lines
#' @param attrs A boolean indicating if the original lines' attributes should be
#'   stored in the final object
#' @param direction A column name indicating authorized travelling direction on
#'   lines. if NULL, then all lines can be used in both directions. Must be the
#'   name of a column otherwise. The values of the column must be "FT" (From -
#'   To), "TF" (To - From) or "Both"
#' @return A list containing the following elements:
#' \itemize{
#'         \item graph: an igraph object;
#'         \item linelist: the dataframe used to build the graph;
#'         \item lines: the original feature collection of lines;
#'         \item spvertices: a feature collection of points representing the vertices
#'         of the graph;
#'         \item digits : the number of digits kept for the coordinates.
#' }
#' @importFrom utils strcapture
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' mtl_network$length <- as.numeric(sf::st_length(mtl_network))
#' mtl_network$direction <- "Both"
#' mtl_network[6, "direction"] <- "TF"
#' mtl_network_directed <- lines_direction(mtl_network, "direction")
#' graph_result <- build_graph_directed(lines = mtl_network_directed,
#'         digits = 2,
#'         line_weight = "length",
#'         direction = "direction",
#'         attrs = TRUE)
#' }
build_graph_directed <- function(lines, digits, line_weight, direction, attrs = FALSE) {

  # doubling the lines if needed
  all_lines <- lines_direction(lines, direction)
  dir <- ifelse(all_lines[[direction]] =="Both", 0,1)
  all_lines <- direct_lines(lines, dir)

  # extracting lines coordinates
  extremites <- lines_extremities(all_lines)
  start_coords <- st_drop_geometry(extremites[extremites$pttype == "start", c("X","Y")])
  end_coords <- st_drop_geometry(extremites[extremites$pttype == "end", c("X", "Y")])
  # extracting the coordinates of the starting and end points start_coords
  start <- sp_char_index(start_coords, digits)
  end <- sp_char_index(end_coords, digits)

  weights <- all_lines[[line_weight]]

  # building the line list
  linelist <- data.frame(start = start,
                         end = end,
                         weight = weights,
                         wkt= st_as_text(st_geometry(all_lines)),
                         graph_id = 1:nrow(all_lines))
  if (attrs) {
    linelist <- cbind(linelist, st_drop_geometry(all_lines))
  }
  graph <- igraph::graph_from_data_frame(linelist, directed = TRUE, vertices = NULL)

  vertices <- igraph::V(graph)
  dfvertices <- data.frame(name = names(vertices), id = as.vector(vertices))
  dfvertices$name <- as.character(dfvertices$name)
  cols <- strcapture("(.*)_(.*)",dfvertices$name,data.frame(x = "", y = ""))
  dfvertices$x <- as.numeric(cols$x)
  dfvertices$y <- as.numeric(cols$y)

  points <- st_as_sf(
    dfvertices, coords = c("x","y"),
    crs = st_crs(lines)
  )
  points$x <- dfvertices$x
  points$y <- dfvertices$y

  ##building a spatial object for the lines
  edge_attrs <- igraph::get.edge.attribute(graph)
  edge_df <- data.frame(
    "edge_id" = as.numeric(igraph::E(graph)),
    "weight" = edge_attrs[[line_weight]],
    "direction" = edge_attrs[[direction]]
  )
  edge_df$wkt <- edge_attrs$wkt

  spedges <- st_as_sf(edge_df,
                      wkt = "wkt",
                      crs = st_crs(lines))
  geoms <- spedges$wkt
  spedges <- st_drop_geometry(spedges)
  sf::st_geometry(spedges) <- geoms
  spedges$wkt <- edge_df$wkt

  vertex_df <- igraph::ends(graph,spedges$edge_id,names = FALSE)
  spedges$start_oid <- vertex_df[,1]
  spedges$end_oid <- vertex_df[,2]

  vertex_df <- igraph::ends(graph,linelist$graph_id,names = FALSE)
  linelist$start_oid <- vertex_df[,1]
  linelist$end_oid <- vertex_df[,2]


  return(list(graph = graph, linelist = linelist, lines = all_lines,
              spvertices = points, digits = digits, spedges = spedges))
}




#' @title Make a network directed
#'
#' @description Function to create complementary lines for a directed network.
#'
#' @param lines The original feature collection of linestrings
#' @param direction A vector of integers. 0 indicates a bidirectional line and 1
#' an unidirectional line
#' @return A feature collection of linestrings with some lines duplicated according to
#' direction
#' @importFrom sf st_reverse
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
direct_lines <- function(lines,direction){
  ##producing all the lines
  lines$jgtmpid <- 1:nrow(lines)
  to_reverse <- lines[direction == 0,]
  to_keep <- lines[direction == 1,]

  reversed <- st_reverse(to_reverse)
  final_lines <- rbind(to_keep, to_reverse, reversed)
  final_lines <- final_lines[order(final_lines$jgtmpid),]

  return(final_lines)
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
plot_graph <- function(graph) { # nocov start
    N <- data.frame(name = names(igraph::V(graph)), id = as.vector(igraph::V(graph)))
    cols <- strcapture("(.*)_(.*)",N$name,data.frame(x = "", y = ""))
    N$x <- as.numeric(cols$x)
    N$y <- as.numeric(cols$y)

    graphics::plot(graph, vertex.size = 0.01,
                   layout = as.matrix(N[c("x", "y")]),
                   vertex.label.cex = 0.1)
}# nocov end

#' @title Topological error
#'
#' @description A utility function to find topological errors in a network.
#'
#' @details This function can be used to check for three common problems in
#'   networks: disconnected components, dangle nodes and close nodes. When a
#'   network has disconnected components, this means that several unconnected
#'   graphs are composing the overall network. This can be caused by topological
#'   errors in the dataset. Dangle nodes are nodes connected to only one other
#'   node. This type of node can be normal at the border of a network, but can
#'   also be caused by topological errors. Close nodes are nodes that are not
#'   coincident, but so close that they probably should be coincident.
#'
#' @param lines A feature collection of linestrings representing the network
#' @param digits An integer indicating the number of digits to retain for
#'   coordinates
#' @param max_search The maximum number of nearest neighbour to search to find
#'   close_nodes
#' @param tol The minimum distance expected between two nodes. If two nodes are
#'   closer, they are returned in the result of the function.
#' @return A list with three elements. The first is a feature collection of points
#'   indicating for each node of the network to which component it belongs. The
#'   second is a feature collection of points with nodes that are too close one of
#'   each other. The third is a feature collection of points with the dangle nodes of
#'   the network.
#' @export
#' @examples
#' \donttest{
#' data(mtl_netowrk)
#' topo_errors <- graph_checking(mtl_network, 2)
#' }
graph_checking <- function(lines,digits, max_search = 5, tol = 0.1){

  ##step1 : adjusting the lines
  lines$length <- as.numeric(st_length(lines))
  lines <- subset(lines, lines$length>0)
  lines$oid <- 1:nrow(lines)

  lines <- simple_lines(lines)
  lines$length <- as.numeric(st_length(lines))

  ##step2 : building the graph
  graph_results <- build_graph(lines, digits, "length", attrs = FALSE)

  ##step3 : identify components
  parts <- igraph::components(graph_results$graph)
  graph_results$spvertices$component <- parts$membership

  ##step4 : identify dangle nodes
  graph_results$spvertices$degree <- igraph::degree(graph_results$graph)
  dangle_nodes <- subset(graph_results$spvertices,
                            graph_results$spvertices$degree==1)

  ##step5 : find nodes that are closer to a tolerance
  xy_nodes <- st_coordinates(graph_results$spvertices)
  #close_dists <- FNN::knn.dist(xy_nodes, k = max_search)
  close_dists <- dbscan::kNN(xy_nodes, k = max_search)$dist
  is_error <- apply(t(close_dists),MARGIN = 2, min) <= tol
  close_nodes <- subset(graph_results$spvertices, is_error)


  return(list("dangle_nodes" = dangle_nodes,
              "close_nodes" = close_nodes,
              "vertex_components" = graph_results$spvertices))

}


#' @title Distance matrix with dupicated
#'
#' @description Function to Create a distance matrix when some vertices are duplicated.
#'
#' @param graph The Graph to use
#' @param start The vertices to use as starting points
#' @param end The vertices to use as ending points
#' @param ... parameters passed to the function igraph::distances
#' @return A matrix with the distances between the vertices
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
dist_mat_dupl <- function(graph, start, end, ...){ # nocov start


  ## we start by calculating the distance matrix with unique values

  dist_mat <- igraph::distances(graph,v = unique(start),
                                to = unique(end), ...)

  ## we must expand the matrix !
  L1 <- length(unique(start)) + length(unique(end))
  L2 <- length(start) + length(end)

  if(L1 != L2){

    pos_row <- match(start, unique(start))
    pos_col <- match(end, unique(end))
    dist_mat2 <- t(sapply(pos_row, function(ii){
      values <- dist_mat[ii,pos_col]
      return(values)
    }))
    dist_mat <- dist_mat2
  }

  return(dist_mat)

} # nocov end




#' @title Split graph components
#'
#' @description Function to split the results of build_graph and build_graph_directed
#' into their sub components
#'
#' @param graph_result A list typically obtained from the function build_graph or build_graph_directed
#' @return A list of lists, the graph_result split for each graph component
#' @export
#' @examples
#' data(mtl_network)
#' mtl_network$length <- as.numeric(sf::st_length(mtl_network))
#' graph_result <- build_graph(mtl_network, 2, "length", attrs = TRUE)
#' sub_elements <- split_graph_components(graph_result)
split_graph_components <- function(graph_result){ # nocov start

  # identifying the components of the graph
  comps <- igraph::components(graph_result$graph)

  # if we have only one component, we return it
  if(comps$no == 1){
    return(graph_result)
  }

  vals <- unique(comps$membership)
  graph_result$spvertices$comp <- comps$membership

  elements <- lapply(vals, function(val){

    # finding the elements in spvertices
    spvertices <- subset(graph_result$spvertices, graph_result$spvertices$comp == val)

    # finding the elements in spedges
    spedges <- subset(graph_result$spedges, graph_result$spedges$start_oid %in% spvertices$id |
                        graph_result$spedges$end_oid %in% spvertices$id)

    # and their corresponding elements in linelist
    linelist <- subset(graph_result$linelist, graph_result$linelist$graph_id %in% spedges$edge_id)

    #finding the subgraph
    graph <- igraph::induced_subgraph (graph_result$graph,
                              igraph::V(graph_result$graph)[comps$membership==val])
    #merging everything
    element <- list(
      "graph" = graph,
      "spvertices" = spvertices,
      "spedges" = spedges,
      "linelist" = linelist,
      "digits" = graph_result$digits
    )

  })

  return(elements)
}# nocov end
