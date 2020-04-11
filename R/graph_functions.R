#' generate an igraph object from a SpatialLinesDataFrame
#'
#' @param lines A SpatialLinesDataFrame
#' @param digits The number of digits to keep from the coordinates
#' @param line_weight The name of a field that represent the cost to use a line
#' @param attrs A boolean indicating if the original lines attributes
#' must be added to the graph lines
#' @return A list containing the folowing elements:
#' \itemize{
#'         \item graph: an igraph object that preserves the original lines
#'         characteristics
#'         \item linelist: the dataframe used to build the graph
#'         \item lines: the original SpatialLinesDataFrame
#'         \item spvertices: a SpatialPointsDataFrame representing the vertices
#'         of the graph
#'         \item digits : the number of digits keeped for the coordinates
#' }
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @examples
#' #This is an internal function, no example provided
build_graph <- function(lines, digits, line_weight, attrs = FALSE) {
    # extracting lines coordinates lines_coords <-
    extremites <- lines_extremities(lines)
    start_coords <- extremites@data[extremites$pttype == "start", c("X","Y")]
    end_coords <- extremites@data[extremites$pttype == "end", c("X", "Y")]
    # extracting the coordinates of the starting and end points start_coords
    start <- sp_char_index(start_coords, digits)
    end <- sp_char_index(end_coords, digits)

    weights <- lines[[line_weight]]

    # building the line list
    linelist <- data.frame(start = start, end = end, weight = weights,
        graph_id = 1:nrow(lines))
    if (attrs) {
        linelist <- cbind(linelist, lines@data)
    }
    graph <- igraph::graph_from_data_frame(linelist, directed = FALSE, vertices = NULL)
    vertices <- igraph::V(graph)
    dfvertices <- data.frame(name = names(vertices), id = as.vector(vertices))
    dfvertices$name <- as.character(dfvertices$name)
    dfvertices <- tidyr::separate(dfvertices, col = "name", into = c("x",
        "y"), sep = "_", remove = F)
    dfvertices$x <- as.numeric(dfvertices$x)
    dfvertices$y <- as.numeric(dfvertices$y)
    points <- SpatialPoints(dfvertices[c("x", "y")])
    points <- SpatialPointsDataFrame(points, dfvertices)
    raster::crs(points) <- raster::crs(lines)
    return(list(graph = graph, linelist = linelist, lines = lines,
                spvertices = points, digits = digits))
}

#' generate an igraph object from a SpatialLinesDataFrame
#'
#' @param lines A SpatialLinesDataFrame
#' @param digits The number of digits to keep from the coordinates
#' @param line_weight The name of a field that represent the cost to use a line
#' @param attrs A boolean indicating if the original lines attributes
#' @param direction A vector of integers. 0 indicates a bidirectional line and 1
#' an unidirectional line
#' must be added to the graph lines
#' @return A list containing the folowing elements:
#' \itemize{
#'         \item graph: an igraph object that preserves the original lines
#'         characteristics
#'         \item linelist: the dataframe used to build the graph
#'         \item lines: the original SpatialLinesDataFrame
#'         \item spvertices: a SpatialPointsDataFrame representing the vertices
#'         of the graph
#'         \item digits : the number of digits keeped for the coordinates
#' }
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @examples
#' #This is an internal function, no example provided
build_graph_directed <- function(lines, digits, line_weight, direction, attrs = FALSE) {
  # doubling the lines if needed
  all_lines <- direct_lines(lines,direction)
  # extracting lines coordinates
  extremites <- lines_extremities(all_lines)
  start_coords <- extremites@data[extremites$pttype == "start", c("X","Y")]
  end_coords <- extremites@data[extremites$pttype == "end", c("X", "Y")]
  # extracting the coordinates of the starting and end points start_coords
  start <- sp_char_index(start_coords, digits)
  end <- sp_char_index(end_coords, digits)

  weights <- all_lines[[line_weight]]

  # building the line list
  linelist <- data.frame(start = start, end = end, weight = weights,
                         graph_id = 1:nrow(all_lines))
  if (attrs) {
    linelist <- cbind(linelist, all_lines@data)
  }
  graph <- igraph::graph_from_data_frame(linelist, directed = T, vertices = NULL)
  vertices <- igraph::V(graph)
  dfvertices <- data.frame(name = names(vertices), id = as.vector(vertices))
  dfvertices$name <- as.character(dfvertices$name)
  dfvertices <- tidyr::separate(dfvertices, col = "name", into = c("x",
                                                                   "y"), sep = "_", remove = F)
  dfvertices$x <- as.numeric(dfvertices$x)
  dfvertices$y <- as.numeric(dfvertices$y)
  points <- SpatialPoints(dfvertices[c("x", "y")])
  points <- SpatialPointsDataFrame(points, dfvertices)
  raster::crs(points) <- raster::crs(lines)
  return(list(graph = graph, linelist = linelist, lines = all_lines,
              spvertices = points, digits = digits))
}


#' function to match some points (SpatialPointsDataFrame) to the vertices of
#' a graph
#'
#' @param spvertices The spatial vertices of a graph (produced whith
#'build_graph)
#' @param points the SpatialPointsDataFrame to match
#' @param digits The number of digits to keep from the coordinates
#' @param tol the maximum distance between a point and a vertex
#' @return A vector of the corresponding vertices id, multiple points may
#' belong to the same vertex
#' @importFrom sp coordinates SpatialPoints SpatialPointsDataFrame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rgeos gIntersects gBuffer
#' @examples
#' #This is an internal function, no example provided
find_vertices <- function(spvertices, points, digits, tol = 0.1) {
    # step1 : calculate the spatialnameindex
    coords <- coordinates(points)
    points$tempoid <- 1:nrow(points)
    points$spIndex <- sp_char_index(coords, digits)
    # step2 : check which points are already well matched
    matching <- dplyr::left_join(points@data, spvertices@data, by = c(spIndex = "name"),
        keep = T)
    matching <- matching %>% dplyr::group_by(tempoid) %>%
      dplyr::summarise_all(dplyr::first)

    if (any(is.na(matching$id)) == F) {
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
                  prepared = T, byid = T))
                if (any(test) == F) {
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


#' function to create complementary lines for a directed network
#'
#' @param lines The original SpatialLinesDataFrame
#' @param direction A vector of integers. 0 indicates a bidirectional line and 1
#' an unidirectional line
#' @return A SpatialLinesDataFrame with some lines dupplicated according to
#' direction
#' @importFrom sp coordinates Line Lines SpatialLinesDataFrame SpatialLines
#' @importFrom raster crs
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
  df <- SpatialLinesDataFrame(splines,data,match.ID = F)
  raster::crs(df) <- raster::crs(lines)
  return(df)
}



#' function to plot a graph (usefull to check connectivity)
#'
#' @param graph A graph object (produced with build_graph)
#' @examples
#' #This is an internal function, no example provided
plot_graph <- function(graph) {
    N <- data.frame(name = names(igraph::V(graph)), id = as.vector(igraph::V(graph)))
    N <- tidyr::separate(N, "name", into = c("x", "y"), sep = "_", remove = F)
    N$x <- as.numeric(N$x)
    N$y <- as.numeric(N$y)

    graphics::plot(graph, vertex.size = 1, layout = as.matrix(N[c("x", "y")]), vertex.label.cex = 0.1)

}
