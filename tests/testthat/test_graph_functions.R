context("testing the graph functions")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TESTING THE FUNCTION TO CHECK IF A GRAPH IS WELL CONNECTED ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the function which check if a graph is valid", {

  # creating a simple case
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)
  result <- graph_checking(all_lines, digits = 2, max_search = 1)

  test1 <- nrow(result$dangle_nodes)== 4
  test2 <- unique(result$vertex_components$component) == 1
  test3 <- nrow(result$close_nodes)== 0

  expect_true(test1 & test2 & test3)

})



test_that("Testing the function which check if a graph is valid, with an invalid graph", {

  # creating a simple case
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 3 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)
  result <- graph_checking(all_lines, digits = 2)

  test1 <- nrow(result$dangle_nodes)== 5
  test2 <- any(!unique(result$vertex_components$component) == c(1,2)) == FALSE

  expect_true(test1 & test2)

})



