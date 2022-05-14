context("testing the spatial indexing functions")
library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### First simple test in comparison with sf intersection ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the spatial indexing functions with a simple case", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)",
    "LINESTRING (1 1, 2 1)",
    "LINESTRING (2 1, 3 1)",
    "LINESTRING (0 2, 1 2)",
    "LINESTRING (1 2, 2 2)",
    "LINESTRING (2 2, 3 2)",
    "LINESTRING (0 3, 1 3)",
    "LINESTRING (1 3, 2 3)",
    "LINESTRING (2 3, 3 3)",
    "LINESTRING (0 0, 0 1)",
    "LINESTRING (0 1, 0 2)",
    "LINESTRING (0 2, 0 3)",
    "LINESTRING (1 0, 1 1)",
    "LINESTRING (1 1, 1 2)",
    "LINESTRING (1 2, 1 3)",
    "LINESTRING (2 0, 2 1)",
    "LINESTRING (2 1, 2 2)",
    "LINESTRING (2 2, 2 3)",
    "LINESTRING (3 0, 3 1)",
    "LINESTRING (3 1, 3 2)",
    "LINESTRING (3 2, 3 3)"
  )

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")
  poly <- st_buffer(all_lines[5,],0.5)
  tree <- build_quadtree(all_lines)
  obtained <- spatial_request(poly, tree, all_lines)

  inter <- lengths(st_intersects(all_lines, poly)) > 0
  expected <- all_lines[inter,]

  test <- any(!(expected$id[order(expected$id)] == obtained$id[order(obtained$id)]))
  expect_false(test)
})

