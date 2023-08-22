context("testing the base kernel implemented")
library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TESTING THE split_by_grid.mc and split_by_grid functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing that split_by_grid.mc and split_by_grid return the same thing", {

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
  buff <- sf::st_buffer(all_lines, dist = 2)

  # definition of one event
  event <- data.frame(x=c(0, 1, 3),
                      y=c(0, 2, 3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))

  # creating the grid
  grid <- build_grid(c(3,3), spatial = list(buff, all_lines, event))


  ## testing with split all
  elements1 <- split_by_grid(grid, event, event, all_lines, bw = 1, tol = 0.1, digits = 2, split_all = T)
  future::plan(future::multisession(workers=1))
  elements2 <- split_by_grid.mc(grid = grid,
                                events = event,
                                samples = event,
                                lines = all_lines,
                                bw = 1,
                                tol = 0.1,
                                digits = 2,
                                split_all = T)
  ord1 <- elements1[order(sapply(elements1, function(i){i$samples$id}))]
  ord2 <- elements2[order(sapply(elements2, function(i){i$samples$id}))]
  tests <- sapply(1:length(ord1), function(i){
    t1 <- any(ord1[[i]]$samples$id != ord2[[i]]$samples$id)
    t2 <- any(ord1[[i]]$events$id != ord2[[i]]$events$id)
    t3 <- any(ord1[[i]]$lines$oid != ord2[[i]]$lines$oid)
    return(t1 & t2 & t3)
  })
  test1 <- (any(tests))

  ## testing without split all
  elements1 <- split_by_grid(grid, event, event, all_lines, bw = 1, tol = 0.1, digits = 2, split_all = FALSE)
  future::plan(future::multisession(workers=1))
  elements2 <- split_by_grid.mc(grid = grid,
                                events = event,
                                samples = event,
                                lines = all_lines,
                                bw = 1,
                                tol = 0.1,
                                digits = 2,
                                split_all = FALSE)
  ord1 <- elements1[order(sapply(elements1, function(i){i$samples$id}))]
  ord2 <- elements2[order(sapply(elements2, function(i){i$samples$id}))]
  tests <- sapply(1:length(ord1), function(i){
    t1 <- any(ord1[[i]]$samples$id != ord2[[i]]$samples$id)
    t2 <- any(ord1[[i]]$events$id != ord2[[i]]$events$id)
    t3 <- any(ord1[[i]]$lines$oid != ord2[[i]]$lines$oid)
    return(t1 & t2 & t3)
  })
  test2 <- (any(tests))

  expect_true(!test1 & !test2)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TESTING THE split_by_grid.mc and split_by_grid functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing that split_by_grid_abw.mc and split_by_grid_abw return the same thing", {

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
  buff <- sf::st_buffer(all_lines, dist = 2)

  # definition of one event
  event <- data.frame(x=c(0, 1, 3),
                      y=c(0, 2, 3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))

  # creating the grid
  grid <- build_grid(c(3,3), spatial = list(buff, all_lines, event))


  ## testing with split all
  elements1 <- split_by_grid_abw(grid, event, all_lines, bw = 1, tol = 0.1, digits = 2)
  future::plan(future::multisession(workers=1))
  elements2 <- split_by_grid_abw.mc(grid, event, all_lines, bw = 1, tol = 0.1, digits = 2)

  tests <- sapply(1:length(elements1), function(i){
    t2 <- any(elements1[[i]]$events$id != elements1[[i]]$events$id)
    t3 <- any(elements2[[i]]$lines$oid != elements2[[i]]$lines$oid)
    return(t2 & t3)
  })
  test1 <- (any(tests))

  expect_false(test1)
})
