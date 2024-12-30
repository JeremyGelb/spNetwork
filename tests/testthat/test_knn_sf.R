context("testing the knn functions")
library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### First simple test ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the knn function with a simple case", {
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

  # definition of one event
  event <- data.frame(x=c(0, 1, 3),
                      y=c(0, 2, 3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))


  ## expected matrices
  exp_dist <- rbind(
    c(3,6),
    c(3,3),
    c(3,6)
  )
  exp_oid <- rbind(
    c(2,3),
    c(1,3),
    c(2,1)
  )

  ## calculating the realvalues
  mats <- network_knn(origins = event,
                      lines = all_lines,
                      k = 2,
                      maxdistance = 100,
                      tol = 0.001)

  test1 <- sum(mats[[1]] != exp_dist)
  test2 <- sum(mats[[2]] != exp_oid)
  expect_true(test1 == 0 & test2 == 0)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Second simple test ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the knn function with a simple case and a directed network", {
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
    # "LINESTRING (0 0, 0 1)",
    # "LINESTRING (0 1, 0 2)",
    # "LINESTRING (0 2, 0 3)",
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
  linesdf$dir <- "Both"
  linesdf[c(8,5),"dir"] <- "FT"
  linesdf[c(13),"dir"] <- "TF"

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of one event
  event <- data.frame(x=c(0, 1, 3),
                      y=c(0, 2, 3))
  event <- st_as_sf(event, coords = c("x","y"))

  # sp::plot(all_lines)
  # sp::plot(event, add = T, pch = 21)
  # sp::plot(subset(all_lines, all_lines$dir != "Both"), c = "red", add = T)

  ## expected matrices
  exp_dist <- rbind(
    c(6,7),
    c(3,3),
    c(3,6)
  )
  exp_oid <- rbind(
    c(3,2),
    c(1,3),
    c(2,1)
  )

  ## calculating the realvalues
  mats <- network_knn(origins = event,
                      lines = all_lines,
                      k = 2,
                      maxdistance = 100,
                      direction = "dir")
  test1 <- sum(mats[[1]] != exp_dist)
  test2 <- sum(mats[[2]] != exp_oid)
  expect_true(test1 == 0 & test2 == 0)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Third simple test ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the knn function with a simple case and specific destinations", {
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

  # definition of some events
  event <- data.frame(x=c(0, 1, 3),
                      y=c(0, 2, 3))
  event <- st_as_sf(event, coords = c("x","y"))

  # definition of some destinations
  destinations <- data.frame(x=c(0, 3),
                      y=c(3, 0))
  destinations <- st_as_sf(destinations, coords = c("x","y"))

  # sp::plot(all_lines)
  # sp::plot(event, add = T, pch = 21)
  # sp::plot(destinations, add = T, pch = 20)

  ## expected matrices
  exp_dist <- rbind(
    c(3,3),
    c(2,4),
    c(3,3)
  )
  exp_oid <- rbind(
    c(1,2),
    c(1,2),
    c(1,2)
  )

  ## calculating the realvalues
  mats <- network_knn(event, all_lines, destinations = destinations, k = 2, maxdistance = 100)
  test1 <- sum(mats[[1]] - exp_dist)
  test2 <- sum(mats[[2]] - exp_oid)

  test3 <- (test1 + test2) == 0
  expect_true(test3)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Fourth simple test ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the knn function with a simple case, specific destinations and directions", {
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
    # "LINESTRING (0 0, 0 1)",
    # "LINESTRING (0 1, 0 2)",
    # "LINESTRING (0 2, 0 3)",
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
  linesdf$dir <- "Both"
  linesdf[c(8,5),"dir"] <- "FT"
  linesdf[c(13),"dir"] <- "TF"

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of some events
  event <- data.frame(x=c(0, 1, 3),
                      y=c(0, 2, 3))
  event <- st_as_sf(event, coords = c("x","y"))

  # definition of some destinations
  destinations <- data.frame(x=c(0, 3),
                             y=c(3, 0))
  destinations <- st_as_sf(destinations, coords = c("x","y"))

  # sp::plot(all_lines)
  # sp::plot(subset(all_lines, all_lines$dir != "Both"), c = "red", add = T)
  # sp::plot(event, add = T, pch = 21)
  # sp::plot(destinations, add = T, pch = 20)

  ## expected matrices
  exp_dist <- rbind(
    c(3,7),
    c(2,4),
    c(3,3)
  )
  exp_oid <- rbind(
    c(2,1),
    c(1,2),
    c(1,2)
  )

  ## calculating the realvalues
  mats <- network_knn(event, all_lines, destinations = destinations, k = 2, maxdistance = 100, direction = "dir")
  test1 <- sum(mats[[1]] - exp_dist)
  test2 <- sum(mats[[2]] - exp_oid)

  test3 <- (test1 + test2) == 0
  expect_true(test3)

})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Fourth simple test BUT MULTICORE ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the knn function with a simple case, specific destinations and directions", {
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
    # "LINESTRING (0 0, 0 1)",
    # "LINESTRING (0 1, 0 2)",
    # "LINESTRING (0 2, 0 3)",
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
  linesdf$dir <- "Both"
  linesdf[c(8,5),"dir"] <- "FT"
  linesdf[c(13),"dir"] <- "TF"

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of some events
  event <- data.frame(x=c(0, 1, 3),
                      y=c(0, 2, 3))
  event <- st_as_sf(event, coords = c("x","y"))

  # definition of some destinations
  destinations <- data.frame(x=c(0, 3),
                             y=c(3, 0))
  destinations <- st_as_sf(destinations, coords = c("x","y"))

  # sp::plot(all_lines)
  # sp::plot(subset(all_lines, all_lines$dir != "Both"), c = "red", add = T)
  # sp::plot(event, add = T, pch = 21)
  # sp::plot(destinations, add = T, pch = 20)

  ## expected matrices
  exp_dist <- rbind(
    c(3,7),
    c(2,4),
    c(3,3)
  )
  exp_oid <- rbind(
    c(2,1),
    c(1,2),
    c(1,2)
  )

  ## calculating the realvalues
  future::plan(future::multisession(workers=1))
  mats <- network_knn.mc(origins = event,
                         lines = all_lines,
                         destinations = destinations,
                         snap_dist = Inf,
                         k = 2,
                         maxdistance = 100,
                         direction = "dir", line_weight = 'length',
                         grid_shape = c(1,1),
                         verbose = TRUE,
                         digits = 2,
                         tol = 0.01)
  test1 <- sum(mats[[1]] - exp_dist)
  test2 <- sum(mats[[2]] - exp_oid)

  test3 <- (test1 + test2) == 0
  expect_true(test3)

})



