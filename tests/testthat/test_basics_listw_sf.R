context("basic tests for the neighbour matrices creation")
library(sf)

test_that("Testing on a simple scenario if the neighbouring functions are producing expected results", {

  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)",
    "LINESTRING (5 0, 5 -5)"
    )

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf,
                         wkt = "wkt")

  all_lines$OID <- 1:nrow(all_lines)

  # definition of three event
  event <- data.frame(x=c(0, 3, -4, 4.8),
                      y=c(3, 0, 0, -4.8))
  event$OID <- 1:nrow(event)

  event <- st_as_sf(event, coords = c("x","y"))

  listw <- network_listw(event,
                         all_lines,
                         maxdistance = 8,
                         dist_func = "squared inverse",
                         matrice_type = "B",
                         grid_shape = c(1,1),
                         verbose = TRUE,
                         mindist = 0.1,
                         digits = 3
  )
  obtainedmat <- spdep::nb2mat(listw$neighbours,style = "B")

  expected <- rbind(
    c(0,1,1,0),
    c(1,0,1,1),
    c(1,1,0,0),
    c(0,1,0,0)
  )

  tottest <- sum(abs(obtainedmat - expected)) == 0
  expect_true(tottest)
})


test_that("Testing on a simple scenario if the neighbouring functions are producing expected results (with polygons)", {

  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)",
    "LINESTRING (5 0, 5 -5)"
  )

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf,
                        wkt = "wkt")
  all_lines$OID <- 1:nrow(all_lines)

  # definition of three polygons
  event <- data.frame(x=c(0, 3, -4, 4.8),
                      y=c(3, 0, 0, -4.8))
  event$OID <- 1:nrow(event)
  event <- st_as_sf(event, coords = c("x","y"))


  eventbuff <- sf::st_buffer(event, dist = 1)


  listw1 <- network_listw(eventbuff,
                         all_lines,
                         method = "pointsalong",
                         point_dist = 0.25,
                         maxdistance = 8,
                         dist_func = "squared inverse",
                         matrice_type = "B",
                         grid_shape = c(1,1),
                         verbose = FALSE,
                         mindist = 0.1,
                         digits = 3
  )

  listw2 <- network_listw(eventbuff,
                          all_lines,
                          method = "pointsalong",
                          point_dist = 0.25,
                          maxdistance = 12,
                          dist_func = "squared inverse",
                          matrice_type = "B",
                          grid_shape = c(1,1),
                          verbose = FALSE,
                          mindist = 0.1,
                          digits = 3
  )

  obtainedmat1 <- spdep::nb2mat(listw1$neighbours,style = "B")
  obtainedmat2 <- spdep::nb2mat(listw2$neighbours,style = "B")

  expected1 <- rbind(
    c(0,1,1,0),
    c(1,0,1,1),
    c(1,1,0,0),
    c(0,1,0,0)
  )

  expected2 <- rbind(
    c(0,1,1,1),
    c(1,0,1,1),
    c(1,1,0,1),
    c(0,1,1,0)
  )

  tottest <- sum(abs(obtainedmat1 - expected1)) == 0 & sum(abs(obtainedmat2 - expected2))
  expect_true(tottest)
})
