context("testing the isochrone creation functions")

test_that("Testing the isochrone function return the expected result on a simple example", {

  #definition of a simple situation
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

  # calculating the isochrone
  observed <- calc_isochrones(all_lines, 3.5, event[1,],
                              mindist = 1)

  # the expected lines in the obtained isochrones
  expected <- wkt_lines <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)",
    "LINESTRING (1 1, 2 1)",
    "LINESTRING (0 2, 1 2)",
    "LINESTRING (0 0, 0 1)",
    "LINESTRING (0 1, 0 2)",
    "LINESTRING (0 2, 0 3)",
    "LINESTRING (1 0, 1 1)",
    "LINESTRING (1 1, 1 2)",
    "LINESTRING (2 0, 2 1)",
    "LINESTRING (2 1, 2.5 1)",
    "LINESTRING (1 2, 1.5 2)",
    "LINESTRING (0 3, 0.5 3)",
    "LINESTRING (1 2, 1 2.5)",
    "LINESTRING (2 1, 2 1.5)",
    "LINESTRING (3 0, 3 0.5)"
  )


  observed$wkt <- st_as_text(st_geometry(observed))
  observed$wkt <- gsub("00","",observed$wkt, fixed = T)
  observed$wkt <- gsub(". "," ",observed$wkt, fixed = T)
  observed$wkt <- gsub(".)",")",observed$wkt, fixed = T)
  observed$wkt <- gsub(".,",",",observed$wkt, fixed = T)
  observed$wkt <- gsub(".50",".5",observed$wkt, fixed = T)

  test <- expected %in% observed$wkt

  expect_false(any(!test))
})


test_that("Testing the isochrone function return the expected result on a simple example with directions", {

  #definition of a simple situation
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
                        id = paste("l",1:length(wkt_lines),sep=""),
                        direction = "Both")
  linesdf[2,"direction"] <- "TF"

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of one event
  event <- data.frame(x=c(0, 1, 3),
                      y=c(0, 2, 3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))

  # calculating the isochrone
  observed <- calc_isochrones(all_lines, 2, event[1,],
                              mindist = 1, direction = "direction")

  # the expected lines in the obtained isochrones
  expected <- wkt_lines <- c(
    "LINESTRING (1 0, 0 0)",
    "LINESTRING (0 1, 0 0)",
    "LINESTRING (1 1, 1 0)",
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (0 2, 0 1)",
    "LINESTRING (0 0, 0 1)",
    "LINESTRING (1 1, 0 1)",
    "LINESTRING (1 0, 1 1)",
    "LINESTRING (0 1, 1 1)",
    "LINESTRING (0 1, 0 2)",
    "LINESTRING (1 1, 1 1)",
    "LINESTRING (1 1, 1 1)",
    "LINESTRING (0 2, 0 2)",
    "LINESTRING (0 2, 0 2)",
    "LINESTRING (1 1, 1 1)",
    "LINESTRING (1 1, 1 1)",
    "LINESTRING (0 2, 0 2)",
    "LINESTRING (0 2, 0 2)"
  )


  #observed$wkt <- st_as_text(st_geometry(observed))
  observed$wkt <- st_as_text(st_geometry(observed))
  observed$wkt <- gsub("00","",observed$wkt, fixed = T)
  observed$wkt <- gsub(". "," ",observed$wkt, fixed = T)
  observed$wkt <- gsub(".)",")",observed$wkt, fixed = T)
  observed$wkt <- gsub(".,",",",observed$wkt, fixed = T)
  observed$wkt <- gsub(".50",".5",observed$wkt, fixed = T)

  test1 <- expected %in% observed$wkt
  test2 <- observed$wkt %in% expected
  test <- c(test1, test2)

  expect_false(any(!test))
})

