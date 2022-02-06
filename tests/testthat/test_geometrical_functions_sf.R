context("Testing the geometrical functions of the package")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the sp_char_index function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_that("formal testing of the sp_char_index function", {

  # creating a spatialPointsDataFrame
  event <- data.frame(x=c(0.1, 1.5, 3.3),
                      y=c(0.1, 2.5, 3.3))
  coords <- cbind(event$x,event$y)

  text <- sp_char_index(coords, digits = 1)

  expected <- c("0.1_0.1", "1.5_2.5", "3.3_3.3")

  expect_identical(text,expected)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the lines_extremities function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("formal testing of the lines_extremities function", {

  # creating a spatialLinesDataFrame
  wkt_lines <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")
  all_lines <- cbind(linesdf$wkt,all_lines)

  # creating the expected spatialPointsDataFrame
  df <- data.frame(
    "wkt" = wkt_lines[c(1,1,2,2,3,3,4,4)],
    "id" = paste("l",c(1,1,2,2,3,3,4,4), sep = ""),
    "pttype" = c("start","end","start","end","start","end","start","end"),
    "X" = c(0,1,1,2,2,3,0,1),
    "Y" = c(0,0,0,0,0,0,1,1)
  )


  # creating the obtained SpatialPointsDataFrame
  obtained <- lines_extremities(all_lines)

  test <- any(st_drop_geometry(obtained) != df)

  expect_false(test)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the lines_center  function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("formal testing of the lines_center  function", {

  # creating a spatialLinesDataFrame
  wkt_lines <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)",
    "LINESTRING (0 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")
  all_lines <- cbind(linesdf$wkt,all_lines)

  # creating the expected spatialPointsDataFrame
  expected <- rbind(
    c(0.5,0),
    c(1.5,0),
    c(2.5,0),
    c(0.5,1),
    c(0,0)
  )
  obtained <- lines_center(all_lines)

  # creating the obtained SpatialPointsDataFrame
  obtained_mat <- st_coordinates(obtained)

  test <- sum(obtained_mat - expected)

  expect_true(test == 0)
})



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the remove_loop_lines function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("formal testing of the remove_loop_lines function", {

  # creating a spatialLinesDataFrame
  wkt_lines <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (1 0, 2 0, 3 0, 3 3, 1 0)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # applying the remove loop function
  new_lines <- remove_loop_lines(all_lines,2)

  # only the line l3 must be missing
  diff <- all_lines$id[!(all_lines$id %in% new_lines$id)]
  test <- diff == c("l3")

  expect_true(test)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the remove_mirror_edges function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("formal testing of the remove_mirror_edges function", {

  # creating a spatialLinesDataFrame
  wkt_lines <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (0 0, 0 5, 1 0)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # applying the remove loop function
  new_lines1 <- remove_mirror_edges(all_lines, keep_shortest = TRUE)
  new_lines2 <- remove_mirror_edges(all_lines, keep_shortest = FALSE, verbose = FALSE)

  # only the line l3 must be missing
  test1 <- "l2" %in% new_lines1$id == FALSE
  test2 <- "l1" %in% new_lines2$id == FALSE

  expect_true(test1 & test2)
})



# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# #### Testing the reverse_lines function ####
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("formal testing of the reverse_lines function", {

  # creating a spatialLinesDataFrame
  wkt_lines <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (1 0, 2 0, 3 0, 3 3, 1 0)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # applying the remove loop function
  new_lines <- reverse_lines(all_lines)

  # getting the new coordinates
  new_coords <- st_coordinates(new_lines)

  # creating the expected result
  expected <- rbind(
    cbind(c(1,0),c(0,0)),
    cbind(c(2,1),c(0,0)),
    cbind(c(1,3,3,2,1),c(0,3,0,0,0)),
    cbind(c(3,2), c(0,0)),
    cbind(c(1,0),c(1,1))
  )

  test<- sum(abs(new_coords[,1:2] - expected)) == 0

  expect_true(test)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Testing the add_vertices_lines function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_that("formal testing of the add_vertices_lines function", {

  # creating a spatialLinesDataFrame
  wkt_lines <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (3 0, 3 3, 1 0)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # creating two points
  event <- data.frame(x=c(0.5, 2.7),
                      y=c(2, 0.4))
  event <- st_as_sf(event, coords = c("x","y"))

  # writing the expected new line coordinates
  expected <- list(
    rbind(c(0,0),c(1,0)),
    rbind(c(1,0),c(2,0)),
    rbind(c(3,0),c(3,0.4),c(3,3),c(1,0)),
    rbind(c(2,0), c(3,0)),
    rbind(c(0,1),c(0.5,1),c(1,1))
  )

  # getting the new_lines coordinates
  new_lines <- add_vertices_lines(all_lines, event, nearest_lines_idx = c(5,3), mindist = 0.1)
  new_coords <- lines_coordinates_as_list(new_lines)

  values <- sapply(1:length(expected), function(i){
    return(sum(abs(expected[[i]] - new_coords[[i]])))
  })

  expect_equal(sum(values),0)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the simple_lines function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("formal testing of the simple_lines function", {

  # creating a spatialLinesDataFrame
  wkt_lines <- c(
    "LINESTRING (0 0, 1 0, 1 1)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (3 0, 3 3, 0 3)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")


  # writing the expected new line coordinates
  expected <- list(
    rbind(c(0,0),c(1,0)),
    rbind(c(1,0),c(1,1)),
    rbind(c(1,0),c(2,0)),
    rbind(c(3,0),c(3,3)),
    rbind(c(3,3),c(0,3)),
    rbind(c(2,0), c(3,0)),
    rbind(c(0,1), c(1,1))
  )

  # getting the new_lines coordinates
  new_lines <- simple_lines(all_lines)
  new_coords <- lines_coordinates_as_list(new_lines)

  values <- sapply(1:length(expected), function(i){
    return(sum(abs(expected[[i]] - new_coords[[i]])))
  })

  expect_equal(sum(values),0)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the add_center_lines function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("formal testing of the simple_lines function", {

  # creating a spatialLinesDataFrame
  wkt_lines <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (1 0, 2 0)",
    "LINESTRING (3 0, 3 3)",
    "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  observed <- add_center_lines(all_lines)

  observed$wkt <- st_as_text(observed$geometry)
  observed$wkt <- gsub("00","",observed$wkt, fixed = T)
  observed$wkt <- gsub(". "," ",observed$wkt, fixed = T)
  observed$wkt <- gsub(".)",")",observed$wkt, fixed = T)
  observed$wkt <- gsub(".,",",",observed$wkt, fixed = T)
  observed$wkt <- gsub(".50",".5",observed$wkt, fixed = T)

  expected <- c("LINESTRING (0 0, 0.5 0, 1 0)",
                "LINESTRING (1 0, 1.5 0, 2 0)",
                "LINESTRING (3 0, 3 1.5, 3 3)",
                "LINESTRING (2 0, 2.5 0, 3 0)",
                "LINESTRING (0 1, 0.5 1, 1 1)")

  test <- c(observed$wkt %in% expected, expected %in% observed$wkt)

  expect_false(any(!test))

})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the lixelize_lines function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


test_that("formal testing of the lixelize_lines function", {

  # creating some line
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # calculating the observed result
  lixels <- lixelize_lines(all_lines, 1, 0.5)

  observed <-st_as_text(lixels$geometry)
  observed <- gsub("00","",observed, fixed = TRUE)
  observed <- gsub(".","",observed, fixed = TRUE)

  # noting the expected results
  expected <- c("LINESTRING (0 5, 0 4)",
  "LINESTRING (0 4, 0 3)",
  "LINESTRING (0 3, 0 2)",
  "LINESTRING (0 2, 0 1)",
  "LINESTRING (0 1, 0 0)",
  "LINESTRING (-5 0, -4 0)",
  "LINESTRING (-4 0, -3 0)",
  "LINESTRING (-3 0, -2 0)",
  "LINESTRING (-2 0, -1 0)",
  "LINESTRING (-1 0, 0 0)" ,
  "LINESTRING (0 -5, 0 -4)",
  "LINESTRING (0 -4, 0 -3)",
  "LINESTRING (0 -3, 0 -2)",
  "LINESTRING (0 -2, 0 -1)",
  "LINESTRING (0 -1, 0 0)",
  "LINESTRING (5 0, 4 0)",
  "LINESTRING (4 0, 3 0)",
  "LINESTRING (3 0, 2 0)",
  "LINESTRING (2 0, 1 0)",
  "LINESTRING (1 0, 0 0)")

  test <- c(expected %in% observed, observed %in% expected)
  expect_false(any(!test))

})


test_that("formal testing of the lixelize_lines.mc function", {

  # creating some line
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # calculating the observed result
  future::plan(future::multisession(workers=1))

  lixels <- lixelize_lines.mc(all_lines, 1, 0.5)
  lixels2 <- lixelize_lines.mc(all_lines, 1, 0.5, verbose = FALSE)

  observed <- st_as_text(lixels$geometry, byid = TRUE)
  observed <- gsub("00","",observed, fixed = TRUE)
  observed <- gsub(".","",observed, fixed = TRUE)

  observed2 <- st_as_text(lixels2$geometry, byid = TRUE)
  observed2 <- gsub("00","",observed2, fixed = TRUE)
  observed2 <- gsub(".","",observed2, fixed = TRUE)

  # noting the expected results
  expected <- c("LINESTRING (0 5, 0 4)",
                "LINESTRING (0 4, 0 3)",
                "LINESTRING (0 3, 0 2)",
                "LINESTRING (0 2, 0 1)",
                "LINESTRING (0 1, 0 0)",
                "LINESTRING (-5 0, -4 0)",
                "LINESTRING (-4 0, -3 0)",
                "LINESTRING (-3 0, -2 0)",
                "LINESTRING (-2 0, -1 0)",
                "LINESTRING (-1 0, 0 0)" ,
                "LINESTRING (0 -5, 0 -4)",
                "LINESTRING (0 -4, 0 -3)",
                "LINESTRING (0 -3, 0 -2)",
                "LINESTRING (0 -2, 0 -1)",
                "LINESTRING (0 -1, 0 0)",
                "LINESTRING (5 0, 4 0)",
                "LINESTRING (4 0, 3 0)",
                "LINESTRING (3 0, 2 0)",
                "LINESTRING (2 0, 1 0)",
                "LINESTRING (1 0, 0 0)")

  test1 <- c(expected %in% observed, observed %in% expected)
  test2 <- c(expected %in% observed, observed %in% expected)
  expect_false(any(!test1)|any(!test2))

})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the healing_edges function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("formal testing of the heal_edges function", {

  # creating some line
  wkt_lines <- c(
    "LINESTRING (0 10, 0 5)",
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (5 0, 0 0)",
    "LINESTRING (5 0, 5 5)",
    "LINESTRING (0 0, 0 -5)"
    )

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  new_edges <- heal_edges(all_lines)

  obtained_lengths <- as.numeric(st_length(new_edges))
  expected_lengths <- c(5,10,10)
  res <- sum(abs(obtained_lengths - expected_lengths))

  n_obs <- nrow(new_edges)

  test1 <- res==0
  test2 <- n_obs == 3
  expect_true(test1 & test2)

})
