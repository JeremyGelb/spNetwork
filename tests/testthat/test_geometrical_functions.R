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

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # creating the expected spatialPointsDataFrame
  df <- data.frame(
    "wkt" = wkt_lines[c(1,1,2,2,3,3,4,4)],
    "id" = paste("l",c(1,1,2,2,3,3,4,4), sep = ""),
    "X" = c(0,1,1,2,2,3,0,1),
    "Y" = c(0,0,0,0,0,0,1,1),
    "pttype" = c("start","end","start","end","start","end","start","end")
  )


  # creating the obtained SpatialPointsDataFrame
  obtained <- lines_extremities(all_lines)

  test <- any(obtained@data != df)

  expect_false(test)
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

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # applying the remove loop function
  new_lines <- remove_loop_lines(all_lines,2)

  # only the line l3 must be missing
  diff <- all_lines$id[!(all_lines$id %in% new_lines$id)]
  test <- diff == c("l3")

  expect_true(test)
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

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = FALSE)

  # applying the remove loop function
  new_lines <- reverse_lines(all_lines)

  # getting the new coordinates
  new_coords <- unlist(sp::coordinates(new_lines), recursive = FALSE)

  # creating the expected result
  expected <- list(
    cbind(c(1,0),c(0,0)),
    cbind(c(2,1),c(0,0)),
    cbind(c(1,3,3,2,1),c(0,3,0,0,0)),
    cbind(c(3,2), c(0,0)),
    cbind(c(1,0),c(1,1))
  )

  vecbool <- sapply(1:length(expected), function(i){
    return(any(expected[[i]] != new_coords[[i]]))
  })

  expect_false(any(vecbool))
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

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = FALSE)

  # creating two points
  event <- data.frame(x=c(0.5, 2.7),
                      y=c(2, 0.4))
  sp::coordinates(event) <- cbind(event$x,event$y)

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
  new_coords <- unlist(sp::coordinates(new_lines), recursive = FALSE)

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

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = FALSE)


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
  new_coords <- unlist(sp::coordinates(new_lines), recursive = FALSE)

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

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = FALSE)

  observed <- add_center_lines(all_lines)
  observed$wkt <- rgeos::writeWKT(observed, byid = T)
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

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # calculating the observed result
  lixels <- lixelize_lines(all_lines, 1, 0.5)
  observed <- rgeos::writeWKT(lixels, byid = TRUE)
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

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # calculating the observed result
  future::plan(future::multisession(workers=1))
  lixels <- lixelize_lines.mc(all_lines, 1, 0.5)
  observed <- rgeos::writeWKT(lixels, byid = TRUE)
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
