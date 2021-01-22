context("testing the kernel functions")

test_that("Testing the simple kernel with a simple case", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 5.1, 0.1 0.1)",
    "LINESTRING (-5.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -5.1, 0.1 0.1)",
    "LINESTRING (5.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1))
  sp::coordinates(event) <- cbind(event$x,event$y)

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                      y=c(0.1))
  sp::coordinates(sp_point) <- cbind(sp_point$x,sp_point$y)

  #sp::proj4string(event) <- sp::CRS("EPSG:32633")
  #sp::proj4string(sp_point) <- sp::CRS("EPSG:32633")
  #sp::proj4string(all_lines) <- sp::CRS("EPSG:32633")

  # real distance is 6, and let us say that the bw is 10
  real_value <- (1/10) * quartic_kernel(6,10)

  #let us calculate the value with our function
  obs_value <- nkde(all_lines,events = event, w = c(1),
       samples = sp_point, check = F,
       kernel_name = "quartic",
       bw = 10, adaptive = F, method = "simple", div = "bw",
       agg = 0.01, verbose = F,tol = 0.0001
       )
  expect_equal(obs_value, real_value)
})


test_that("Testing the discontinuous kernel with a simple case", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 5.1, 0.1 0.1)",
    "LINESTRING (-5.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -5.1, 0.1 0.1)",
    "LINESTRING (5.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1))
  sp::coordinates(event) <- cbind(event$x,event$y)

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                         y=c(0.1))
  sp::coordinates(sp_point) <- cbind(sp_point$x,sp_point$y)

  #sp::proj4string(event) <- sp::CRS("+init=EPSG:32633")
  #sp::proj4string(sp_point) <- sp::CRS("+init=EPSG:32633")
  #sp::proj4string(all_lines) <- sp::CRS("+init=EPSG:32633")

  # real distance is 6, and let us say that the bw is 10
  real_value <- (1/10) * quartic_kernel(6,10) * 1/3

  #let us calculate the value with our function
  obs_value <- nkde(all_lines,events = event, w = c(1),
                    samples = sp_point, check = F,
                    kernel_name = "quartic",
                    bw = 10, adaptive = F, method = "discontinuous", div = "bw",
                    agg = 0.01, verbose = F,tol = 0.0001
  )
  expect_equal(obs_value, real_value)
})


test_that("Testing the continuous kernel with a simple case", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 10.1, 0.1 0.1)",
    "LINESTRING (-10.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -10.1, 0.1 0.1)",
    "LINESTRING (10.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1))
  sp::coordinates(event) <- cbind(event$x,event$y)

  # definition of one sampling point
  sp_point <- data.frame(x=c(0.1),
                         y=c(1.1))
  sp::coordinates(sp_point) <- cbind(sp_point$x,sp_point$y)

  #sp::proj4string(event) <- sp::CRS("+init=EPSG:32633")
  #sp::proj4string(sp_point) <- sp::CRS("+init=EPSG:32633")
  #sp::proj4string(all_lines) <- sp::CRS("+init=EPSG:32633")

  bw <- 5
  # real distance is 2, and let us say that the bw is 5
  real_value <- quartic_kernel(2,bw) - ((1/2) * quartic_kernel(4,bw))

  #let us calculate the value with our function
  obs_value <- nkde(all_lines,events = event, w = c(1),
                    samples = sp_point,check = F,
                    kernel_name = "quartic",
                    bw = bw, adaptive = F, method = "continuous", div = "none",
                    agg = 0.01, verbose = F,tol = 0.01,digits = 3,sparse = T,
  )
  expect_equal(obs_value, real_value)
})

