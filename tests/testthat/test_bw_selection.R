context("testing the bandwidth selection functions")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR BW SELECTION WITH CV LIKELIHOOD ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the bw selection function with CV likelihood and simple kernel", {
  ## creating the simple situation
  # start with de definition of some lines
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

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3))
  sp::coordinates(event) <- cbind(event$x,event$y)


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  s1 <- (quartic_kernel(6,10) + quartic_kernel(6,10)) *(1/10)

  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- 3*log(s1)


  #let us calculate the value with our function
  obs_value <- bw_cv_likelihood_calc(c(10,10),1,
                                     lines = all_lines,
                                     events = event, w = c(1,1,1),
                                     check = F,
                                     kernel_name = "quartic",
                                     method = "simple",
                                     digits = 1,
                                     agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[1,2], total)
})



test_that("Testing the bw selection function with CV likelihood and discontinuous kernel", {
  ## creating the simple situation
  # start with de definition of some lines
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

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3))
  sp::coordinates(event) <- cbind(event$x,event$y)


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  # for each point, we split the kernel at the intersection (1/3)
  s1 <- (quartic_kernel(6,10) * 1/3 + quartic_kernel(6,10) * 1/3) *(1/10)

  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- 3*log(s1)


  #let us calculate the value with our function
  obs_value <- bw_cv_likelihood_calc(c(10,10),1,
                    lines = all_lines,
                    events = event, w = c(1,1,1),
                    check = F,
                    kernel_name = "quartic",
                    method = "discontinuous",
                    digits = 1,
                    agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[1,2], total)
})


test_that("Testing the bw selection function with CV likelihood and continuous kernel", {
  ## creating the simple situation
  # start with de definition of some lines
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

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3))
  sp::coordinates(event) <- cbind(event$x,event$y)


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  # for each point, we split the kernel at the intersection and their is some backfire
  # lets calculate the impact of en event on another
  n <- 4
  imp <- quartic_kernel(6,10) * (2.0/n)

  s1 <- (imp + imp) * (1/10)

  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- 3*log(s1)


  #let us calculate the value with our function
  obs_value <- bw_cv_likelihood_calc(c(10,10),1,
                                     lines = all_lines,
                                     events = event, w = c(1,1,1),
                                     check = F,
                                     kernel_name = "quartic",
                                     method = "continuous",
                                     digits = 1,
                                     agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[1,2], total)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR BW SELECTION WITH Van Lieshout's Criterion ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the bw selection function with Van Lieshout's Criterion and simple kernel", {
  ## creating the simple situation
  # start with de definition of some lines
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

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3))
  sp::coordinates(event) <- cbind(event$x,event$y)


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  # so the density at one event is
  s1 <- (quartic_kernel(6,10) + quartic_kernel(6,10) + quartic_kernel(0,10)) *(1/10)

  #so the score value is
  Wl <- rgeos::gLength(all_lines) /1000
  score <- (Wl - ((1/s1)*3))**2


  #let us calculate the value with our function
  obs_value <- bw_cvl_calc(c(10,10),1,
                                     lines = all_lines,
                                     events = event, w = c(1,1,1),
                                     check = F,
                                     kernel_name = "quartic",
                                     method = "simple",
                                     digits = 1,
                                     agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[1,2], score)
})


test_that("Testing the bw selection function with Van Lieshout's Criterion and discontinuous kernel", {
  ## creating the simple situation
  # start with de definition of some lines
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

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3))
  sp::coordinates(event) <- cbind(event$x,event$y)


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  # for each point, we split the kernel at the intersection (1/3)
  # so the density at one event is
  s1 <- (quartic_kernel(6,10)*1/3 + quartic_kernel(6,10)*1/3 + quartic_kernel(0,10)) *(1/10)

  #so the score value is
  Wl <- rgeos::gLength(all_lines) /1000
  score <- (Wl - ((1/s1)*3))**2


  #let us calculate the value with our function
  obs_value <- bw_cvl_calc(c(10,10),1,
                           lines = all_lines,
                           events = event, w = c(1,1,1),
                           check = F,
                           kernel_name = "quartic",
                           method = "discontinuous",
                           digits = 1,
                           agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[1,2], score)
})


test_that("Testing the bw selection function with Van Lieshout's Criterion and continuous kernel", {
  ## creating the simple situation
  # start with de definition of some lines
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

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3))
  sp::coordinates(event) <- cbind(event$x,event$y)


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  # for each point, we split the kernel at the intersection and their is some backfire
  # lets calculate the impact of an event on another
  n <- 4
  imp <- quartic_kernel(6,10) * (2.0/n)

  #and now the impact on itself (backfire)
  itself <- quartic_kernel(0,10) + (((2.0-n)/n)*quartic_kernel(6,10))

  s1 <- (imp + imp + itself) * (1/10)
  #so the score value is
  Wl <- rgeos::gLength(all_lines) /1000
  score <- (Wl - ((1/s1)*3))**2


  #let us calculate the value with our function
  obs_value <- bw_cvl_calc(c(10,10),1,
                           lines = all_lines,
                           events = event, w = c(1,1,1),
                           check = F,
                           kernel_name = "quartic",
                           method = "continuous",
                           digits = 1,
                           agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[1,2], score)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Comparing multicore and single core ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing that bw selection by cv-likelihood gives the same score in single and multicore", {

  networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
  eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
  mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
  bike_accidents <- rgdal::readOGR(eventsgpkg,layer="bike_accidents", verbose=FALSE)

  case1 <- bike_accidents[25,]
  buff <- rgeos::gBuffer(case1, width = 1000)

  inter_idx1 <-  gIntersects(buff, bike_accidents, byid=TRUE)
  inter_idx2 <-  gIntersects(buff, mtl_network, byid=TRUE)

  # subset
  events <- bike_accidents[as.vector(inter_idx1),]
  lines <- mtl_network[as.vector(inter_idx2),]

  ## multicore cv score
  future::plan(future::multisession(workers=1))
  cv_scores.mc <- bw_cv_likelihood_calc.mc(c(200,400),100,
                                 lines, events,
                                 rep(1,nrow(events)),
                                 "quartic", "discontinuous",
                                 diggle_correction = FALSE, study_area = NULL,
                                 max_depth = 8,
                                 digits=2, tol=0.1, agg=5,
                                 sparse=TRUE, grid_shape=c(2,2),
                                 sub_sample = 1, verbose=FALSE, check=TRUE)

  ## single core cv score
  cv_scores <- bw_cv_likelihood_calc(c(200,400),100,
                                        lines, events,
                                        rep(1,nrow(events)),
                                        "quartic", "discontinuous",
                                        diggle_correction = FALSE, study_area = NULL,
                                        max_depth = 8,
                                        digits=2, tol=0.1, agg=5,
                                        sparse=TRUE, grid_shape=c(2,2),
                                        sub_sample = 1, verbose=FALSE, check=TRUE)
  test <- cv_scores == cv_scores.mc
  expect_false(any(!test))

})


test_that("Testing that bw selection with Van Lieshout's Criterion gives the same score in single and multicore", {

  networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
  eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
  mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
  bike_accidents <- rgdal::readOGR(eventsgpkg,layer="bike_accidents", verbose=FALSE)

  case1 <- bike_accidents[25,]
  buff <- rgeos::gBuffer(case1, width = 1000)

  inter_idx1 <-  gIntersects(buff, bike_accidents, byid=TRUE)
  inter_idx2 <-  gIntersects(buff, mtl_network, byid=TRUE)

  # subset
  events <- bike_accidents[as.vector(inter_idx1),]
  lines <- mtl_network[as.vector(inter_idx2),]

  ## multicore cv score
  future::plan(future::multisession(workers=1))
  cv_scores.mc <- bw_cvl_calc.mc(c(200,400),100,
                                           lines, events,
                                           rep(1,nrow(events)),
                                           "quartic", "discontinuous",
                                           diggle_correction = FALSE, study_area = NULL,
                                           max_depth = 8,
                                           digits=2, tol=0.1, agg=5,
                                           sparse=TRUE, grid_shape=c(2,2),
                                           sub_sample = 1, verbose=FALSE, check=TRUE)

  ## single core cv score
  cv_scores <- bw_cvl_calc(c(200,400),100,
                                     lines, events,
                                     rep(1,nrow(events)),
                                     "quartic", "discontinuous",
                                     diggle_correction = FALSE, study_area = NULL,
                                     max_depth = 8,
                                     digits=2, tol=0.1, agg=5,
                                     sparse=TRUE, grid_shape=c(2,2),
                                     sub_sample = 1, verbose=FALSE, check=TRUE)
  test <- cv_scores == cv_scores.mc
  expect_false(any(!test))

})

