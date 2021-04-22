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
  expect_equal(obs_value[1,2], real_value)
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
  expect_equal(obs_value[1,2], real_value)
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
  expect_equal(obs_value[1,2], real_value)
})
