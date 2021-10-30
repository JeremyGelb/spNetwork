context("testing the bandwidth selection functions for tnkde")

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
  event$Time <- c(5,7,6)
  sp::coordinates(event) <- cbind(event$x,event$y)


  # we can admit a bw_net of 10 here
  bw_net <- 10
  bw_time <- 6

  #at e1, the time distance to e2 is 2, with a bw_time of 6
  loo1 <- sum(log(
    (quartic_kernel(2,bw_time) * quartic_kernel(6,bw_net) +
    quartic_kernel(1,bw_time) * quartic_kernel(6,bw_net)) * (1/(bw_net*bw_time))
    ))

  loo2 <- sum(log(
    (quartic_kernel(2,bw_time) * quartic_kernel(6,bw_net) +
       quartic_kernel(1,bw_time) * quartic_kernel(6,bw_net)) * (1/(bw_net*bw_time))
     ))

  loo3 <- sum(log(
    (quartic_kernel(1,bw_time) * quartic_kernel(6,bw_net) +
       quartic_kernel(1,bw_time) * quartic_kernel(6,bw_net)) * (1/(bw_net*bw_time))
  ))


  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- loo1 + loo2 + loo3


  #let us calculate the value with our function
  obs_value <-bws_tnkde_cv_likelihood_calc(bw_net_range = c(10,15),
                                           bw_net_step = 5,
                                           bw_time_range = c(6,7),
                                           bw_time_step = 1,
                                           lines = all_lines,
                                           events = event,
                                           time_field = "Time",
                                           w = c(1,1,1),
                                           kernel_name = "quartic",
                                           method = "simple",
                                           diggle_correction = FALSE,
                                           study_area = NULL,
                                           max_depth = 15,
                                           digits=5,
                                           tol=0.001,
                                           agg=NULL,
                                           sparse=TRUE,
                                           grid_shape=c(1,1),
                                           sub_sample=1,
                                           verbose=TRUE,
                                           check=FALSE)

  expect_equal(obs_value[1,1], total)
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
  event$Time <- c(5,7,6)
  sp::coordinates(event) <- cbind(event$x,event$y)


  # we can admit a bw_net of 10 here
  bw_net <- 10
  bw_time <- 6

  #at e1, the time distance to e2 is 2, with a bw_time of 6
  loo1 <- sum(log(
    (quartic_kernel(2,bw_time) * (quartic_kernel(6,bw_net)*(1/3)) +
       quartic_kernel(1,bw_time) * (quartic_kernel(6,bw_net)*(1/3))) * (1/(bw_net*bw_time))
  ))

  loo2 <- sum(log(
    (quartic_kernel(2,bw_time) * (quartic_kernel(6,bw_net)*(1/3)) +
       quartic_kernel(1,bw_time) * (quartic_kernel(6,bw_net)*(1/3))) * (1/(bw_net*bw_time))
  ))

  loo3 <- sum(log(
    (quartic_kernel(1,bw_time) * (quartic_kernel(6,bw_net)*(1/3)) +
       quartic_kernel(1,bw_time) * (quartic_kernel(6,bw_net)*(1/3))) * (1/(bw_net*bw_time))
  ))


  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- loo1 + loo2 + loo3


  #let us calculate the value with our function
  obs_value <-bws_tnkde_cv_likelihood_calc(bw_net_range = c(10,15),
                                           bw_net_step = 5,
                                           bw_time_range = c(6,7),
                                           bw_time_step = 1,
                                           lines = all_lines,
                                           events = event,
                                           time_field = "Time",
                                           w = c(1,1,1),
                                           kernel_name = "quartic",
                                           method = "discontinuous",
                                           diggle_correction = FALSE,
                                           study_area = NULL,
                                           max_depth = 15,
                                           digits=5,
                                           tol=0.001,
                                           agg=NULL,
                                           sparse=TRUE,
                                           grid_shape=c(1,1),
                                           sub_sample=1,
                                           verbose=TRUE,
                                           check=FALSE)

  expect_equal(obs_value[1,1], total)
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
  event$Time <- c(5,7,6)
  sp::coordinates(event) <- cbind(event$x,event$y)


  # we can admit a bw_net of 10 here
  bw_net <- 12
  bw_time <- 6

  #at e1, the time distance to e2 is 2, with a bw_time of 6
  # calculating the kernel density provoqued by an other event
  n <- 4

  direct_effect <- quartic_kernel(6,bw_net)*(2/n)
  back_fire <- quartic_kernel(10, bw_net) * ()

  loo1 <- sum(log(
    (quartic_kernel(2,bw_time) * (quartic_kernel(6,bw_net)*(1/3)) +
       quartic_kernel(1,bw_time) * (quartic_kernel(6,bw_net)*(1/3))) * (1/(bw_net*bw_time))
  ))

  loo2 <- sum(log(
    (quartic_kernel(2,bw_time) * (quartic_kernel(6,bw_net)*(1/3)) +
       quartic_kernel(1,bw_time) * (quartic_kernel(6,bw_net)*(1/3))) * (1/(bw_net*bw_time))
  ))

  loo3 <- sum(log(
    (quartic_kernel(1,bw_time) * (quartic_kernel(6,bw_net)*(1/3)) +
       quartic_kernel(1,bw_time) * (quartic_kernel(6,bw_net)*(1/3))) * (1/(bw_net*bw_time))
  ))


  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- loo1 + loo2 + loo3


  #let us calculate the value with our function
  obs_value <-bws_tnkde_cv_likelihood_calc(bw_net_range = c(10,15),
                                           bw_net_step = 5,
                                           bw_time_range = c(6,7),
                                           bw_time_step = 1,
                                           lines = all_lines,
                                           events = event,
                                           time_field = "Time",
                                           w = c(1,1,1),
                                           kernel_name = "quartic",
                                           method = "discontinuous",
                                           diggle_correction = FALSE,
                                           study_area = NULL,
                                           max_depth = 15,
                                           digits=5,
                                           tol=0.001,
                                           agg=NULL,
                                           sparse=TRUE,
                                           grid_shape=c(1,1),
                                           sub_sample=1,
                                           verbose=TRUE,
                                           check=FALSE)

  expect_equal(obs_value[1,1], total)
})



