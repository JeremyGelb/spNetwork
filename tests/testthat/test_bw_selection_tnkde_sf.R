context("testing the bandwidth selection functions for tnkde")

library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR BW SELECTION WITH CV LIKELIHOOD ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the bw selection function with CV likelihood and simple kernel", {
  ## creating the simple sf::st_as_sf
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3))
  event$Time <- c(5,7,6)
  event <- st_as_sf(event, coords = c("x","y"))


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
  total <- (loo1 + loo2 + loo3) / 3


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
                                           grid_shape=c(3,3),
                                           sub_sample=1,
                                           verbose=FALSE,
                                           check=FALSE)

  expect_equal(obs_value[1,1], total)
})



test_that("Testing the bw selection function with CV likelihood and simple kernel and a 0 density", {
  ## creating the simple sf::st_as_sf
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-5))
  event$Time <- c(5,7,6)
  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw_net of 10 here
  bw_net <- 7
  bw_time <- 6

  #at e1, the time distance to e2 is 2, with a bw_time of 6
  loo1 <- sum(log(
    (quartic_kernel(2,bw_time) * quartic_kernel(6,bw_net)) * (1/(bw_net*bw_time))
  ))

  loo2 <- sum(log(
    (quartic_kernel(2,bw_time) * quartic_kernel(6,bw_net)) * (1/(bw_net*bw_time))
  ))


  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  # but the last event is too far and we set zero_strat to "remove"
  total <- (loo1 + loo2) / 2


  #let us calculate the value with our function
  obs_value <- bws_tnkde_cv_likelihood_calc(bw_net_range = c(7,15),
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
                                           zero_strat = "remove",
                                           grid_shape=c(3,3),
                                           sub_sample=1,
                                           verbose=FALSE,
                                           check=FALSE)

  expect_equal(obs_value[1,1], total)
})



test_that("Testing the bw selection function with CV likelihood and discontinuous kernel", {
  ## creating the simple sf::st_as_sf
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3))
  event$Time <- c(5,7,6)
  event <- st_as_sf(event, coords = c("x","y"))


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
  total <- (loo1 + loo2 + loo3)/3


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
  ## creating the simple sf::st_as_sf
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)",
    "LINESTRING (-5 5, 0 5)",
    "LINESTRING (0 5, 5 5)"
    )

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-4))
  event$Time <- c(5,7,6)
  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw_net of 11 here
  bw_net <- 11
  bw_time <- 6

  # at e2 and e3, the network densities will be the same, there is no backfire at
  # the end of the lines
  knet1 <- quartic_kernel(6,bw_net)*(2/4) + (((-1/3)*(2/4)) * quartic_kernel(10,bw_net))
  knet2 <- quartic_kernel(7,bw_net)*(2/4)
  loo2 <-  sum(log(
    (quartic_kernel(2,bw_time) * (knet1) +
       quartic_kernel(1,bw_time) * (knet2)) * (1/(bw_net*bw_time))
  ))


  loo3 <-  sum(log(
    (quartic_kernel(1,bw_time) * (quartic_kernel(7,bw_net)*(2/4)) +
       quartic_kernel(1,bw_time) * (quartic_kernel(7,bw_net)*(2/4))) * (1/(bw_net*bw_time))
  ))

  #at e1, there is some backfire
  alpha1 <- 2/4
  alpha2 <- alpha1 * ((2-3)/3)

  net_kernel1 <- (quartic_kernel(6,bw_net)*(alpha1)) + alpha2 * quartic_kernel(10,bw_net)
  net_kernel2 <- (quartic_kernel(7,bw_net)*(alpha1))

  time_kernel1 <- quartic_kernel(2,bw_time)
  time_kernel2 <- quartic_kernel(1,bw_time)

  loo1 <- sum(log(
    ((time_kernel1 * net_kernel1) +
      (time_kernel2 * net_kernel2)) * (1/(bw_net*bw_time))
  ))



  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- (loo1 + loo2 + loo3)/3


  #let us calculate the value with our function
  obs_value <-bws_tnkde_cv_likelihood_calc(bw_net_range = c(11,12),
                                           bw_net_step = 1,
                                           bw_time_range = c(6,7),
                                           bw_time_step = 1,
                                           lines = all_lines,
                                           events = event,
                                           time_field = "Time",
                                           w = c(1,1,1),
                                           kernel_name = "quartic",
                                           method = "continuous",
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


test_that("Testing the bw selection function with CV likelihood in multicore and simple core", {
  ## creating the simple sf::st_as_sf
  # start with de definition of some lines

  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)",
    "LINESTRING (-5 5, 0 5)",
    "LINESTRING (0 5, 5 5)"
  )

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-4))
  event$Time <- c(5,7,6)
  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw_net of 11 here
  bw_net <- 11
  bw_time <- 6

  # at e2 and e3, the network densities will be the same, there is no backfire at
  # the end of the lines
  knet1 <- quartic_kernel(6,bw_net)*(2/4) + (((-1/3)*(2/4)) * quartic_kernel(10,bw_net))
  knet2 <- quartic_kernel(7,bw_net)*(2/4)
  loo2 <-  sum(log(
    (quartic_kernel(2,bw_time) * (knet1) +
       quartic_kernel(1,bw_time) * (knet2)) * (1/(bw_net*bw_time))
  ))


  loo3 <-  sum(log(
    (quartic_kernel(1,bw_time) * (quartic_kernel(7,bw_net)*(2/4)) +
       quartic_kernel(1,bw_time) * (quartic_kernel(7,bw_net)*(2/4))) * (1/(bw_net*bw_time))
  ))

  #at e1, there is some backfire
  alpha1 <- 2/4
  alpha2 <- alpha1 * ((2-3)/3)

  net_kernel1 <- (quartic_kernel(6,bw_net)*(alpha1)) + alpha2 * quartic_kernel(10,bw_net)
  net_kernel2 <- (quartic_kernel(7,bw_net)*(alpha1))

  time_kernel1 <- quartic_kernel(2,bw_time)
  time_kernel2 <- quartic_kernel(1,bw_time)

  loo1 <- sum(log(
    ((time_kernel1 * net_kernel1) +
       (time_kernel2 * net_kernel2)) * (1/(bw_net*bw_time))
  ))



  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- (loo1 + loo2 + loo3)/3


  #let us calculate the value with our function
  obs_value <-bws_tnkde_cv_likelihood_calc(bw_net_range = c(11,12),
                                           bw_net_step = 1,
                                           bw_time_range = c(6,7),
                                           bw_time_step = 1,
                                           lines = all_lines,
                                           events = event,
                                           time_field = "Time",
                                           w = c(1,1,1),
                                           kernel_name = "quartic",
                                           method = "continuous",
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

  future::plan(future::multisession(workers=2))
  obs_value2 <-bws_tnkde_cv_likelihood_calc.mc(bw_net_range = c(11,12),
                                           bw_net_step = 1,
                                           bw_time_range = c(6,7),
                                           bw_time_step = 1,
                                           lines = all_lines,
                                           events = event,
                                           time_field = "Time",
                                           w = c(1,1,1),
                                           kernel_name = "quartic",
                                           method = "continuous",
                                           diggle_correction = FALSE,
                                           study_area = NULL,
                                           max_depth = 15,
                                           digits=5,
                                           tol=0.001,
                                           agg=NULL,
                                           sparse=TRUE,
                                           grid_shape=c(3,3),
                                           sub_sample=1,
                                           verbose=TRUE,
                                           check=FALSE)

  test1 <- sum(obs_value - obs_value2) == 0
  test2 <- obs_value[1,1] == total
  expect_true(test1 & test2)
})

