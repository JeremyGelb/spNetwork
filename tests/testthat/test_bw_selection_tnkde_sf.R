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
  obs_value <- bw_tnkde_cv_likelihood_calc(bws_net = seq(10,15,5),
                                           bws_time = seq(6,7,1),
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
  obs_value <- bw_tnkde_cv_likelihood_calc(bws_net = seq(7,15,5),
                                           bws_time = seq(6,7,1),
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
  obs_value <-bw_tnkde_cv_likelihood_calc(bws_net = seq(10,15,5),
                                           bws_time = seq(6,7,1),
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
  obs_value <-bw_tnkde_cv_likelihood_calc(bws_net = seq(11,12,1),
                                           bws_time = seq(6,7,1),
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
  obs_value <- bw_tnkde_cv_likelihood_calc(bws_net = seq(11,12,1),
                                           bws_time = seq(6,7,1),
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
  obs_value2 <-bw_tnkde_cv_likelihood_calc.mc(bws_net = seq(11,12,1),
                                           bws_time = seq(6,7,1),
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR BW SELECTION WITH CV LIKELIHOOD WITH ADAPTIVE BW ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## TEST FONCTIONNEL POUR LA VERSION SIMPLE
test_that("Testing the bw selection function with CV likelihood and simple kernel AND adaptive BW", {
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
  event <- data.frame(x=c(0,2,0,2),
                      y=c(3,0,-3,0))
  event$Time <- c(5,5,6,6)
  event <- st_as_sf(event, coords = c("x","y"))

  # tm_shape(all_lines) + tm_lines('black') + tm_shape(events) + tm_dots('red', size = 0.5)

  bws_net <- c(10, 15, 20)
  bws_time <- c(6,7)

  # first : for the pair of bws 10,6
  # I need to calculate the local densities
  # and use it to caculate the local bandwidths
  bw_time <- 6
  bw_net <- 10
  net_dist_mat <- rbind(
    c(0,5,6,5),
    c(5,0,5,0),
    c(6,5,0,5),
    c(5,0,5,0)
  )
  time_dist_mat <- rbind(
    c(0,0,1,1),
    c(0,0,1,1),
    c(1,1,0,0),
    c(1,1,0,0)
  )
  n <- nrow(time_dist_mat)
  hf0 <- sapply(1:n, function(i){

    # ids <- setdiff(c(1:n),i)
    #
    # sum(quartic_kernel(net_dist_mat[i,ids],10) *
    #       quartic_kernel(time_dist_mat[i,ids],6))

    sum(quartic_kernel(net_dist_mat[i,],bw_net) *
      quartic_kernel(time_dist_mat[i,],bw_time))

  }) * (1/(10*6))

  gamma_val <- calc_gamma(hf0)
  abws_net <- bw_net * (1/sqrt(hf0)) * (1/gamma_val)
  abws_time <- bw_time * (1/sqrt(hf0)) * (1/gamma_val)

  # second, I can use the local densities to calculate the loo values
  loo_values <- sapply(1:n, function(i){
    ids <- setdiff(c(1:n),i)
    sum( ((quartic_kernel(net_dist_mat[i,ids],abws_net[ids])) *
            quartic_kernel(time_dist_mat[i,ids],abws_time[ids])) * (1/(abws_net[ids]*abws_time[ids]))

    )
  })
  loo_value <- sum(log(loo_values)) / n

  #let us calculate the value with our function
  obs_value <- bw_tnkde_cv_likelihood_calc(bws_net = seq(10,20,5),
                                          bws_time = seq(6,7,1),
                                          lines = all_lines,
                                          events = event,
                                          time_field = "Time",
                                          w = c(1,1,1,1),
                                          kernel_name = "quartic",
                                          method = "simple",
                                          diggle_correction = FALSE,
                                          study_area = NULL,
                                          adaptive = TRUE,
                                          trim_net_bws = c(20,30,40),
                                          trim_time_bws = c(12,14),
                                          max_depth = 15,
                                          digits=5,
                                          tol=0.001,
                                          agg=NULL,
                                          sparse=TRUE,
                                          grid_shape=c(1,1),
                                          sub_sample=1,
                                          verbose=FALSE,
                                          check=FALSE)

  expect_equal(obs_value[1,1], loo_value)
})

## TEST POUR LA VERSION DISCONTINUE
test_that("Testing the bw selection function with CV likelihood and discontinuous kernel AND adaptive BW", {
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
  event <- data.frame(x=c(0,2,0,2),
                      y=c(3,0,-3,0))
  event$Time <- c(5,5,6,6)
  event <- st_as_sf(event, coords = c("x","y"))

  # tm_shape(all_lines) + tm_lines('black') + tm_shape(events) + tm_dots('red', size = 0.5)

  bws_net <- c(10, 15, 20)
  bws_time <- c(6,7)

  # first : for the pair of bws 10,6
  # I need to calculate the local densities
  # and use it to caculate the local bandwidths
  bw_time <- 6
  bw_net <- 10
  net_dist_mat <- rbind(
    c(0,5,6,5),
    c(5,0,5,0),
    c(6,5,0,5),
    c(5,0,5,0)
  )
  time_dist_mat <- rbind(
    c(0,0,1,1),
    c(0,0,1,1),
    c(1,1,0,0),
    c(1,1,0,0)
  )

  alpha_mat <- rbind(
    c(1,3,3,3),
    c(3,1,3,1),
    c(3,3,1,3),
    c(3,1,3,1)
  )

  n <- nrow(time_dist_mat)
  hf0 <- sapply(1:n, function(i){

    sum( (quartic_kernel(net_dist_mat[i,],bw_net) / alpha_mat[i,]) *
          quartic_kernel(time_dist_mat[i,],bw_time))

  }) * (1/(10*6))

  gamma_val <- calc_gamma(hf0)
  abws_net <- bw_net * (1/sqrt(hf0)) * (1/gamma_val)
  abws_time <- bw_time * (1/sqrt(hf0)) * (1/gamma_val)

  # second, I can use the local densities to calculate the loo values

  loo_values <- sapply(1:n, function(i){
    ids <- setdiff(c(1:n),i)
    sum( ((quartic_kernel(net_dist_mat[i,ids],abws_net[ids]) /  alpha_mat[i,ids]) *
          quartic_kernel(time_dist_mat[i,ids],abws_time[ids])) * (1/(abws_net[ids]*abws_time[ids]))

         )
  })
  loo_value <- sum(log(loo_values)) / n

  #let us calculate the value with our function
  obs_value <- bw_tnkde_cv_likelihood_calc(bws_net = seq(10,20,5),
                                           bws_time = seq(6,7,1),
                                           lines = all_lines,
                                           events = event,
                                           time_field = "Time",
                                           w = c(1,1,1,1),
                                           kernel_name = "quartic",
                                           method = "discontinuous",
                                           diggle_correction = FALSE,
                                           study_area = NULL,
                                           adaptive = TRUE,
                                           trim_net_bws = c(20,30,40),
                                           trim_time_bws = c(12,14),
                                           max_depth = 15,
                                           digits=5,
                                           tol=0.001,
                                           agg=NULL,
                                           sparse=TRUE,
                                           grid_shape=c(1,1),
                                           sub_sample=1,
                                           verbose=FALSE,
                                           check=FALSE)

  expect_equal(obs_value[1,1], loo_value)
})

## TEST POUR LA VERSION CONTINUE
# PAS TERMINE
test_that("Testing the bw selection function with CV likelihood and continuous kernel AND adaptive BW", {
  ## creating the simple sf::st_as_sf
  # start with de definition of some lines

  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)"
    # "LINESTRING (-5 5, 0 5)",
    # "LINESTRING (-5 -5, 0 -5)",
    # "LINESTRING (0 5, 5 5)",
    # "LINESTRING (0 -5, 5 -5)"
    )

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,2,0,2),
                      y=c(3,0,-3,0))
  event$Time <- c(5,5,6,6)
  event <- st_as_sf(event, coords = c("x","y"))

  # tm_shape(all_lines) + tm_lines('black') + tm_shape(events) + tm_dots('red', size = 0.5)

  bws_net <- c(10, 15, 20)
  bws_time <- c(6,7)

  time_dist_mat <- rbind(
    c(0,0,1,1),
    c(0,0,1,1),
    c(1,1,0,0),
    c(1,1,0,0)
  )

  # first : for the pair of bws 10,6
  # I need to calculate the local densities
  # and use it to caculate the local bandwidths
  bw_time <- 6
  bw_net <- 10
  # I will have to calculate everything by hand because vector based
  # calculating is not possible in this complex case

  # for the first event
  dd1 <- quartic_kernel(0, bw_net) # self effect
  dd1 <- dd1 + (quartic_kernel(6,bw_net) * (-2/4))# backfire from below
  dd2 <- quartic_kernel(5, bw_net) * (2/4) # effect of event 2
  dd3 <- quartic_kernel(6, bw_net) * (2/4) # effect of event 3
  dd4 <- quartic_kernel(5, bw_net) * (2/4) # effect of event 4
  dens1_net <- c(dd1, dd2, dd3, dd4)

  dens1_time <- c(quartic_kernel(time_dist_mat[1,], bw_time))
  dens1 <- sum((dens1_net * dens1_time)) * (1/(bw_time*bw_net))

  # for the second event
  dens2_net <- c(
    (quartic_kernel(5, bw_net) * (2/4)),  # effect from event 1
    quartic_kernel(0, bw_net) + # self direct density
                   ( (-2/4) * (quartic_kernel(4, bw_net))), #backfire from center
    (quartic_kernel(5, bw_net) * (2/4)),  # effect from event 3
     quartic_kernel(0, bw_net) + # effect from event 4 (same location)
       ( (-2/4) * (quartic_kernel(4, bw_net))) #backfire from center
  )

  dens2_time <- c(quartic_kernel(time_dist_mat[2,], bw_time))
  dens2 <- sum((dens2_net * dens2_time)) * (1/(bw_time*bw_net))

  # for the third event
  dens3_net <- c(
    (quartic_kernel(6, bw_net) * (2/4)),  # effect from event 1
    (quartic_kernel(5, bw_net) * (2/4)),  # effect from event 2
    (quartic_kernel(0, bw_net) + # self direct density
                   ( (-2/4) * (quartic_kernel(6, bw_net)))), #backfire from center
    (quartic_kernel(5, bw_net) * (2/4))  # effect from event 4
  )

  dens3_time <- c(quartic_kernel(time_dist_mat[3,], bw_time))
  dens3 <- sum((dens3_net * dens3_time)) * (1/(bw_time*bw_net))

  # for the fourth event
  dens4_net <- c(
    (quartic_kernel(5, bw_net) * (2/4)),  # effect from event 1
    (quartic_kernel(0, bw_net) + # effect from event 2 (same location)
       ( (-2/4) * (quartic_kernel(4, bw_net)))), #backfire from center
    (quartic_kernel(5, bw_net) * (2/4)),  # effect from event 3
    (quartic_kernel(0, bw_net) + # self direct density
                   ( (-2/4) * (quartic_kernel(4, bw_net)))) #backfire from center
  )

  dens4_time <- c(quartic_kernel(time_dist_mat[4,], bw_time))
  dens4 <- sum((dens4_net * dens4_time)) *  (1/(bw_time*bw_net))

  hf0 <- c(dens1, dens2, dens3, dens4)

  n <- nrow(time_dist_mat)

  gamma_val <- calc_gamma(hf0)
  abws_net <- bw_net * (1/sqrt(hf0)) * (1/gamma_val)
  abws_time <- bw_time * (1/sqrt(hf0)) * (1/gamma_val)

  # second, I can use the local bandwidths to calculate the loo values
  # for the first event
  dd2 <- quartic_kernel(5, abws_net[[2]]) * (2/4) # effect of event 2
  dd3 <- quartic_kernel(6, abws_net[[3]]) * (2/4) # effect of event 3
  dd4 <- quartic_kernel(5, abws_net[[4]]) * (2/4) # effect of event 4
  dens1_net <- c(dd2, dd3, dd4)

  dens1_time <- c(quartic_kernel(time_dist_mat[1,c(2,3,4)], abws_time[c(2,3,4)]))
  dens1 <- sum((dens1_net * dens1_time) * (1/(abws_net[c(2,3,4)]*abws_time[c(2,3,4)])))

  # for the second event
  dens2_net <- c(
    (quartic_kernel(5, abws_net[[1]]) * (2/4)),  # effect from event 1
    (quartic_kernel(5, abws_net[[3]]) * (2/4)),  # effect from event 3
    quartic_kernel(0, abws_net[[4]]) + # effect from event 4 (same location)
      ( (-2/4) * (quartic_kernel(4, abws_net[[2]]))) #backfire from center
  )

  dens2_time <- c(quartic_kernel(time_dist_mat[2,c(1,3,4)], abws_time[c(1,3,4)]))
  dens2 <- sum((dens2_net * dens2_time) * (1/(abws_net[c(1,3,4)]*abws_time[c(1,3,4)])))

  # for the third event
  dens3_net <- c(
    (quartic_kernel(6, abws_net[[1]]) * (2/4)),  # effect from event 1
    (quartic_kernel(5, abws_net[[2]]) * (2/4)),  # effect from event 2
    (quartic_kernel(5, abws_net[[4]]) * (2/4))  # effect from event 4
  )

  dens3_time <- c(quartic_kernel(time_dist_mat[3,c(1,2,4)], abws_time[c(1,2,4)]))
  dens3 <- sum((dens3_net * dens3_time) * (1/(abws_time[c(1,2,4)]*abws_net[c(1,2,4)])))

  # for the fourth event
  dens4_net <- c(
    (quartic_kernel(5, abws_net[[1]]) * (2/4)),  # effect from event 1
    (quartic_kernel(0, abws_net[[2]]) + # effect from event 2 (same location)
       ( (-2/4) * (quartic_kernel(4, abws_net[[4]])))), #backfire from center
    (quartic_kernel(5, abws_net[[3]]) * (2/4))  # effect from event 3
  )

  dens4_time <- c(quartic_kernel(time_dist_mat[4,c(1,2,3)], abws_time[c(1,2,3)]))
  dens4 <- sum((dens4_net * dens4_time) *  (1/(abws_time[c(1,2,3)]*abws_net[c(1,2,3)])))




  loo_value <- sum(log(c(dens1, dens2, dens3, dens4))) / 4

  #let us calculate the value with our function
  # TODO : MAKE THE TEST WORK WITH A GRID !
  obs_value <- bw_tnkde_cv_likelihood_calc(bws_net = seq(10,20,5),
                                           bws_time = seq(6,7,1),
                                           lines = all_lines,
                                           events = event,
                                           time_field = "Time",
                                           w = c(1,1,1,1),
                                           kernel_name = "quartic",
                                           method = "continuous",
                                           diggle_correction = FALSE,
                                           study_area = NULL,
                                           adaptive = TRUE,
                                           trim_net_bws = c(20,30,40),
                                           trim_time_bws = c(12,14),
                                           max_depth = 15,
                                           digits=5,
                                           tol=0.001,
                                           agg=NULL,
                                           sparse=TRUE,
                                           grid_shape=c(3,3),
                                           sub_sample=1,
                                           verbose=FALSE,
                                           check=FALSE)

  expect_equal(obs_value[1,1], loo_value)
})



