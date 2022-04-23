context("testing the bandwidth selection functions")
library(sf)

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

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3),
                      id = c(1,2,3))

  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  s1 <- (quartic_kernel(6,10) + quartic_kernel(6,10)) *(1/10)

  #so the CV likelihood is the mean of the loo
  total <- log(s1)


  #let us calculate the value with our function
  obs_value <- bw_cv_likelihood_calc(bw_range = c(8,10),1,
                                     lines = all_lines,
                                     events = event, w = c(1,1,1),
                                     check = F,
                                     kernel_name = "quartic",
                                     method = "simple",
                                     digits = 1,
                                     agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[3,2], total)
})


test_that("Testing the bw selection function with CV likelihood and simple kernel and a 0 density discarded", {
  ## creating the simple situation

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
                      y=c(3,0,-5),
                      id = c(1,2,3))

  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw of 7 here
  # the network distance between two points is 6
  bw <- 7
  s1 <- quartic_kernel(6,bw) *(1/bw)
  s2 <- quartic_kernel(6,bw) *(1/bw)

  # NOTE s3 is discared because too far

  #so the CV likelihood is the mean of the loo
  total <- (log(s1) + log(s2)) / 2


  #let us calculate the value with our function
  obs_value <- bw_cv_likelihood_calc(bw_range = c(6,8),1,
                                     lines = all_lines,
                                     events = event, w = c(1,1,1),
                                     check = F,
                                     kernel_name = "quartic",
                                     method = "simple",
                                     digits = 1, zero_strat = "remove",
                                     agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[2,2], total)
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


  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  # for each point, we split the kernel at the intersection (1/3)
  s1 <- (quartic_kernel(6,10) * 1/3 + quartic_kernel(6,10) * 1/3) *(1/10)

  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- log(s1)


  #let us calculate the value with our function
  obs_value <- bw_cv_likelihood_calc(c(8,10),1,
                    lines = all_lines,
                    events = event, w = c(1,1,1),
                    check = F,
                    kernel_name = "quartic",
                    method = "discontinuous",
                    digits = 1,
                    agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[3,2], total)
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

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  # for each point, we split the kernel at the intersection and their is some backfire
  # lets calculate the impact of en event on another
  n <- 4
  imp <- quartic_kernel(6,10) * (2.0/n)

  s1 <- (imp + imp) * (1/10)

  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- log(s1)


  #let us calculate the value with our function
  obs_value <- bw_cv_likelihood_calc(c(8,10),1,
                                     lines = all_lines,
                                     events = event, w = c(1,1,1),
                                     check = F,
                                     kernel_name = "quartic",
                                     method = "continuous",
                                     digits = 1,
                                     agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[3,2], total)
})



test_that("Testing the bw selection function with CV likelihood and continuous kernel", {
  ## creating the simple situation
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
  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw_net of 11 here
  bw_net <- 11

  # at e2 and e3, the network densities will be the same, there is no backfire at
  # the end of the lines
  knet1 <- quartic_kernel(6,bw_net)*(2/4) + (((-1/3)*(2/4)) * quartic_kernel(10,bw_net))
  knet2 <- quartic_kernel(7,bw_net)*(2/4)
  loo2 <-  sum(log(
    ((knet1) +(knet2)) * (1/(bw_net))
  ))


  loo3 <-  sum(log(
    ((quartic_kernel(7,bw_net)*(2/4)) +
       (quartic_kernel(7,bw_net)*(2/4))) * (1/(bw_net))
  ))

  #at e1, there is some backfire
  alpha1 <- 2/4
  alpha2 <- alpha1 * ((2-3)/3)

  net_kernel1 <- (quartic_kernel(6,bw_net)*(alpha1)) + alpha2 * quartic_kernel(10,bw_net)
  net_kernel2 <- (quartic_kernel(7,bw_net)*(alpha1))

  loo1 <- sum(log(
    ((net_kernel1) +
       (net_kernel2)) * (1/(bw_net))
  ))



  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  total <- (loo1 + loo2 + loo3) / 3


  #let us calculate the value with our function
  obs_value <- bw_cv_likelihood_calc(bw_range = c(11,12),
                                           bw_step = 1,
                                           lines = all_lines,
                                           events = event,
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

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  # so the density at one event is
  s1 <- (quartic_kernel(6,10) + quartic_kernel(6,10) + quartic_kernel(0,10)) *(1/10)

  #so the score value is
  Wl <- sum(as.numeric(st_length(all_lines)))
  score <- ((1/s1+1/s1+1/s1) - Wl)**2


  #let us calculate the value with our function
  obs_value <- bw_cvl_calc(bw_range = c(8,10),
                            bw_step = 1,
                            lines = all_lines,
                            events = event,
                            w = c(1,1,1),
                            check = F,
                            kernel_name = "quartic",
                            method = "simple",
                            digits = 1,
                            agg = NULL,
                            verbose = F,
                            tol = 0.00001
  )
  expect_equal(obs_value[3,2], score)
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

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw of 10 here
  # the network distance between two points is 6
  # for each point, we split the kernel at the intersection (1/3)
  # so the density at one event is
  s1 <- (quartic_kernel(6,10)*1/3 + quartic_kernel(6,10)*1/3 + quartic_kernel(0,10)) *(1/10)

  #so the score value is
  Wl <- sum(as.numeric(st_length(all_lines)))
  score <- (Wl - ((1/s1)*3))**2


  #let us calculate the value with our function
  obs_value <- bw_cvl_calc(c(8,10),1,
                           lines = all_lines,
                           events = event, w = c(1,1,1),
                           check = F,
                           kernel_name = "quartic",
                           method = "discontinuous",
                           digits = 1,
                           agg = NULL, verbose = F,tol = 0.00001
  )
  expect_equal(obs_value[3,2], score)
})



test_that("Testing the bw selection function with Van Lieshout's Criterion and continuous kernel", {
  ## creating the simple situation
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
                      y=c(3,0,-4),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))


  # we can admit a bw_net of 11 here
  bw_net <- 11

  # calculating the density as e2
  knet0 <- quartic_kernel(0, bw_net) + (((2-4)/4) *  quartic_kernel(6, bw_net)) # self density + self backfire

  knet1 <- quartic_kernel(6,bw_net)*(2/4) + # the direct impact
    (((-1/3)*(2/4)) * quartic_kernel(10,bw_net)) # the impact after coming back from top

  knet2 <- quartic_kernel(7,bw_net)*(2/4) # the direct impact

  dens2 <-  sum(
    ((knet1) +(knet2) + knet0) * (1/(bw_net))
  )


  # calculating the density as e3
  knet0 <- quartic_kernel(0, bw_net) + (((2-4)/4) *  quartic_kernel(8, bw_net)) # self density + self backfire
  knet1 <- quartic_kernel(7,bw_net)*(2/4) # the direct impact, backfire is too far
  knet2 <- quartic_kernel(7,bw_net)*(2/4) # the direct impact, backfire is too far

  dens3 <-  sum(
    (knet0 + knet1 + knet2) * (1/(bw_net))
  )

  #at e1, there is some backfire

  knet0 <- quartic_kernel(0, bw_net) + # self direct
    (((2-4)/4) *  quartic_kernel(6, bw_net)) + # return after center
    (((2-3)/3) *  quartic_kernel(4, bw_net)) + # return after top
    (((2-3)/3) * ((2-4)/4) *  quartic_kernel(10, bw_net))+ # return after top and center (from top)
  (((2-3)/3) * ((2-4)/4) *  quartic_kernel(10, bw_net)) # return after top and center (from center)

  knet1 <- quartic_kernel(6,bw_net) * (2/4) + # direct impact
    quartic_kernel(10,bw_net) * (2/4) * ((2-3)/3) # and backfire

  knet2 <- quartic_kernel(7,bw_net) * (2/4) # direct impact


  dens1 <- sum(
    (knet0 + knet1 +knet2) * (1/(bw_net))
  )



  #so the CV likelihood is the sum for each point loo
  #in other words, three times the sum of two kerneffect
  wl <- sum(as.numeric(st_length(all_lines)))
  s1 <- sum(c(dens1,dens2,dens3)**(-1.0))
  total <- (s1 - wl) * (s1-wl)


  #let us calculate the value with our function
  obs_value <- bw_cvl_calc(c(8,11),1,
                            lines = all_lines,
                            events = event, w = c(1,1,1),
                            check = F,
                            kernel_name = "quartic",
                            method = "continuous",
                            digits = 1,
                            zero_strat = "min_double",
                            agg = NULL, verbose = TRUE,
                            tol = 0.00001
  )

  expect_equal(obs_value[4,2], total)
})




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Comparing multicore and single core ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing that bw selection by cv-likelihood gives the same score in single and multicore", {

  data(mtl_network)
  data(bike_accidents)

  case1 <- bike_accidents[25,]
  buff <- st_buffer(case1, dist = 1000)

  inter_idx1 <-  sf::st_intersects(bike_accidents, buff, sparse = FALSE)[,1]
  inter_idx2 <-  sf::st_intersects(mtl_network, buff, sparse = FALSE)[,1]

  # subset
  events <- subset(bike_accidents, inter_idx1)
  lines <- subset(mtl_network, inter_idx2)

  ## multicore cv score
  future::plan(future::multisession(workers=2))

  cv_scores.mc <- bw_cv_likelihood_calc.mc(c(200,400),100,
                                 lines, events,
                                 rep(1,nrow(events)),
                                 "scaled gaussian", "discontinuous",
                                 diggle_correction = FALSE, study_area = NULL,
                                 max_depth = 8,
                                 digits=2, tol=0.1, agg=5,
                                 sparse=TRUE, grid_shape=c(2,2),
                                 sub_sample = 1, verbose=FALSE, check=TRUE)

  cv_scores.mc2 <- bw_cv_likelihood_calc.mc(c(200,400),100,
                                           lines, events,
                                           rep(1,nrow(events)),
                                           "scaled gaussian", "discontinuous",
                                           diggle_correction = FALSE, study_area = NULL,
                                           max_depth = 8,
                                           digits=2, tol=0.1, agg=5,
                                           sparse=TRUE, grid_shape=c(2,2),
                                           sub_sample = 1, verbose=TRUE, check=TRUE)

  ## single core cv score
  cv_scores <- bw_cv_likelihood_calc(bw_range = c(200,400),
                                     bw_step = 100,
                                     lines = lines,
                                     events = events,
                                     w = rep(1,nrow(events)),
                                     kernel_name = "scaled gaussian",
                                     method = "discontinuous",
                                     diggle_correction = FALSE,
                                     study_area = NULL,
                                     max_depth = 8,
                                     digits=2,
                                     tol=0.1,
                                     agg=5,
                                     sparse=TRUE,
                                     grid_shape=c(2,2),
                                     sub_sample = 1,
                                     verbose=TRUE,
                                     check=TRUE)
  test1 <- cv_scores == cv_scores.mc
  test2 <- cv_scores.mc2 == cv_scores.mc
  expect_false(any(!test1) | any(!test2))

})


test_that("Testing that bw selection with Van Lieshout's Criterion gives the same score in single and multicore", {

  data(mtl_network)
  data(bike_accidents)

  case1 <- bike_accidents[25,]
  buff <- st_buffer(case1, dist = 1000)

  inter_idx1 <-  sf::st_intersects(bike_accidents, buff, sparse = FALSE)[,1]
  inter_idx2 <-  sf::st_intersects(mtl_network, buff, sparse = FALSE)[,1]

  # subset
  events <- subset(bike_accidents, inter_idx1)
  lines <- subset(mtl_network, inter_idx2)

  ## multicore cv score
  future::plan(future::multisession(workers=1))
  cv_scores.mc <- bw_cvl_calc.mc(c(200,400),100,
                                           lines, events,
                                           rep(1,nrow(events)),
                                           "gaussian", "discontinuous",
                                           diggle_correction = FALSE, study_area = NULL,
                                           max_depth = 8,
                                           digits=2, tol=0.1, agg=5,
                                           sparse=TRUE, grid_shape=c(2,2),
                                           sub_sample = 1, verbose=FALSE, check=TRUE)

  cv_scores.mc2 <- bw_cvl_calc.mc(c(200,400),100,
                                 lines, events,
                                 rep(1,nrow(events)),
                                 "gaussian", "discontinuous",
                                 diggle_correction = FALSE, study_area = NULL,
                                 max_depth = 8,
                                 digits=2, tol=0.1, agg=5,
                                 sparse=TRUE, grid_shape=c(2,2),
                                 sub_sample = 1, verbose=FALSE, check=TRUE)

  ## single core cv score
  cv_scores <- bw_cvl_calc(c(200,400),100,
                                     lines, events,
                                     rep(1,nrow(events)),
                                     "gaussian", "discontinuous",
                                     diggle_correction = FALSE, study_area = NULL,
                                     max_depth = 8,
                                     digits=2, tol=0.1, agg=5,
                                     sparse=TRUE, grid_shape=c(1,1),
                                     sub_sample = 1, verbose=FALSE, check=TRUE)

  diff1 <- cv_scores[,2]/100000 -  cv_scores.mc[,2]/100000
  s1 <- sum(abs(round(diff1)))
  diff2 <- cv_scores[,2]/100000 -  cv_scores.mc2[,2]/100000
  s2 <- sum(abs(round(diff2)))
  print(cv_scores)
  print(cv_scores.mc)
  print(cv_scores.mc2)

  test1 <- s1 == 0
  test2 <- s2 == 0

  expect_true(test1 & test2)

})

