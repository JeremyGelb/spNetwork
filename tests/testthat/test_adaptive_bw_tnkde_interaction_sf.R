context("testing functions for adptative bandwidth with interaction in TNKDE")
library(sf)

test_that("Testing the adaptive_bw_tnkde function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf,
                        wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1),
                      y=c(3,0,0),
                      id = 1:3,
                      times = c(0,0,2))

  event <- st_as_sf(event, coords = c("x","y"))

  # definition of one sample point
  sp_points <- data.frame(x=c(0,-1),
                      y=c(0,0),
                      id = 1:2)

  sp_points <- st_as_sf(sp_points, coords = c("x","y"))

  # the bandwidth will be adapted arround the events
  # based on the density at the event, so we must first
  # calculate that density. we have here a case of
  # discontinuous NKDE with bw_net = 3 et kernel quartic
  # et bw_time = 3

  n1 <- quartic_kernel(0,3) * quartic_kernel(0,3) #(spatial alone * temporal alone)

  n2 <- (quartic_kernel(0,3) * quartic_kernel(0,3)) +
    (quartic_kernel(2,3) * quartic_kernel(2,3))

  n3 <- (quartic_kernel(0,3) * quartic_kernel(0,3)) +
    (quartic_kernel(2,3) * quartic_kernel(2,3))

  hf0 <- c(n1,n2,n3)*(1/9)

  h0 <- 3
  gamma_val <- exp(sum(log(1/sqrt(hf0)))/3)
  abws_net <- h0 * (1/sqrt(hf0)) * (1/gamma_val)
  abws_time <- h0 * (1/sqrt(hf0)) * (1/gamma_val)

  observed <- tnkde(
    events = event,
    time_field = "times",
    w = rep(1,nrow(event)),
    samples_loc = sp_points,
    samples_time = c(1,2,3),
    lines = all_lines,
    adaptive = TRUE,
    adaptive_separate = FALSE,
    bw_net = 3, bw_time = 3,
    trim_bw_net = 150,
    trim_bw_time = 150,
    method = "discontinuous",
    kernel_name = "quartic",
    max_depth = 8,
    tol = 0.1,
    digits = 2,
    sparse = TRUE,
    verbose = TRUE,
    check = FALSE
  )
  diff1 <- sum(abs(abws_net - observed$events$bw_net))
  diff2 <- sum(abs(abws_time - observed$events$bw_time))

  diff <- round(sum(diff1,diff2),12)
  expect_equal(diff, 0)
})


test_that("Testing the adaptive_bw_tnkde.mc function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf,
                        wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1),
                      y=c(3,0,0),
                      id = 1:3,
                      times = c(0,0,2))

  event <- st_as_sf(event, coords = c("x","y"))

  # definition of two sample point
  sp_points <- data.frame(x=c(0,-1),
                          y=c(0,0),
                          id = 1:2)

  sp_points <- st_as_sf(sp_points, coords = c("x","y"))

  # the bandwidth will be adapted arround the events
  # based on the density at the event, so we must first
  # calculate that density. we have here a case of
  # continuous NKDE with bw_net = 3 et kernel quartic
  # et bw_time = 3

  n1 <- quartic_kernel(0,3) * quartic_kernel(0,3) #(spatial alone * temporal alone)

  n3 <- (quartic_kernel(0,3) * quartic_kernel(0,3)) +
    (quartic_kernel(2,3) * quartic_kernel(2,3))

  ## n2 is slightly harder

  nk1 <- quartic_kernel(0,3) # this is the base self-weight for network
  n <- 4
  nk2 <- quartic_kernel(2,3) * ((2.0-n)/n) # this is the backfire on network
  n2 <- ((nk1+nk2) * quartic_kernel(0,3)) + (quartic_kernel(2,3) * quartic_kernel(2,3))

  hf0 <- c(n1,n3,n2)*(1/9)

  h0 <- 3
  gamma_val <- exp(sum(log(1/sqrt(hf0)))/3)
  abws_net <- h0 * (1/sqrt(hf0)) * (1/gamma_val)
  abws_time <- h0 * (1/sqrt(hf0)) * (1/gamma_val)

  future::plan(future::multisession(workers=1))
  observed <- tnkde.mc(
    events = event,
    time_field = "times",
    w = rep(1,nrow(event)),
    samples_loc = sp_points,
    samples_time = c(1,2,3),
    lines = all_lines,
    adaptive = TRUE,
    adaptive_separate = FALSE,
    bw_net = 3, bw_time = 3,
    trim_bw_net = 150,
    trim_bw_time = 150,
    method = "continuous",
    kernel_name = "quartic",
    max_depth = 8,
    tol = 0.1,
    digits = 2,
    sparse = TRUE,
    verbose = TRUE,
    check = FALSE
  )
  diff1 <- sum(abs(abws_net - observed$events$bw_net))
  diff2 <- sum(abs(abws_time - observed$events$bw_time))

  diff <- round(sum(diff1,diff2),12)
  expect_equal(diff, 0)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TESTING THE ADAPTIVE BW FUNCTION WHEN SEVERAL BWS ARE GIVEN ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the adaptive_bw_tnkde function with multiple bws (simple)", {

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
  event$Time <- c(5,5,6,5)
  event <- st_as_sf(event, coords = c("x","y"))

  bws_net <- c(10, 15, 20)
  bws_time <- c(6,7)

  # we will calculate the expected kernel densities at each location
  # for the simple TNKDE

  bw_net <- 10
  bw_time <- 6

  ## location 1
  k1 <- quartic_kernel(0,bw_net) * quartic_kernel(0,bw_time)
  k1 <- k1 / (bw_net * bw_time)

  k2 <- quartic_kernel(5,bw_net) * quartic_kernel(0,bw_time)
  k2 <- (k2 * 2) / (bw_net * bw_time)

  k3 <- quartic_kernel(6,bw_net) * quartic_kernel(1,bw_time)
  k3 <- k3 / (bw_net * bw_time)

  dens1 <- k1 + k2 + k3


  ## location 2
  k1 <- quartic_kernel(0,bw_net) * quartic_kernel(0,bw_time)
  k1 <- (k1 * 2) / (bw_net * bw_time)

  k2 <- quartic_kernel(5,bw_net) * quartic_kernel(0,bw_time)
  k2 <- k2 / (bw_net * bw_time)

  k3 <- quartic_kernel(5,bw_net) * quartic_kernel(1,bw_time)
  k3 <- k3 / (bw_net * bw_time)

  dens2 <- k1 + k2 + k3

  ## location 3
  k1 <- quartic_kernel(0,bw_net) * quartic_kernel(0,bw_time)
  k1 <- (k1) / (bw_net * bw_time)

  k2 <- quartic_kernel(5,bw_net) * quartic_kernel(1,bw_time)
  k2 <- (k2 * 2) / (bw_net * bw_time)

  k3 <- quartic_kernel(6,bw_net) * quartic_kernel(1,bw_time)
  k3 <- k3 / (bw_net * bw_time)

  dens3 <- k1 + k2 + k3

  # la première valeur attendue avec les bw 10,6
  hf0 <- c(dens1, dens2, dens3, dens2)
  gamma_val <- calc_gamma(hf0)
  abws_net1 <- bw_net * (1/sqrt(hf0)) * (1/gamma_val)
  abws_time1 <- bw_time * (1/sqrt(hf0)) * (1/gamma_val)


  # deuxieme comparaison
  bw_net <- 15
  bw_time <- 7

  ## location 1
  k1 <- quartic_kernel(0,bw_net) * quartic_kernel(0,bw_time)
  k1 <- k1 / (bw_net * bw_time)

  k2 <- quartic_kernel(5,bw_net) * quartic_kernel(0,bw_time)
  k2 <- (k2 * 2) / (bw_net * bw_time)

  k3 <- quartic_kernel(6,bw_net) * quartic_kernel(1,bw_time)
  k3 <- k3 / (bw_net * bw_time)

  dens1 <- k1 + k2 + k3


  ## location 2
  k1 <- quartic_kernel(0,bw_net) * quartic_kernel(0,bw_time)
  k1 <- (k1 * 2) / (bw_net * bw_time)

  k2 <- quartic_kernel(5,bw_net) * quartic_kernel(0,bw_time)
  k2 <- k2 / (bw_net * bw_time)

  k3 <- quartic_kernel(5,bw_net) * quartic_kernel(1,bw_time)
  k3 <- k3 / (bw_net * bw_time)

  dens2 <- k1 + k2 + k3

  ## location 3
  k1 <- quartic_kernel(0,bw_net) * quartic_kernel(0,bw_time)
  k1 <- (k1) / (bw_net * bw_time)

  k2 <- quartic_kernel(5,bw_net) * quartic_kernel(1,bw_time)
  k2 <- (k2 * 2) / (bw_net * bw_time)

  k3 <- quartic_kernel(6,bw_net) * quartic_kernel(1,bw_time)
  k3 <- k3 / (bw_net * bw_time)

  dens3 <- k1 + k2 + k3

  # la première valeur attendue avec les bw 10,6
  hf0 <- c(dens1, dens2, dens3, dens2)
  gamma_val <- calc_gamma(hf0)
  abws_net2 <- bw_net * (1/sqrt(hf0)) * (1/gamma_val)
  abws_time2 <- bw_time * (1/sqrt(hf0)) * (1/gamma_val)

  ## on se lance dans la preparation des donnees brutes
  events <- event
  time_field <- "Time"
  w <- c(1,1,1,1)

  samples <- events

  events$time <- events[[time_field]]
  events$weight <- w
  div <- "bw"
  events$wid <- 1:nrow(events)
  bw_net_range <- c(10,20)
  bw_net_step <- 5
  bw_time_range <- c(6,7)
  bw_time_step <- 1
  agg <- NULL
  digits <- 2
  tol <- 0.001
  lines <- all_lines
  grid_shape <- c(3,3)

  data <- prepare_data(samples, lines, events, w , digits,tol,agg)
  lines <- data$lines
  events_loc <- data$events
  idx <- closest_points(events, events_loc)
  events$goid <- events_loc$goid[idx]
  grid <- build_grid(grid_shape,list(lines,samples,events))
  net_bws <- seq(from = bw_net_range[[1]], to = bw_net_range[[2]], by = bw_net_step)
  time_bws <- seq(from = bw_time_range[[1]], to = bw_time_range[[2]], by = bw_time_step)

  all_bws <- adaptive_bw_tnkde(grid = grid,
                               events_loc = events_loc,
                               events = events,
                               lines = lines,
                               bw_net = net_bws,
                               bw_time = time_bws,
                               trim_bw_net = net_bws * 2,
                               trim_bw_time = time_bws * 2,
                               method = "simple",
                               kernel_name = "quartic",
                               max_depth = 8,
                               div = "bw",
                               tol = tol,
                               digits = digits,
                               sparse = TRUE,
                               verbose = TRUE)

  # comparaison 1 :
  obt_net_1 <- all_bws[[1]][1,1,]
  obt_time_1 <- all_bws[[2]][1,1,]

  expect_equal(obt_net_1, abws_net1)
  expect_equal(obt_time_1, abws_time1)

  # comparaison 2
  obt_net_2 <- all_bws[[1]][2,2,]
  obt_time_2 <- all_bws[[2]][2,2,]

  expect_equal(obt_net_2, abws_net2)
  expect_equal(obt_time_2, abws_time2)

})



test_that("Testing the adaptive_bw_tnkde function with multiple bws (simple)", {

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
  event$Time <- c(5,5,6,5)
  event <- st_as_sf(event, coords = c("x","y"))

  bws_net <- c(10, 15, 20)
  bws_time <- c(6,7)

  # we will calculate the expected kernel densities at each location
  # for the simple TNKDE

  bw_net <- 10
  bw_time <- 6

  ## location 1
  k1 <- quartic_kernel(0,bw_net) * quartic_kernel(0,bw_time)
  k1 <- k1 / (bw_net * bw_time)

  k2 <- (quartic_kernel(5,bw_net) * (1/3)) * quartic_kernel(0,bw_time)
  k2 <- (k2 * 2) / (bw_net * bw_time)

  k3 <- (quartic_kernel(6,bw_net) * (1/3)) * quartic_kernel(1,bw_time)
  k3 <- k3 / (bw_net * bw_time)

  dens1 <- k1 + k2 + k3


  ## location 2
  k1 <- quartic_kernel(0,bw_net) * quartic_kernel(0,bw_time)
  k1 <- (k1 * 2) / (bw_net * bw_time)

  k2 <- (quartic_kernel(5,bw_net) * (1/3)) * quartic_kernel(0,bw_time)
  k2 <- k2 / (bw_net * bw_time)

  k3 <- (quartic_kernel(5,bw_net)* (1/3)) * quartic_kernel(1,bw_time)
  k3 <- k3 / (bw_net * bw_time)

  dens2 <- k1 + k2 + k3

  ## location 3
  k1 <- quartic_kernel(0,bw_net) * quartic_kernel(0,bw_time)
  k1 <- (k1) / (bw_net * bw_time)

  k2 <- (quartic_kernel(5,bw_net)* (1/3)) * quartic_kernel(1,bw_time)
  k2 <- (k2 * 2) / (bw_net * bw_time)

  k3 <- (quartic_kernel(6,bw_net)* (1/3)) * quartic_kernel(1,bw_time)
  k3 <- k3 / (bw_net * bw_time)

  dens3 <- k1 + k2 + k3

  # la première valeur attendue avec les bw 10,6
  hf0 <- c(dens1, dens2, dens3, dens2)
  gamma_val <- calc_gamma(hf0)
  abws_net1 <- bw_net * (1/sqrt(hf0)) * (1/gamma_val)
  abws_time1 <- bw_time * (1/sqrt(hf0)) * (1/gamma_val)


  ## on se lance dans la preparation des donnees brutes
  events <- event
  time_field <- "Time"
  w <- c(1,1,1,1)

  samples <- events

  events$time <- events[[time_field]]
  events$weight <- w
  div <- "bw"
  events$wid <- 1:nrow(events)
  bw_net_range <- c(10,20)
  bw_net_step <- 5
  bw_time_range <- c(6,7)
  bw_time_step <- 1
  agg <- NULL
  digits <- 2
  tol <- 0.001
  lines <- all_lines
  grid_shape <- c(3,3)

  data <- prepare_data(samples, lines, events, w , digits,tol,agg)
  lines <- data$lines
  events_loc <- data$events
  idx <- closest_points(events, events_loc)
  events$goid <- events_loc$goid[idx]
  grid <- build_grid(grid_shape,list(lines,samples,events))
  net_bws <- seq(from = bw_net_range[[1]], to = bw_net_range[[2]], by = bw_net_step)
  time_bws <- seq(from = bw_time_range[[1]], to = bw_time_range[[2]], by = bw_time_step)

  all_bws <- adaptive_bw_tnkde(grid = grid,
                               events_loc = events_loc,
                               events = events,
                               lines = lines,
                               bw_net = net_bws,
                               bw_time = time_bws,
                               trim_bw_net = net_bws * 2,
                               trim_bw_time = time_bws * 2,
                               method = "discontinuous",
                               kernel_name = "quartic",
                               max_depth = 8,
                               div = "bw",
                               tol = tol,
                               digits = digits,
                               sparse = TRUE,
                               verbose = TRUE)

  # comparaison 1 :
  obt_net_1 <- all_bws[[1]][1,1,]
  obt_time_1 <- all_bws[[2]][1,1,]

  expect_equal(obt_net_1, abws_net1)
  expect_equal(obt_time_1, abws_time1)


})



