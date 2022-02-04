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

  # definition of one sample point
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

  hf0 <- c(n1,n3, n2)*(1/9)

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


