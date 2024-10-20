context("testing functions for adptative bandwidth in NKDE")
library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE adaptive_bw FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the adaptive_bw function", {

  #skip_on_cran()
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
                      w = c(1,1,2), # the last event has a weight of 2
                      id = 1:3)
  event <- st_as_sf(event, coords = c("x","y"))

  # the bandwidth will be adapted arround the events
  # based on the density at the event, so we must first
  # calculate that density. we have here a case of
  # discontinuous NKDE with bw = 3 et kernel quartic

  n1 <- quartic_kernel(0,3)
  n3 <- quartic_kernel(0,3)*2 + quartic_kernel(2,3)
  n2 <- quartic_kernel(0,3) + quartic_kernel(2,3)*2

  hf0 <- c(n1,n2,n3)

  h0 <- 3
  gamma_val <- exp(sum(log(1/sqrt(hf0)))/3)
  abws <- h0 * (1/sqrt(hf0)) * (1/gamma_val)

  grid <- build_grid(c(1,1),list(all_lines))

  observed <- nkde(
    events = event,
    w = event$w,
    samples = event,
    lines = all_lines,
    adaptive = TRUE,
    bw = 3,
    trim_bw = 5,
    method = "discontinuous",
    kernel_name = "quartic",
    max_depth = 8,
    tol = 0.1,
    digits = 2,
    sparse = TRUE,
    verbose = TRUE,
    check = FALSE,
    grid_shape = c(1,1)
  )

  diff <- sum(round(abs(abws - observed$events$bw),6))
  expect_equal(diff, 0)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE adaptive_bw.mc FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the adaptive_bw.mc function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1),
                      y=c(3,0,0),
                      w = c(1,1,2))
  event <- st_as_sf(event, coords = c("x","y"))

  event$sp_id <- sp_char_index(st_coordinates(event),1)

  # the bandwidth will be adapte arround the events
  # based on the density at the event, so we must first
  # calculate that density. we have here a case of
  # discontinuous NKDE with bw = 3 et kernel quartic

  n1 <- quartic_kernel(0,3)
  n2 <- quartic_kernel(0,3) + quartic_kernel(2,3) * 2
  n3 <- quartic_kernel(0,3) * 2 + quartic_kernel(2,3)

  hf0 <- c(n1,n2,n3)

  h0 <- 3
  gamma_val <- exp(sum(log(1/sqrt(hf0)))/3)
  abws <- h0 * (1/sqrt(hf0)) * (1/gamma_val)

  grid <- build_grid(c(1,1),list(all_lines))


  future::plan(future::multisession(workers=1))

  print(paste0("precalculated bws : ", paste0(abws, collapse = ';')))
  observed <- nkde.mc(
    events = event,
    w = event$w,
    samples = event,
    lines = all_lines,
    adaptive = TRUE,
    bw = 3,
    trim_bw = 5,
    method = "discontinuous",
    kernel_name = "quartic",
    max_depth = 8,
    tol = 0.1,
    digits = 2,
    sparse = TRUE,
    verbose = TRUE,
    grid_shape = c(3,3),
    check = FALSE
  )

  observed2 <- nkde.mc(
    events = event,
    w = event$w,
    samples = event,
    lines = all_lines,
    adaptive = TRUE,
    bw = 3,
    trim_bw = 5,
    method = "discontinuous",
    kernel_name = "quartic",
    max_depth = 8,
    tol = 0.1,
    digits = 2,
    sparse = TRUE,
    verbose = FALSE,
    grid_shape = c(3,3),
    check = FALSE
  )

  abws_2 <- observed$events$bw[match(event$sp_id, observed$events$spid)]
  abws_3 <- observed2$events$bw[match(event$sp_id, observed2$events$spid)]

  diff <- sum(round(abs(abws - abws_2),6))
  diff2 <- sum(round(abs(abws - abws_3),6))
  expect_equal(diff, 0, diff2)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE adaptive_bw FUNCTION WITH TWO POINTS AT THE SAME LOCATION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the adaptive_bw function WITH TWO POINTS AT THE SAME LOCATION", {

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
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,3),
                      w = c(1,1,1,1), # the last event has a weight of 2
                      id = 1:4)
  event <- st_as_sf(event, coords = c("x","y"))

  # the bandwidth will be adapted arround the events
  # based on the density at the event, so we must first
  # calculate that density. we have here a case of
  # discontinuous NKDE with bw = 3 et kernel quartic

  n1 <- quartic_kernel(0,3) + quartic_kernel(0,3)
  n2 <- quartic_kernel(0,3) + quartic_kernel(2,3)
  n3 <- quartic_kernel(0,3) + quartic_kernel(2,3)
  n4 <- n1

  hf0 <- c(n1,n2,n3)

  h0 <- 3
  gamma_val <- exp(sum(log(1/sqrt(hf0)))/3)
  abws <- h0 * (1/sqrt(hf0)) * (1/gamma_val)

  grid <- build_grid(c(1,1),list(all_lines))

  observed <- nkde(
    events = event,
    w = event$w,
    samples = event,
    lines = all_lines,
    adaptive = TRUE,
    bw = 3,
    trim_bw = 5,
    method = "discontinuous",
    kernel_name = "quartic",
    max_depth = 8,
    tol = 0.1,
    digits = 2,
    sparse = TRUE,
    verbose = TRUE,
    check = FALSE,
    grid_shape = c(3,3)
  )

  diff <- sum(round(abs(abws - observed$events$bw),6))
  expect_equal(diff, 0)
})

