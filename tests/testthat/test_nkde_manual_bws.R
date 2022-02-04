context("testing the kernel functions with manually changed bandwidths")
library(sf)

#### ONLY FOR SPATIAL CASE ####

test_that("Testing the discontinuous kernel with a simple case with manually selected bandwidths", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 5.1, 0.1 0.1)",
    "LINESTRING (-5.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -5.1, 0.1 0.1)",
    "LINESTRING (5.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf,
                        wkt = "wkt")

  # definition of one event
  event <- data.frame(x=c(0.1, -3.1),
                      y=c(3.1, 0.1),
                      id = 1)
  event <- st_as_sf(event,
                    coords = c("x","y"))

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                         y=c(0.1),
                         id = 1)

  sp_point <- st_as_sf(sp_point,
                       coords = c("x","y"))


  # real distance is 6, and let us say that the bw is 10 for the first point and 12 for the second
  real_value <- ((1/10) * quartic_kernel(6,10) * 1/3) + ((1/12) * quartic_kernel(6.2,12) * 1/3)

  #let us calculate the value with our function
  obs_value <- nkde(lines = all_lines,
                    events = event,
                    w = c(1,1),
                    samples = sp_point,
                    check = F,
                    kernel_name = "quartic",
                    trim_bw  = 20,
                    bw = c(10,12),
                    adaptive = FALSE,
                    method = "discontinuous",
                    div = "bw",
                    agg = 0.01,
                    verbose = F,
                    tol = 0.0001
  )
  expect_equal(obs_value, real_value)
})


#### FOR THE SPATIO-TEMPORAL CASE ###
test_that("Testing the discontinuous tnkde with manually selected bandwidths", {

  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 5.1, 0.1 0.1)",
    "LINESTRING (-5.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -5.1, 0.1 0.1)",
    "LINESTRING (5.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of one event
  event <- data.frame(x=c(0.1, -3.1),
                      y=c(3.1, 0.1))
  event <- st_as_sf(event, coords = c("x","y"))
  event$time <- c(5, 5)

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                         y=c(0.1))
  sp_point <- st_as_sf(sp_point, coords = c("x","y"))
  sample_times <- c(1,2,3,4,5)


  # real distance is 6, and let us say that the bw_net is 10 and bw_times is 5 for event 1
  # and 12 and 7 for event 2
  net_density1 <- (1/10) * quartic_kernel(6,10) * 1/3
  time_densities1 <- (1/5) * quartic_kernel(abs(sample_times-event$time[[1]]),5)

  net_density2 <- (1/12) * quartic_kernel(6.2,12) * 1/3
  time_densities2 <- (1/7) * quartic_kernel(abs(sample_times-event$time[[2]]),7)

  result_mat1 <- sapply(time_densities1, function(x){
    x*net_density1
  })

  result_mat2 <- sapply(time_densities2, function(x){
    x*net_density2
  })

  result_mat <- result_mat2 + result_mat1

  #let us calculate the value with our function
  obs_value <- tnkde(
    lines = all_lines,
    events = event,
    w = c(1),
    time_field = "time",
    samples_loc = sp_point,
    samples_time = sample_times,
    check = F,
    kernel_name = "quartic",
    bw_net  = c(10,12),
    bw_time = c(5,7),
    adaptive = F,
    method = "discontinuous",
    div = "bw",
    agg = 0.01,
    verbose = F,
    tol = 0.0001,
    digits = 2,
    grid_shape = c(1,1),
    sparse = TRUE,
    max_depth = 10
  )

  test <- sum(round(abs(obs_value - result_mat),18)) == 0

  expect_true(test)
})
