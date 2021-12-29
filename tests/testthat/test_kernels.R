context("testing the kernel functions")
library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the kernels with a SPARSE matrix ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the simple kernel with a simple case", {
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
  event <- data.frame(x=c(0.1),
                      y=c(3.1),
                      oid = 1)
  event <- st_as_sf(event, coords = c("x","y"))

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                      y=c(0.1),
                      oid = 1)
  sp_point <- st_as_sf(sp_point, coords = c("x","y"))

  # real distance is 6, and let us say that the bw is 10
  real_value <- (1/10) * quartic_kernel(6,10)

  #let us calculate the value with our function
  obs_value <- nkde(
       all_lines,
       events = event,
       w = c(1),
       samples = sp_point,
       check = FALSE,
       kernel_name = "quartic",
       bw = 10,
       adaptive = F,
       method = "simple",
       div = "bw",
       agg = 0.01,
       verbose =TRUE,
       tol = 0.0001
       )
  expect_equal(obs_value, real_value)
})


test_that("Testing the discontinuous kernel with a simple case", {
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
  event <- data.frame(x=c(0.1),
                      y=c(3.1),
                      id = 1)
  event <- st_as_sf(event,
                    coords = c("x","y"))

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                         y=c(0.1),
                         id = 1)

  sp_point <- st_as_sf(sp_point,
                    coords = c("x","y"))


  # real distance is 6, and let us say that the bw is 10
  real_value <- (1/10) * quartic_kernel(6,10) * 1/3

  #let us calculate the value with our function
  obs_value <- nkde(all_lines,events = event, w = c(1),
                    samples = sp_point, check = F,
                    kernel_name = "quartic",
                    bw = 10, adaptive = F, method = "discontinuous", div = "bw",
                    agg = 0.01, verbose = F,tol = 0.0001
  )
  expect_equal(obs_value, real_value)
})


test_that("Testing the continuous kernel with a simple case", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 10.1, 0.1 0.1)",
    "LINESTRING (-10.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -10.1, 0.1 0.1)",
    "LINESTRING (10.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf,
                        wkt = "wkt")

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1),
                      id = 1)

  event <- st_as_sf(event,
                       coords = c("x","y"))

  # definition of one sampling point
  sp_point <- data.frame(x=c(0.1),
                         y=c(1.1),
                         id = 1)

  sp_point <- st_as_sf(sp_point,
                    coords = c("x","y"))


  bw <- 5
  # real distance is 2, and let us say that the bw is 5
  real_value <- quartic_kernel(2,bw) - ((1/2) * quartic_kernel(4,bw))

  #let us calculate the value with our function
  obs_value <- nkde(all_lines,events = event, w = c(1),
                    samples = sp_point,check = F,
                    kernel_name = "quartic",
                    bw = bw, adaptive = F, method = "continuous", div = "none",
                    agg = 0.01, verbose = F,tol = 0.01,digits = 3,sparse = T,
  )
  expect_equal(obs_value, real_value)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### comparing results between sparse and integer matrix ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("comparing the simple kernel obtained for a sparse and integer matrix", {

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
  event <- data.frame(x=c(0.1),
                      y=c(3.1),
                      id = 1)

  event <- st_as_sf(event,
                       coords = c("x","y"))

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                         y=c(0.1),
                         id = 1)

  sp_point <- st_as_sf(sp_point,
                    coords = c("x","y"))


  #let us calculate the value with our function
  sparse_value <- nkde(all_lines,events = event, w = c(1),
                    samples = sp_point, check = F,
                    kernel_name = "quartic",
                    bw = 10, adaptive = F, method = "simple", div = "bw",
                    sparse = TRUE,
                    agg = 0.01, verbose = F,tol = 0.0001
  )

  integer_value <- nkde(all_lines,events = event, w = c(1),
                       samples = sp_point, check = F,
                       kernel_name = "quartic",
                       bw = 10, adaptive = F, method = "simple", div = "bw",
                       sparse = FALSE,
                       agg = 0.01, verbose = F,tol = 0.0001
  )

  expect_equal(sparse_value, integer_value)
})


test_that("comparing the discontinuous kernel obtained for a sparse and integer matrix", {

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
  event <- data.frame(x=c(0.1),
                      y=c(3.1),
                      id = 1)

  event <- st_as_sf(event,
                       coords = c("x","y"))

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                         y=c(0.1),
                         id = 1)

  sp_point <- st_as_sf(sp_point,
                    coords = c("x","y"))


  #let us calculate the value with our function
  sparse_value <- nkde(all_lines,events = event, w = c(1),
                       samples = sp_point, check = F,
                       kernel_name = "quartic",
                       bw = 10, adaptive = F, method = "discontinuous", div = "bw",
                       sparse = TRUE,
                       agg = 0.01, verbose = F,tol = 0.0001
  )

  integer_value <- nkde(all_lines,events = event, w = c(1),
                        samples = sp_point, check = F,
                        kernel_name = "quartic",
                        bw = 10, adaptive = F, method = "discontinuous", div = "bw",
                        sparse = FALSE,
                        agg = 0.01, verbose = F,tol = 0.0001
  )

  expect_equal(sparse_value, integer_value)
})


test_that("comparing the continuous kernel obtained for a sparse and integer matrix", {

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
  event <- data.frame(x=c(0.1),
                      y=c(3.1),
                      id = 1)
  event <- st_as_sf(event,
                       coords = c("x","y"))

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                         y=c(0.1),
                         id = 1)
  sp_point <- st_as_sf(sp_point,
                    coords = c("x","y"))


  #let us calculate the value with our function
  sparse_value <- nkde(all_lines,events = event, w = c(1),
                       samples = sp_point, check = F,
                       kernel_name = "quartic",
                       bw = 7, adaptive = F, method = "continuous", div = "none",
                       sparse = TRUE,
                       agg = 0.01, verbose = F,tol = 0.0001
  )

  integer_value <- nkde(all_lines,events = event, w = c(1),
                        samples = sp_point, check = F,
                        kernel_name = "quartic",
                        bw = 7, adaptive = F, method = "continuous", div = "none",
                        sparse = FALSE,
                        agg = 0.01, verbose = F,tol = 0.0001
  )

  expect_equal(sparse_value, integer_value)
})

