context("testing function used to apply border correction in NKDE")
library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR CORRECTION FACTOR WITH SIMPLE NKDE ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the correction_factor function for a simple NKDE", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  wkt_polygon <- "POLYGON ((-5 -4, 5 -4, 5 5, -5 5, -5 -4))"
  polydf <- data.frame(wkt = wkt_polygon,
                        id = 1)
  polygons <- st_as_sf(polydf, wkt = "wkt")


  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3),
                      id = c(1,2,3))

  event <- st_as_sf(event, coords = c("x","y"))

  # evaluating the correction for the simple case
  # here, we only need to apply the correction to 1 event
  # in this case, the BW is 2, thus, a total length of 1 is
  # out of the area
  # let us calculate how much weight it represents
  out_density <- cubintegrate(triangle_kernel,lower=1,upper=2,
                              bw=2, relTol = 1e-15)$integral

  expected_val <- corr_factor <- 1/(1-out_density)

  expected_vals <- c(1,1,expected_val)
  # let us calculate the observed value
  observed <- correction_factor(polygons,event,all_lines,
                    "simple", 2, "triangle", 0.1, 2, 8, sparse=TRUE)

  expect_equal(sum(expected_vals - observed),0)


})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR CORRECTION FACTOR WITH DISCONTINUOUS NKDE ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the correction_factor function for a discontinuous NKDE", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)",
    "LINESTRING (-5 -5, 0 -5)",
    "LINESTRING (0 -5, 5 -5)"
    )

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  wkt_polygon <- "POLYGON ((-5 -4, 5 -4, 5 5, -5 5, -5 -4))"
  polydf <- data.frame(wkt = wkt_polygon,
                       id = 1)
  polygons <- st_as_sf(polydf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))

  # evaluating the correction for the discontinuous case
  # here, we only need to apply the correction to 1 event
  # in this case, the BW is 3
  # so we have a firts part outside on one line
  # and two other parts separated (correction factor)
  # we can calculate the density of the three parts
  dens1 <- cubintegrate(uniform_kernel,lower=1,upper=2,
                              bw=3, relTol = 1e-15)$integral

  dens2 <- cubintegrate(uniform_kernel,lower=2,upper=3,
                        bw=3, relTol = 1e-15)$integral * 0.5

  out_density <- dens1 + (2*dens2)

  expected_val <- corr_factor <- 1/(1-out_density)

  expected_vals <- c(1,1,expected_val)
  # let us calculate the observed value
  study_area = polygons
  events = event
  lines = all_lines
  method = "discontinuous"
  bws = c(3,3,3)
  kernel_name = "uniform"
  tol = 0.1
  digits = 2
  max_depth = 8
  sparse=TRUE
  observed <- correction_factor(study_area = polygons,
                                events = event,
                                lines = all_lines,
                                method = "discontinuous",
                                bws = c(3,3,3),
                                kernel_name = "uniform",
                                tol = 0.1,
                                digits = 2,
                                max_depth = 8,
                                sparse=TRUE)

  observed2 <- correction_factor(polygons,event,all_lines,
                                "discontinuous", c(3,3,3), "uniform", 0.1, 2, 8, sparse=FALSE)

  test1 <- round(sum(expected_vals - observed),10) == 0
  test2 <- round(sum(expected_vals - observed2),10) == 0

  expect_true(test1 & test2)

})


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# #### TEST FOR CORRECTION FACTOR WITH CONTINUOUS NKDE ####
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
test_that("Testing the correction_factor function for a continuous NKDE", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)",
    "LINESTRING (-5 -5, 0 -5)",
    "LINESTRING (0 -5, 5 -5)"
  )

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  wkt_polygon <- "POLYGON ((-5 -4, 5 -4, 5 5, -5 5, -5 -4))"
  polydf <- data.frame(wkt = wkt_polygon,
                       id = 1)

  polygons <- st_as_sf(polydf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,0),
                      y=c(3,0,-3),
                      id = c(1,2,3))
  event <- st_as_sf(event, coords = c("x","y"))

  # evaluating the correction for the discontinuous case
  # here, we only need to apply the correction to 1 event
  # in this case, the BW is 3
  # so we have a firts part outside on one line
  # and two other parts separated (correction factor)
  # we can calculate the density of the three parts

  n <- 3

  dens1 <- cubintegrate(cosine_kernel,lower=1,upper=2,
                        bw=3, relTol = 1e-15)$integral

  corr_dens <- cubintegrate(cosine_kernel,lower=2,upper=3,
                            bw=3, relTol = 1e-15)$integral * ((2.0-n)/n)

  dens2 <- cubintegrate(cosine_kernel,lower=2,upper=3,
                        bw=3, relTol = 1e-15)$integral * (2.0/n)

  out_density <- (dens1+corr_dens) + (2*dens2)

  expected_val <- corr_factor <- 1/(1-out_density)

  expected_vals <- c(1,1,expected_val)
  # let us calculate the observed value
  observed <- correction_factor(polygons,event,all_lines,
                                "continuous", c(3,3,3), "cosine", 0.1, 2, 8, sparse=TRUE)

  observed2 <- correction_factor(polygons,event,all_lines,
                                "continuous", c(3,3,3), "cosine", 0.1, 2, 8, sparse=FALSE)

  test1 <- round(sum(expected_vals - observed),10) == 0
  test2 <- round(sum(expected_vals - observed2),10) == 0

  expect_true(test1 & test2)

})

