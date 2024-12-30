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
  observed <- correction_factor(
    study_area = polygons,
    events = event,
    lines = all_lines,
    method = "simple",
    bws = 2,
    kernel_name = "triangle",
    tol = 0.1,
    digits = 2,
    max_depth = 8,
    sparse=TRUE)

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR CORRECTION FACTOR WITH SIMPLE TNKDE AND ADAPTIVE BW ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


test_that("Testing the correction_factor function for a simple TNKDE and adaptive BWS", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 10, 0 0)",
    "LINESTRING (-10 0, 0 0)",
    "LINESTRING (0 -10, 0 0)",
    "LINESTRING (10 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  wkt_polygon <- "POLYGON ((-5 -4, 5 -4, 5 5, -5 5, -5 -4))"
  polydf <- data.frame(wkt = wkt_polygon,
                       id = 1)
  polygons <- st_as_sf(polydf, wkt = "wkt")

  time_bounds <- c(0,9)

  bws_net <- c(5,6)
  bws_time <- c(3,4)

  # definition of three events
  event <- data.frame(x=c(0,2,0),
                      y=c(2,0,-2),
                      id = c(1,2,3),
                      time = c(5,5,7))

  event <- st_as_sf(event, coords = c("x","y"))

  #tm_shape(all_lines) + tm_lines('black') + tm_shape(event) + tm_dots("red", size = 2) + tm_shape(polygons) + tm_borders("blue")

  # the first step is to evaluate the local densities
  # and then to obtain the local bandwidths
  loc_densities <- array(0, c(length(bws_net),length(bws_time),nrow(event)))
  loc_bws_net <- array(0, c(length(bws_net),length(bws_time),nrow(event)))
  loc_bws_time <- array(0, c(length(bws_net),length(bws_time),nrow(event)))


  for(i in 1:length(bws_net)){
    for(j in 1:length(bws_time)){
      bw_net <- bws_net[[i]]
      bw_time <- bws_time[[j]]

      # pour le point 1
      k1 <- (triangle_kernel(0,bw_net) *  triangle_kernel(0,bw_time)) / (bw_net * bw_time)
      k2 <- (triangle_kernel(4,bw_net) * triangle_kernel(0,bw_time)) / (bw_net * bw_time)
      k3 <- (triangle_kernel(4,bw_net) * triangle_kernel(2,bw_time)) / (bw_net * bw_time)
      dens1 <- k1 + k2 + k3

      # pour le point 2
      k1 <- (triangle_kernel(0,bw_net) *  triangle_kernel(0,bw_time)) / (bw_net * bw_time)
      k2 <- (triangle_kernel(4,bw_net) * triangle_kernel(0,bw_time)) / (bw_net * bw_time)
      k3 <- (triangle_kernel(4,bw_net) * triangle_kernel(2,bw_time)) / (bw_net * bw_time)
      dens2 <- k1 + k2 + k3

      # pour le point 3
      k1 <- (triangle_kernel(0,bw_net) *  triangle_kernel(0,bw_time)) / (bw_net * bw_time)
      k2 <- (triangle_kernel(4,bw_net) * triangle_kernel(2,bw_time)) / (bw_net * bw_time)
      k3 <- (triangle_kernel(4,bw_net) * triangle_kernel(2,bw_time)) / (bw_net * bw_time)
      dens3 <- k1 + k2 + k3

      # we can store de densities
      densities <- c(dens1, dens2, dens3)
      loc_densities[i,j,] <- densities

      # and now calculate the local bws for this pair of bws
      gamma_val <- exp(sum(log(1/sqrt(densities)))/3)
      loc_bws_net_vec <- bw_net * (1/sqrt(densities)) * (1/gamma_val)
      loc_bws_time_vec <- bw_time * (1/sqrt(densities)) * (1/gamma_val)
      loc_bws_net[i,j,] <- loc_bws_net_vec
      loc_bws_time[i,j,] <- loc_bws_time_vec

    }
  }

  # great, we will focus on two cases to make the evaluation easier
  # the last event is located at (0,-2, 7) (x,y,time)
  # let me check for the first bws in both space and time and the last
  e1_bw_net <- loc_bws_net[1,1,3]
  e1_bw_time <- loc_bws_time[1,1,3]

  # for this event, we know that a distance of two is separating it from the limit in space
  out_density_net <- cubintegrate(triangle_kernel,lower=2,upper=e1_bw_net,
                              bw=e1_bw_net, relTol = 1e-15)$integral

  # the time bounds are (0,9), so we also have some density outside
  out_density_time <- cubintegrate(triangle_kernel,lower=2, upper = e1_bw_time,
                                  bw=e1_bw_time, relTol = 1e-15)$integral

  # we can now calculate the in_density and the expected correction factor
  in_density <- (1-out_density_net) * (1-out_density_time)
  expected_val <- 1 / in_density
  expected_val2 <- (1/(1-out_density_net)) * (1/(1-out_density_time))

  event$goid <- 1:nrow(event)
  event$oid <- 1:nrow(event)

  # let us calculate the observed value
  observed <- bw_tnkde_corr_factor_arr(net_bws = loc_bws_net,
                                       time_bws = loc_bws_time,
                                       diggle_correction = TRUE,
                                       study_area = polygons,
                                       time_limits = time_bounds,
                                       events_loc = event,
                                       events = event,
                                       lines = all_lines,
                                       method = "simple",
                                       kernel_name = "triangle",
                                       tol = 0.0001,
                                       digits = 2,
                                       max_depth = 8)


  expect_equal(expected_val2, observed[1,1,3])

})

