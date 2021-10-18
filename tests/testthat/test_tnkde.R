context("testing the kernel functions")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the kernels with a SPARSE matrix ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## -- TEST FOR THE SIMPLE TNKDE -- ##

test_that("Testing the simple tnkde with a simple case", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 5.1, 0.1 0.1)",
    "LINESTRING (-5.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -5.1, 0.1 0.1)",
    "LINESTRING (5.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1))
  sp::coordinates(event) <- cbind(event$x,event$y)
  event$time <- 5

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                      y=c(0.1))
  sp::coordinates(sp_point) <- cbind(sp_point$x,sp_point$y)
  sample_times <- c(1,2,3,4,5)


  # real distance is 6, and let us say that the bw_net is 10 and bw_times is 5
  net_density <- (1/10) * quartic_kernel(6,10)
  time_densities <- (1/5) * quartic_kernel(abs(sample_times-event$time),5)

  result_mat <-sapply(time_densities, function(x){
    x*net_density
  })

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
       bw_net  = 10,
       bw_time = 5,
       adaptive = F,
       method = "simple",
       div = "bw",
       agg = 0.01,
       verbose = F,
       tol = 0.0001,
       digits = 2
       )
  comp <- sum(round(abs(obs_value - result_mat),9))
  expect_equal(comp, 0)
})


## -- TEST FOR THE DISCONTINUOUS TNKDE -- ##

test_that("Testing the discontinuous tnkde with a simple case", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 5.1, 0.1 0.1)",
    "LINESTRING (-5.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -5.1, 0.1 0.1)",
    "LINESTRING (5.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1))
  sp::coordinates(event) <- cbind(event$x,event$y)
  event$time <- 5

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                         y=c(0.1))
  sp::coordinates(sp_point) <- cbind(sp_point$x,sp_point$y)
  sample_times <- c(1,2,3,4,5)


  # real distance is 6, and let us say that the bw_net is 10 and bw_times is 5
  net_density <- (1/10) * quartic_kernel(6,10) * 1/3
  time_densities <- (1/5) * quartic_kernel(abs(sample_times-event$time),5)

  result_mat <-sapply(time_densities, function(x){
    x*net_density
  })

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
    bw_net  = 10,
    bw_time = 5,
    adaptive = F,
    method = "discontinuous",
    div = "bw",
    agg = 0.01,
    verbose = F,
    tol = 0.0001,
    digits = 2
  )
  comp <- sum(round(abs(obs_value - result_mat),9))
  expect_equal(comp, 0)
})



## -- TEST FOR THE CONTINUOUS TNKDE -- ##

test_that("Testing the continuous tnkde with a simple case", {

  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 10.1, 0.1 0.1)",
    "LINESTRING (-10.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -10.1, 0.1 0.1)",
    "LINESTRING (10.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1))
  sp::coordinates(event) <- cbind(event$x,event$y)
  event$time <- 5

  # definition of one sampling point
  sp_point <- data.frame(x=c(0.1),
                         y=c(1.1))
  sp::coordinates(sp_point) <- cbind(sp_point$x,sp_point$y)

  # real distance is 2, and let us say that the network bw is 5 and bw_times is 5
  net_density <- (1/5) * (quartic_kernel(2,5) - ((1/2) * quartic_kernel(4,5)))
  time_densities <- (1/5) * quartic_kernel(abs(sample_times-event$time),5)

  result_mat <-sapply(time_densities, function(x){
    x*net_density
  })

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
    bw_net  = 5,
    bw_time = 5,
    adaptive = F,
    method = "continuous",
    div = "bw",
    agg = 0.01,
    verbose = F,
    tol = 0.0001,
    digits = 2
  )
  comp <- sum(round(abs(obs_value - result_mat),9))
  expect_equal(comp, 0)

})



