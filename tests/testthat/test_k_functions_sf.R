context("testing functions for k and g function analysis")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE SIMPLE K AND G FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the simple k function", {

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
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))
  event <- st_as_sf(event, coords = c("x","y"))

  # calculating the observed values
  observed <- kfunctions(
                         lines = all_lines,
                         points = event,
                         start = 0,
                         end = 6,
                         step = 1,
                         width = 2,
                         nsim = 5,
                         conf_int = 0.05,
                         digits = 2,
                         tol = 0.1,
                         calc_g_func = TRUE,
                         resolution = NULL,
                         agg = NULL,
                         verbose = TRUE)

  observed2 <- kfunctions(
    lines = all_lines,
    points = event,
    start = 0,
    end = 6,
    step = 1,
    width = 2,
    nsim = 5,
    conf_int = 0.05,
    digits = 2,
    tol = 0.1,
    calc_g_func = FALSE,
    resolution = 1,
    agg = NULL,
    verbose = TRUE)

  # here we calculate the true value manually for the distance 3
  Lt <- sum(as.numeric(st_length(all_lines)))
  n <- nrow(event)
  t1 <- 1.0/((n-1)/Lt)

  # for the k value
  dist_matrix <- as.matrix(dist(st_coordinates(event), method = 'manhattan', upper = T))
  k <- mean((rowSums(dist_matrix <= 3) - 1)) * t1

  # for the g alue
  g <- mean((rowSums( (dist_matrix <= 4) & (dist_matrix >= 2) ))) * t1

  expected_vals <- c(k,g)
  test <- round(observed$values[4,c("obs_k","obs_g")],3) == round(expected_vals,3)
  test2 <- (any(!(round(observed$values$obs_k,3) == round(observed2$values$obs_k,3))) == FALSE)
  test <- c(test,test2)
  expect_equal(sum(test), 3)

})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE SIMPLE K AND G FUNCTION (multicore) ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the simple k function (multicore)", {

  skip_on_cran()

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
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))
  event <- st_as_sf(event, coords = c("x","y"))

  # calculating the observed values
  observed <- kfunctions.mc(
    lines = all_lines,
    points = event,
    start = 0,
    end = 6,
    step = 1,
    width = 2,
    nsim = 5,
    conf_int = 0.05,
    digits = 2,
    tol = 0.1,
    calc_g_func = TRUE,
    resolution = NULL,
    agg = NULL,
    verbose = TRUE)

  observed2 <- kfunctions(
    lines = all_lines,
    points = event,
    start = 0,
    end = 6,
    step = 1,
    width = 2,
    nsim = 5,
    conf_int = 0.05,
    digits = 2,
    tol = 0.1,
    calc_g_func = FALSE,
    resolution = 1,
    agg = NULL,
    verbose = TRUE)

  # here we calculate the true value manually for the distance 3
  Lt <- sum(as.numeric(st_length(all_lines)))
  n <- nrow(event)
  t1 <- 1.0/((n-1)/Lt)

  # for the k value
  dist_matrix <- as.matrix(dist(st_coordinates(event), method = 'manhattan', upper = T))
  k <- mean((rowSums(dist_matrix <= 3) - 1)) * t1

  # for the g alue
  g <- mean((rowSums( (dist_matrix <= 4) & (dist_matrix >= 2) ))) * t1

  expected_vals <- c(k,g)
  test <- round(observed$values[4,c("obs_k","obs_g")],3) == round(expected_vals,3)
  test2 <- (any(!(round(observed$values$obs_k,3) == round(observed2$values$obs_k,3))) == FALSE)
  test <- c(test,test2)
  expect_equal(sum(test), 3)

})



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE cross K AND G FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the cross k function", {

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
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))

  event <- st_as_sf(event, coords = c("x","y"))

  As <- event[c(1,2),]
  Bs <- event[c(3,4),]

  # calculating the observed values

  observed <- cross_kfunctions(
    lines = all_lines,
    pointsA = As,
    pointsB = Bs,
    start = 0,
    end = 6,
    step = 1,
    width = 2,
    nsim = 5,
    conf_int = 0.05,
    digits = 2,
    tol = 0.1,
    calc_g_func = TRUE,
    resolution = NULL,
    agg = NULL,
    verbose = TRUE)

  observed2 <- cross_kfunctions(
    lines = all_lines,
    pointsA = As,
    pointsB = Bs,
    start = 0,
    end = 6,
    step = 1,
    width = 2,
    nsim = 5,
    conf_int = 0.05,
    digits = 2,
    tol = 0.1,
    calc_g_func = FALSE,
    resolution = NULL,
    agg = NULL,
    verbose = TRUE)

  # here we calculate the true value manually for the distance 3
  Lt <- sum(as.numeric(st_length(all_lines)))
  n <- nrow(event)
  t1 <- 1.0/((nrow(As))/Lt)

  # for the k value
  dist_matrix <- rbind(c(4,2), c(2,4))

  k <- mean((rowSums(dist_matrix <= 3))) * t1

  # for the g alue
  g <- mean((rowSums( (dist_matrix <= 4) & (dist_matrix >= 2) ))) * t1

  expected_vals <- c(k,g)
  test <- round(observed$values[4,c("obs_k","obs_g")],3) == round(expected_vals,3)
  test2 <- (any(!(observed$values$obs_k == observed2$values$obs_k)) == FALSE)
  test <- c(test,test2)
  expect_equal(sum(test), 3)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE cross K AND G FUNCTION (MULTICORE) ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the cross k function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  skip_on_cran()

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))

  event <- st_as_sf(event, coords = c("x","y"))

  As <- event[c(1,2),]
  Bs <- event[c(3,4),]

  # calculating the observed values

  observed <- cross_kfunctions.mc(
    lines = all_lines,
    pointsA = As,
    pointsB = Bs,
    start = 0,
    end = 6,
    step = 1,
    width = 2,
    nsim = 5,
    conf_int = 0.05,
    digits = 2,
    tol = 0.1,
    calc_g_func = TRUE,
    resolution = NULL,
    agg = NULL,
    verbose = TRUE,
    grid_shape = c(1,1)
    )

  observed2 <- cross_kfunctions(
    lines = all_lines,
    pointsA = As,
    pointsB = Bs,
    start = 0,
    end = 6,
    step = 1,
    width = 2,
    nsim = 5,
    conf_int = 0.05,
    digits = 2,
    tol = 0.1,
    calc_g_func = TRUE,
    resolution = NULL,
    agg = NULL,
    verbose = TRUE)

  # here we calculate the true value manually for the distance 3
  Lt <- sum(as.numeric(st_length(all_lines)))
  n <- nrow(event)
  t1 <- 1.0/((nrow(As))/Lt)

  # for the k value
  dist_matrix <- rbind(c(4,2), c(2,4))

  k <- mean((rowSums(dist_matrix <= 3))) * t1

  # for the g alue
  g <- mean((rowSums( (dist_matrix <= 4) & (dist_matrix >= 2) ))) * t1

  expected_vals <- c(k,g)
  test <- round(observed$values[4,c("obs_k","obs_g")],3) == round(expected_vals,3)
  test2 <- (any(!(observed$values$obs_k == observed2$values$obs_k)) == FALSE)
  test <- c(test,test2)
  expect_equal(sum(test), 3)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE SIMPLE space-time K AND G FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the simple k function in space-time", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- sf::st_as_sf(linesdf, wkt = "wkt")

  # definition of four events
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4),
                      time = c(1,1,3,2))

  event <- st_as_sf(event, coords = c("x","y"))

  # calculating the observed values
  future::plan(future::multisession(workers=1))
  observed <- k_nt_functions(
                        lines = all_lines,
                         points = event,
                         points_time = event$time,
                         start_net = 0,
                         end_net = 6,
                         step_net = 1,
                         width_net = 2,
                         start_time = 0,
                         end_time = 3,
                         step_time = 1,
                         width_time = 2,
                         nsim = 5,
                         conf_int = 0.05,
                         digits = 2,
                         tol = 0.1,
                         resolution = NULL,
                         agg = NULL,
                         verbose = TRUE,
                         calc_g_func = TRUE
                         )

  # after checking on a paper with a pen, the observed k and g values at  network
  # distance 3 and time distance 1 must be :
  Lt <- sum(st_length(all_lines))
  Tt <- max(event$time) - min(event$time)
  n <- nrow(event)
  p <- (n-1)/(Lt * Tt);

  k_counts <- c(0,1,1,2)
  g_counts <- c(2,2,3,3)

  exp_k <- (1/p) * mean(k_counts)
  exp_g <- (1/p) * mean(g_counts)
  expected_vals <- c(exp_k, exp_g)

  obtained <- c(observed$obs_k[4,2], observed$obs_g[4,2])

  test <- sum(round(obtained - expected_vals,5)) == 0
  expect_true(test)

})



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE SIMPLE space-time K AND G FUNCTION (multicoce)  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the simple k function in space-time multicore", {

  skip_on_cran()

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of four events
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4),
                      time = c(1,1,3,2))

  event <- st_as_sf(event, coords = c("x","y"))

  # calculating the observed values
  observed <- k_nt_functions(lines = all_lines,
                             points = event,
                             points_time = event$time,
                             start_net = 0,
                             end_net = 6,
                             step_net = 0.5,
                             width_net = 2,
                             start_time = 0,
                             end_time = 6,
                             step_time = 0.5,
                             width_time = 2,
                             nsim = 5,
                             conf_int = 0.05,
                             digits = 2,
                             tol = 0.1,
                             resolution = NULL,
                             agg = NULL,
                             verbose = TRUE)


  observed.mc <- k_nt_functions.mc(lines = all_lines,
                             points = event,
                             points_time = event$time,
                             start_net = 0,
                             end_net = 6,
                             step_net = 0.5,
                             width_net = 2,
                             start_time = 0,
                             end_time = 6,
                             step_time = 0.5,
                             width_time = 2,
                             nsim = 5,
                             conf_int = 0.05,
                             digits = 2,
                             tol = 0.1,
                             resolution = NULL,
                             agg = NULL,
                             verbose = TRUE,
                             grid_shape = c(2,2))


  t1 <- sum(round(observed$obs_g,3) - round(observed.mc$obs_g,3))
  t2 <- sum(round(observed$obs_k,3) - round(observed.mc$obs_k,3))

  test <- (t1 + t2) == 0

  expect_true(test)

})


