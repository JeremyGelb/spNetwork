context("testing functions for adptative bandwidth in NKDE")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE adaptive_bw FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the adaptive_bw function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # definition of three events
  event <- data.frame(x=c(0,3,1),
                      y=c(3,0,0))
  sp::coordinates(event) <- cbind(event$x,event$y)

  # the bandwidth will be adapte arround the events
  # based on the density at the event, so we must first
  # calculate that density. we have here a case of
  # discontinuous NKDE with bw = 3 et kernel quartic

  n1 <- quartic_kernel(0,3)
  n2 <- quartic_kernel(0,3) + quartic_kernel(2,3)
  n3 <- quartic_kernel(0,3) + quartic_kernel(2,3)

  hf0 <- c(n1,n2,n3)

  h0 <- 3
  gamma_val <- exp(sum(log(1/sqrt(hf0)))/3)
  abws <- h0 * (1/sqrt(hf0)) * (1/gamma_val)

  grid <- build_grid(c(1,1),list(all_lines))

  observed <- nkde(
    events = event,
    w = rep(1,nrow(event)),
    samples = event,
    lines = all_lines,
    adaptive = TRUE,
    bw = 3, trim_bw = 5, method = "discontinuous",
    kernel_name = "quartic", max_depth = 8,
    tol = 0.1, digits = 2, sparse = TRUE, verbose = TRUE,
    check = FALSE
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

  geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
    txt <- as.character(linesdf[i,]$wkt)
    geom <- rgeos::readWKT(txt,id=i)
    return(geom)
  }))

  all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

  # definition of three events
  event <- data.frame(x=c(0,3,1),
                      y=c(3,0,0))
  sp::coordinates(event) <- cbind(event$x,event$y)

  # the bandwidth will be adapte arround the events
  # based on the density at the event, so we must first
  # calculate that density. we have here a case of
  # discontinuous NKDE with bw = 3 et kernel quartic

  n1 <- quartic_kernel(0,3)
  n2 <- quartic_kernel(0,3) + quartic_kernel(2,3)
  n3 <- quartic_kernel(0,3) + quartic_kernel(2,3)

  hf0 <- c(n1,n2,n3)

  h0 <- 3
  gamma_val <- exp(sum(log(1/sqrt(hf0)))/3)
  abws <- h0 * (1/sqrt(hf0)) * (1/gamma_val)

  grid <- build_grid(c(1,1),list(all_lines))


  future::plan(future::multisession(workers=1))
  observed <- nkde.mc(
    events = event,
    w = rep(1,nrow(event)),
    samples = event,
    lines = all_lines,
    adaptive = TRUE,
    bw = 3, trim_bw = 5, method = "discontinuous",
    kernel_name = "quartic", max_depth = 8,
    tol = 0.1, digits = 2, sparse = TRUE, verbose = TRUE,
    check = FALSE
  )

  diff <- sum(round(abs(abws - observed$events$bw),6))
  expect_equal(diff, 0)
})

