context("symetry and comparison between simple and multicore")
library(sf)

test_that("A listw object returned by the function network_listw must be symetric and identic as one returned by network_listw.mc", {

  skip_on_cran()

  data(small_mtl_network)

  listw <- network_listw(origins = small_mtl_network,
                         lines = small_mtl_network,
                         method="centroid",
                         maxdistance = 300,
                         dist_func = "squared inverse",
                         matrice_type = "B",
                         grid_shape = c(1,1),
                         verbose = FALSE,
                         mindist = 10,
                         digits = 3,
                         direction = NULL,
                         point_dist = NULL,
                         snap_dist = Inf,
                         line_weight = 'length',
                         tol = 0.1
                         )
  # test if symetric
  test1 <- spdep::is.symmetric.nb(listw$neighbours)

  # now test if identical
  future::plan(future::multisession(workers = 1))
  listwmc <- network_listw.mc(small_mtl_network,small_mtl_network,
                              method="centroid",
                              maxdistance = 300,
                              dist_func = "squared inverse",
                              matrice_type = "B",
                              grid_shape = c(1,1),
                              verbose = FALSE,
                              mindist = 10,
                              digits = 3
  )
  CompareNeighbours <- function(nb1, nb2){
    diff <- sapply(1:length(nb1),function(i){
      return(any((nb1[[i]]==nb2[[i]])==F))
    })
    Value <- sum(diff) == 0
    return(Value)
  }

  CompareWeights <- function(nb1,nb2){
    diff <- sapply(1:length(nb1),function(i){
      return(sum(abs(nb1[[i]]-nb2[[i]])))
    })
    Value <- round(sum(diff),8) == 0
    return(Value)
  }

  test2 <- CompareNeighbours(listw$neighbours, listwmc$neighbours)
  test3 <- CompareWeights(listw$weights, listwmc$weights)
  tottest <- test1 & test2 & test3

  expect_true(tottest)
})



test_that("The distances returned by network_listw are valid ", {


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

  listw <- network_listw(origins = event,
                         lines = all_lines,
                         method="centroid",
                         maxdistance = 300,
                         dist_func = "identity",
                         matrice_type = "I",
                         grid_shape = c(1,1),
                         verbose = FALSE,
                         mindist = 0.1,
                         digits = 3,
                         direction = NULL,
                         point_dist = NULL,
                         snap_dist = Inf,
                         line_weight = 'length',
                         tol = 0.1
  )

  expected <- rbind(
    c(6,4,2),
    c(6,2,4),
    c(4,2,2),
    c(2,4,2)
  )

  obtained <- do.call(rbind, listw$weights)

  diff <- sum(round(expected - obtained,2))

  expect_equal(diff, 0)

})




