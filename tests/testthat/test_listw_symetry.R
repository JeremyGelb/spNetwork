context("symetry and comparison between simple and multicore")

test_that("A listw object returned by the function network_listw must be symetric", {
  data(small_mtl_network)
  listw <- network_listw(small_mtl_network,small_mtl_network,
                                  method="centroid",
                                  maxdistance = 300,
                                  dist_func = "squared inverse",
                                  matrice_type = "B",
                                  grid_shape = c(1,1),
                                  verbose = "silent",
                                  mindist = 10,
                                  digits = 3
                                  )
  Value <- is.symmetric.nb(listw$neighbours)
  expect_true(Value)
})


test_that("A listw object returned by the function network_listw must be symetric", {
  data(small_mtl_network)
  listw <- network_listw(small_mtl_network,small_mtl_network,
                         method="ends",
                         maxdistance = 300,
                         dist_func = "squared inverse",
                         matrice_type = "B",
                         grid_shape = c(1,1),
                         verbose = "silent",
                         mindist = 10,
                         digits = 3
  )
  Value <- is.symmetric.nb(listw$neighbours)
  expect_true(Value)
})


test_that("A listw object returned by the function network_listw must be the same as a listw returned by network_listw.mc", {
  data(small_mtl_network)
  listw <- network_listw(small_mtl_network,small_mtl_network,
                         method="ends",
                         maxdistance = 300,
                         dist_func = "squared inverse",
                         matrice_type = "W",
                         grid_shape = c(1,1),
                         verbose = "silent",
                         mindist = 10,
                         digits = 3
  )
  future::plan(future::multiprocess(workers = 2))
  listwmc <- network_listw.mc(small_mtl_network,small_mtl_network,
                           method="ends",
                           maxdistance = 300,
                           dist_func = "squared inverse",
                           matrice_type = "W",
                           grid_shape = c(1,1),
                           verbose = "silent",
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
  test1 <- CompareNeighbours(listw$neighbours, listwmc$neighbours)
  test2 <- CompareWeights(listw$weights, listwmc$weights)
  tottest <- test1 & test2
  expect_true(tottest)
})

