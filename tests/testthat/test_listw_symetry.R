context("symetry and comparison between simple and multicore")

test_that("A listw object returned by the function network_listw must be symetric and identic as one returned by network_listw.mc", {
  #data(small_mtl_network)
  networkgpkg <- system.file("extdata", "networks.gpkg",
                             package = "spNetwork", mustWork = TRUE)
  small_mtl_network <- rgdal::readOGR(networkgpkg,layer="small_mtl_network", verbose = FALSE)
  listw <- network_listw(small_mtl_network,small_mtl_network,
                                  method="centroid",
                                  maxdistance = 300,
                                  dist_func = "squared inverse",
                                  matrice_type = "B",
                                  grid_shape = c(1,1),
                                  verbose = FALSE,
                                  mindist = 10,
                                  digits = 3
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
