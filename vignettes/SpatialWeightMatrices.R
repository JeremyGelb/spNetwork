## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
backup_option <- options()
base_wd <- getwd()

## ----include=FALSE------------------------------------------------------------
load(system.file("extdata", "results_vignette_wmat.rda",
                           package = "spNetwork", mustWork = TRUE))
library(spdep)

## ----message=FALSE, warning=FALSE, eval = FALSE, message=FALSE----------------
#  options("rgdal_show_exportToProj4_warnings"="none")
#  
#  library(spNetwork)
#  library(sf)
#  library(dplyr)
#  library(spdep)
#  
#  data(mtl_network)
#  data(bike_accidents)
#  
#  

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  lixels <- lixelize_lines(mtl_network,200,mindist = 50)
#  
#  # defining and oid for the lixels
#  lixels$oid <- 1:nrow(lixels)
#  
#  # snapping the points on the lines and counting
#  snapped_acc <- snapPointsToLines2(bike_accidents,lixels, idField ="oid")
#  counts <- table(snapped_acc$nearest_line_id)
#  counts_df <- data.frame("oid" = as.numeric(as.character(names(counts))),
#                          "count" = as.numeric(counts))
#  
#  # merging the results
#  lixels$nbAccident <- left_join(lixels,counts_df, by="oid")$count
#  nbAccident <- ifelse(is.na(lixels$nbAccident),0,lixels$nbAccident)

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  netlistw <- network_listw(lixels,mtl_network,
#                             method = "ends",
#                             mindist = 10,
#                             maxdistance = 300,
#                             dist_func = "squared inverse",
#                             line_weight = 'length',
#                             matrice_type = 'W',
#                             grid_shape = c(1,1),
#                             verbose=FALSE)

## ----message=FALSE, warning=FALSE---------------------------------------------
Test <- moran.test(nbAccident, netlistw, zero.policy = T)
print(round(Test$estimate,4))

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  my_conv_func <- function(x){
#    if (x>=300){
#      return(0)
#    }else{
#      return(1/x**3)
#    }
#  }
#  
#  netlistw2 <- network_listw(lixels,mtl_network,
#                               method = "ends",
#                               mindist = 10,
#                               maxdistance = 300,
#                               dist_func = my_conv_func,
#                               line_weight = 'length',
#                               matrice_type = 'W',
#                               grid_shape = c(1,1),
#                               verbose=FALSE)
#  

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  # setting the multiprocess plan
#  future::plan(future::multisession(workers=2))
#  
#  netlistw3 <- network_listw.mc(lixels,lixels,
#                               method = "ends",
#                               mindist = 10,
#                               maxdistance = 300,
#                               dist_func = my_conv_func,
#                               line_weight = 'length',
#                               matrice_type = 'W',
#                               grid_shape = c(2,2),
#                               verbose=FALSE)
#  
#  if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)

## ----include = FALSE----------------------------------------------------------
# reset all the user parameters
options(backup_option)
setwd(base_wd)
oldpar <- par(mfrow = c(1,2))
par(oldpar)

