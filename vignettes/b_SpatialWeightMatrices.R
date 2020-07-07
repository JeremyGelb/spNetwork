## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")
library(spNetwork)
library(maptools)
library(dplyr)
library(spdep)

data("mtl_network")
data("bike_accidents")


# set this function to silent in the vignette
lixels <- lixelize_lines(mtl_network,200,mindist = 50)

## ----message=FALSE, warning=FALSE---------------------------------------------
# defining and oid for the lines
lixels$oid <- 1:nrow(lixels)


snapped_acc <- snapPointsToLines(bike_accidents,lixels, idField ="oid")
counts <- table(snapped_acc$nearest_line_id)
counts_df <- data.frame("oid" = as.numeric(as.character(names(counts))),
                        "count" = as.numeric(counts))

lixels$nbAccident <- left_join(lixels@data,counts_df, by="oid")$count
lixels$nbAccident <- ifelse(is.na(lixels$nbAccident),0,lixels$nbAccident)

## ----message=FALSE, warning=FALSE---------------------------------------------
netlistw <- network_listw(lixels,mtl_network,
                           method = "ends",
                           mindist = 10,
                           maxdistance = 300,
                           dist_func = "squared inverse",
                           line_weight = 'length',
                           matrice_type = 'W',
                           grid_shape = c(1,1),
                           verbose=FALSE)

Test <- moran.test(lixels$nbAccident, netlistw, zero.policy = T)
print(round(Test$estimate,4))

## ----message=FALSE, warning=FALSE---------------------------------------------
my_conv_func <- function(x){
  if (x>=300){
    return(0)
  }else{
    return(1/x**3)
  }
}


netlistw2 <- network_listw(lixels,mtl_network,
                             method = "ends",
                             mindist = 10,
                             maxdistance = 300,
                             dist_func = my_conv_func,
                             line_weight = 'length',
                             matrice_type = 'W',
                             grid_shape = c(1,1),
                             verbose=FALSE)

Test2 <- moran.test(lixels$nbAccident, netlistw2, zero.policy = T)
print(round(Test2$estimate,4))

## ----message=FALSE, warning=FALSE---------------------------------------------

# setting the multiprocess plan
# future::plan(future::multiprocess(workers=4))
# 
# netlistw3 <- network_listw.mc(lixels,lixels,
#                              method = "ends",
#                              mindist = 10,
#                              maxdistance = 300,
#                              dist_func = my_conv_func,
#                              line_weight = 'length',
#                              matrice_type = 'W',
#                              grid_shape = c(3,3),
#                              verbose=FALSE)
# 
# if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)

