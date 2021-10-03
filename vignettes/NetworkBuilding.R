## ----message=FALSE, warning=FALSE---------------------------------------------
# first load data and packages
library(sp)
library(maptools)
library(rgeos)
library(spNetwork)
library(raster)
library(tmap)
library(FNN)

networkgpkg <- system.file("extdata", "networks.gpkg",
                           package = "spNetwork", mustWork = TRUE)
eventsgpkg <- system.file("extdata", "events.gpkg",
                          package = "spNetwork", mustWork = TRUE)

mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network",verbose = FALSE)
bike_accidents <- rgdal::readOGR(eventsgpkg,layer="bike_accidents", verbose = FALSE)

# then plotting the data
plot(mtl_network)
plot(bike_accidents,add=T,col='red',pch = 19)

