## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------

#first load data and packages
#devtools::load_all("..")
library(spNetwork)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
data(mtl_network)
data(bike_accidents)

#then plotting the data
plot(mtl_network)
plot(bike_accidents,add=T,col='red')


# #then applying the NKDE
# lixels <- nkde_grided(mtl_network,bike_accidents,
#             snap_dist = 150,
#             lx_length=200,mindist=50,
#             kernel_range = 800, kernel='quartic',
#             weights=NULL, grid_shape = c(5,5),
#             show_progress = F)
# 
# #the generated lixels have the desired length
# hist(gLength(lixels,byid=T),breaks=50)

