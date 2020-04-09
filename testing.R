######## Testing before building

library(rgdal)
library(rgeos)
library(sp)
library(maptools)
library(tidyr)
library(dplyr)
library(sf)
library(raster)
library(igraph)


load('data/mtl_network.rda')
load('data/bike_accidents.rda')

source('R/graph_functions.R')
source('R/geometrical_functions.R')
source('R/nkde.R')


lixels <- nkde_grided(mtl_network,bike_accidents,snap_dist = 150, lixel_length = 100,
            kernel_range = 500,kerne='quartic',
            tol=0.1, digits = 3,grid_shape = c(5,5), show_progress = T)

future::plan(future::multiprocess(workers=4))
lixels2 <- nkde_grided.mc(mtl_network,bike_accidents,snap_dist = 150, lixel_length = 100,
                      kernel_range = 500,kerne='quartic',
                      tol=0.1, digits = 3,grid_shape = c(5,5), show_progress = T)

writeOGR(lixels,"mylixels.gpkg",layer="lixels",driver="GPKG")
