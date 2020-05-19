## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------

#first load data and packages
library(spNetwork)
library(sp)
library(maptools)
library(rgeos)
data(mtl_network)
data(bike_accidents)

#then plotting the data
plot(mtl_network)
plot(bike_accidents,add=T,col='red')


#then applying the NKDE
lixels <- nkde(mtl_network,bike_accidents,
            snap_dist = 150,
            lx_length=200,mindist=50,
            kernel_range = 300, kernel='quartic',
            weights=NULL, grid_shape = c(5,5),
            verbose = "silent")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(classInt)

## a little function to help the mapping
get_colors <-function(x,colors,breaks){
  mycolors <- sapply(x,function(xi){
    for(i in seq(1:length(colors))){
      if (xi<breaks$brks[[i+1]]){
        return(colors[[i]])
      }
    }
    return(colors[[length(colors)]])
  })
  return(mycolors)
}

#setting an easier scale for mapping
lixels$mapdensity <- round(lixels$density*1000,4)

#using a discretization method
breaks <- classIntervals(lixels$mapdensity, n = 7, style = "fisher", intervalClosure = "right")
PaletteCouleur <- brewer.pal(n = 7, name = "Spectral")
PaletteCouleur <- rev(PaletteCouleur)
lixels$class <- get_colors(lixels$mapdensity,PaletteCouleur, breaks)

#and finally map with ggplot
labels <- names(print(breaks))
MapData <- fortify(lixels,id="lxid")
MapData$id <- as.numeric(MapData$id)
MapData <- left_join(MapData, lixels@data, by=c("id"="lxid"))
ggplot(MapData) + 
  geom_path(aes(x=long,y=lat,group=group,color=class))+
  scale_color_manual("density",
    breaks = PaletteCouleur, values = PaletteCouleur, 
    label = labels)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_fixed()+
  ggtitle("bike accident density by kilometers in 2016",
          subtitle = "within a radius of 300 meters")
  


