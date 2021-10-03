## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
backup_option <- options()
base_wd <- getwd()

## ----message=FALSE, warning=FALSE, fig.height = 3, fig.width  = 3-------------
library(spNetwork)
library(rgeos)
library(sp)
library(maptools)

# start with de definition of some lines
wkt_lines <- c(
  "LINESTRING (0.0 0.0, 5.0 0.0)",
  "LINESTRING (0.0 -5.0, 5.0 -5.0)",
  "LINESTRING (5.0 0.0, 5.0 5.0)",
  "LINESTRING (5.0 -5.0, 5.0 -10.0)",
  "LINESTRING (5.0 0.0, 5.0 -5.0)",
  "LINESTRING (5.0 0.0, 10.0 0.0)",
  "LINESTRING (5.0 -5.0, 10.0 -5.0)",
  "LINESTRING (10.0 0, 10.0 -5.0)",
  "LINESTRING (10.0 -10.0, 10.0 -5.0)",
  "LINESTRING (15.0 -5.0, 10.0 -5.0)",
  "LINESTRING (10.0 0.0, 15.0 0.0)",
  "LINESTRING (10.0 0.0, 10.0 5.0)")

linesdf <- data.frame(wkt = wkt_lines,
                      id = paste("l",1:length(wkt_lines),sep=""))

geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
      txt <- as.character(linesdf[i,]$wkt)
      geom <- rgeos::readWKT(txt,id=i)
      return(geom)
    }))

all_lines <- SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

# and the definition of one event
event <- data.frame(x=c(5),
                    y=c(-2.5))
coordinates(event) <- cbind(event$x,event$y)

# and map the situation
par(mar=c(0.1,0.1,0.1,0.1))
sp::plot(all_lines)
sp::plot(event,add=T,col="red",pch = 19)

## ----message=FALSE, warning=FALSE---------------------------------------------
samples_pts <- lines_points_along(all_lines,0.01)

simple_kernel <- nkde(all_lines, event, w = 1,
                  samples = samples_pts, kernel_name = "quartic", 
                  bw = 10, method = "simple", div = "bw", 
                  digits = 3, tol = 0.001, grid_shape = c(1,1),
                  check = FALSE,
                  verbose = FALSE)

## ----message=FALSE, warning=FALSE, dev='png', dpi=300, out.width = "50%"------
library(plot3D)

zfactor <- 1000

par(mar=c(0.1,0.1,0.1,0.1))

# step1 : plot the lines
lines_coords <- coordinates(all_lines)
x0 <- c()
y0 <- c()
x1 <- c()
y1 <- c()
for(i in 1:length(lines_coords)){
  ci <- lines_coords[[i]][[1]]
  x0 <- c(x0,ci[1,1])
  y0 <- c(y0,ci[1,2])
  x1 <- c(x1,ci[2,1])
  y1 <- c(y1,ci[2,2])
}
z0 <- rep(0,length(x0))
z1 <- rep(0,length(x1))
segments3D(x0,y0,z0,x1,y1,z1,
           zlim = c(0,20),
           phi = 15,
           theta = 30,
           axes = FALSE,
           border = NA,
           box = FALSE,
           r = 0
           )

# step2 : add the event
coords_events <- coordinates(event)
scatter3D(x = coords_events[,1],
          y = coords_events[,2],
          z = rep(0,nrow(event)),
          col="red",
          add=T, pch = 19)

# step3 : add the samples
coords_samples <- coordinates(samples_pts)
scatter3D(x = coords_samples[,1],
          y = coords_samples[,2],
          z = simple_kernel * zfactor,
          col="blue",
          cex=0.1,
          pch = 19,
          add=T)

