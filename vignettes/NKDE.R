## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
backup_option <- options()
base_wd <- getwd()

## ----include=FALSE------------------------------------------------------------
load(system.file("extdata", "results_vignette_nkde.rda",
                           package = "spNetwork", mustWork = TRUE))

## ----message=FALSE, warning=FALSE---------------------------------------------
library(ggplot2)
library(reshape2)
library(kableExtra)
library(spNetwork)
library(RColorBrewer)

x <- seq(-15.01,15.01,by = 0.01)
df <- data.frame(
  x = x,
  gaussian = gaussian_kernel(x,15),
  epanechnikov = epanechnikov_kernel(x,15),
  quartic = quartic_kernel(x,15),
  triangle = triangle_kernel(x,15),
  tricube = tricube_kernel(x,15),
  triweight = triweight_kernel(x,15),
  cosine = cosine_kernel(x,15),
  uniform = uniform_kernel(x,15))

df2 <- melt(df, id.vars = "x")
names(df2) <- c("x","kernel","y")

ggplot(df2) + 
  geom_line(aes(x=x,y=y,color=kernel),size=1)


## ----message=FALSE, warning=FALSE---------------------------------------------
Funs <- c(gaussian_kernel, gaussian_kernel_scaled,epanechnikov_kernel,
          quartic_kernel,triangle_kernel, uniform_kernel,
          tricube_kernel,triweight_kernel,
          cosine_kernel)
Names <- c("gaussian", "scaled gaussian", "epanechnikov", "quartic",
           "triangle","uniform", "tricube", "triweight", "cosine")
integrals <- sapply(Funs,function(f){
  return(round(integrate(f,upper=15,lower=-15,bw=15)$value,3))
         })
df <- data.frame("kernel"=Names,
                 "integrals" = integrals)

kable(df)


## ----message=FALSE, warning=FALSE---------------------------------------------
x <- seq(-15.01,15.01,by = 0.01)
df <- data.frame(
  x = x,
  gaussian = gaussian_kernel(x,15),
  gaussian_scaled = gaussian_kernel_scaled(x,15),
  epanechnikov = epanechnikov_kernel(x,15),
  quartic = quartic_kernel(x,15)
  )

df2 <- melt(df, id.vars = "x")
names(df2) <- c("x","kernel","y")

ggplot(df2) + 
  geom_line(aes(x=x,y=y,color=kernel),size=1)

## ----message=FALSE, warning=FALSE---------------------------------------------
# first load data and packages
library(sf)
library(spNetwork)
library(tmap)

data(mtl_network) 
data(bike_accidents)

# then plotting the data
tm_shape(mtl_network) + 
  tm_lines("black") + 
  tm_shape(bike_accidents) + 
  tm_dots("red", size = 0.2)

# then calculating some lixels to use as sampling points
lixels <- lixelize_lines(mtl_network,200,mindist = 50)
samples <- lines_center(lixels)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  # then applying the NKDE
#  densities <- nkde(mtl_network,
#                    events = bike_accidents,
#                    w = rep(1,nrow(bike_accidents)),
#                    samples = samples,
#                    kernel_name = "quartic",
#                    bw = 300, div= "bw",
#                    method = "discontinuous", digits = 1, tol = 1,
#                    grid_shape = c(1,1), max_depth = 8,
#                    agg = 5, #we aggregate events within a 5m radius (faster calculation)
#                    sparse = TRUE,
#                    verbose = FALSE)
#  

