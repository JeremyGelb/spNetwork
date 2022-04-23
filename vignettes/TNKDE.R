## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
backup_option <- options()
base_wd <- getwd()
library(ggplot2)

## ----include=FALSE------------------------------------------------------------
load(system.file("extdata", "results_vignette_tnkde.rda",
                           package = "spNetwork", mustWork = TRUE))

## ----message=FALSE, warning=FALSE, out.width = "70%", message=FALSE-----------
# first load data and packages
library(spNetwork)
library(tmap)
library(sf)

data(bike_accidents)

# converting the Date field to a numeric field (counting days)
bike_accidents$Time <- as.POSIXct(bike_accidents$Date, format = "%Y/%m/%d")
start <- as.POSIXct("2016/01/01", format = "%Y/%m/%d")
bike_accidents$Time <- difftime(bike_accidents$Time, start, units = "days")
bike_accidents$Time <- as.numeric(bike_accidents$Time)

months <- as.character(1:12)
months <- ifelse(nchar(months)==1, paste0("0", months), months)
months_starts_labs <- paste("2016/",months,"/01", sep = "")
months_starts_num <- as.POSIXct(months_starts_labs, format = "%Y/%m/%d")
months_starts_num <- difftime(months_starts_num, start, units = "days")
months_starts_num <- as.numeric(months_starts_num)
months_starts_labs <- gsub("2016/", "", months_starts_labs, fixed = TRUE)

ggplot(bike_accidents) + 
  geom_histogram(aes(x = Time), bins = 30, color = "white") + 
  scale_x_continuous(breaks = months_starts_num, labels = months_starts_labs)

## ----message=FALSE, warning=FALSE---------------------------------------------
bike_accidents <- subset(bike_accidents, bike_accidents$Time >= 90)

## ----message=FALSE, warning=FALSE, out.width = "80%"--------------------------
w <- rep(1,nrow(bike_accidents))
samples <- seq(0, max(bike_accidents$Time), 0.5)

time_kernel_values <- data.frame(
  bw_10 = tkde(bike_accidents$Time, w = w, samples = samples, bw = 10, kernel_name = "quartic"),
  bw_20 = tkde(bike_accidents$Time, w = w, samples = samples, bw = 20, kernel_name = "quartic"),
  bw_30 = tkde(bike_accidents$Time, w = w, samples = samples, bw = 30, kernel_name = "quartic"),
  bw_40 = tkde(bike_accidents$Time, w = w, samples = samples, bw = 40, kernel_name = "quartic"),
  bw_50 = tkde(bike_accidents$Time, w = w, samples = samples, bw = 50, kernel_name = "quartic"),
  bw_60 = tkde(bike_accidents$Time, w = w, samples = samples, bw = 60, kernel_name = "quartic"),
  time = samples
)

df_time <- reshape2::melt(time_kernel_values,id.vars = "time")
df_time$variable <- as.factor(df_time$variable)

ggplot(data = df_time) + 
  geom_line(aes(x = time, y = value)) +
  scale_x_continuous(breaks = months_starts_num, labels = months_starts_labs) +
  facet_wrap(vars(variable), ncol=2, scales = "free") + 
  theme(axis.text = element_text(size = 5))

## ----message=FALSE, warning=FALSE, out.width = "80%"--------------------------
bw1 <- bw.bcv(bike_accidents$Time, nb = 1000, lower = 1, upper = 80)
bw2 <- bw.ucv(bike_accidents$Time, nb = 1000, lower = 1, upper = 80)
bw3 <- bw.SJ(bike_accidents$Time, nb = 1000, lower = 1, upper = 80)


time_kernel_values <- data.frame(
  bw_bcv = tkde(bike_accidents$Time, w = w, samples = samples, bw = bw1, kernel_name = "quartic"),
  bw_ucv = tkde(bike_accidents$Time, w = w, samples = samples, bw = bw2, kernel_name = "quartic"),
  bw_SJ = tkde(bike_accidents$Time, w = w, samples = samples, bw = bw3, kernel_name = "quartic"),
  time = samples
)

df_time <- reshape2::melt(time_kernel_values,id.vars = "time")
df_time$variable <- as.factor(df_time$variable)

ggplot(data = df_time) + 
  geom_line(aes(x = time, y = value)) +
  scale_x_continuous(breaks = months_starts_num, labels = months_starts_labs) +
  facet_wrap(vars(variable), ncol=2, scales = "free")  + 
  theme(axis.text = element_text(size = 5))


## ----message=FALSE, warning=FALSE, out.width = "75%", fig.align='center', message = FALSE----
# loading the road network
data(mtl_network)

tm_shape(mtl_network) + 
  tm_lines(col = "black") + 
  tm_shape(bike_accidents) + 
  tm_dots(col = "red", size = 0.1)


## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  # creating sample points
#  lixels <- lixelize_lines(mtl_network, 50)
#  sample_points <- lines_center(lixels)
#  
#  # calculating the densities
#  nkde_densities <- nkde(lines = mtl_network,
#                         events = bike_accidents,
#                         w = rep(1,nrow(bike_accidents)),
#                         samples = sample_points,
#                         kernel_name = "quartic",
#                         bw = 450,
#                         adaptive = TRUE, trim_bw = 900,
#                         method = "discontinuous",
#                         div = "bw",
#                         max_depth = 10,
#                         digits = 2, tol = 0.1, agg = 5,
#                         grid_shape = c(1,1),
#                         verbose = FALSE)
#  
#  sample_points$density <- nkde_densities$k * 1000
#  

## ----message=FALSE, warning=FALSE, out.width="75%", fig.align='center'--------
tm_shape(sample_points) + 
  tm_dots(col = "density", style = "kmeans", n = 8, palette = "viridis", size = 0.05) + 
  tm_layout(legend.outside = TRUE)

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  cv_scores <- bw_tnkde_cv_likelihood_calc(
#    bw_net_range = c(100,1000),
#    bw_net_step = 100,
#    bw_time_range = c(10,60),
#    bw_time_step = 10,
#    lines = mtl_network,
#    events = bike_accidents,
#    time_field = "Time",
#    w = rep(1, length(bike_accidents)),
#    kernel_name = "quartic",
#    method = "discontinuous",
#    diggle_correction = FALSE,
#    study_area = NULL,
#    max_depth = 10,
#    digits = 2,
#    tol = 0.1,
#    agg = 15,
#    sparse=TRUE,
#    grid_shape=c(1,1),
#    sub_sample=1,
#    verbose = FALSE,
#    check = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
knitr::kable(cv_scores)

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  # choosing sample in times (every 10 days)
#  sample_time <- seq(0, max(bike_accidents$Time), 10)
#  
#  # calculating densities
#  tnkde_densities <- tnkde(lines = mtl_network,
#                     events = bike_accidents,
#                     time_field = "Time",
#                     w = rep(1, nrow(bike_accidents)),
#                     samples_loc = sample_points,
#                     samples_time = sample_time,
#                     kernel_name = "quartic",
#                     bw_net = 700, bw_time = 60,
#                     adaptive = TRUE,
#                     trim_bw_net = 900,
#                     trim_bw_time = 80,
#                     method = "discontinuous",
#                     div = "bw", max_depth = 10,
#                     digits = 2, tol = 0.01,
#                     agg = 15, grid_shape = c(1,1),
#                     verbose  = FALSE)
#  
#  # creating a color palette for all the densities
#  library(classInt)
#  library(viridis)
#  all_densities <- c(tnkde_densities$k)
#  color_breaks <- classIntervals(all_densities, n = 10, style = "kmeans")
#  
#  # generating a map at each sample time
#  all_maps <- lapply(1:ncol(tnkde_densities$k), function(i){
#    time <- sample_time[[i]]
#    date <- as.Date(start) + time
#  
#    sample_points$density <- tnkde_densities$k[,i]
#    map1 <- tm_shape(sample_points) +
#    tm_dots(col = "density", size = 0.01,
#            breaks = color_breaks$brks, palette = viridis(10)) +
#      tm_layout(legend.show=FALSE, main.title = as.character(date), main.title.size = 0.5)
#    return(map1)
#  })
#  
#  # creating a gif with all the maps
#  tmap_animation(all_maps, filename = "images/animated_map.gif",
#                 width = 1000, height = 1000, dpi = 300, delay = 50)

## ----message=FALSE, warning=FALSE, out.width = "75%", fig.align='center'------
knitr::include_graphics("images/animated_map.gif")

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  tnkde_densities <- tnkde(lines = mtl_network,
#                     events = bike_accidents,
#                     time_field = "Time",
#                     w = rep(1, nrow(bike_accidents)),
#                     samples_loc = sample_points,
#                     samples_time = sample_time,
#                     kernel_name = "quartic",
#                     bw_net = 700, bw_time = 60,
#                     adaptive = TRUE,
#                     adaptive_separate = FALSE,
#                     trim_bw_net = 900,
#                     trim_bw_time = 80,
#                     method = "discontinuous",
#                     div = "bw", max_depth = 10,
#                     digits = 2, tol = 0.01,
#                     agg = 15, grid_shape = c(1,1),
#                     verbose  = FALSE)

