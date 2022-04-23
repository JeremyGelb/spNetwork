## ----message=FALSE, warning=FALSE, message=FALSE------------------------------
# first load data and packages
library(sf)
library(spNetwork)
library(tmap)
library(dbscan)

networkgpkg <- system.file("extdata", "networks.gpkg",
                           package = "spNetwork", mustWork = TRUE)
eventsgpkg <- system.file("extdata", "events.gpkg",
                          package = "spNetwork", mustWork = TRUE)

data(mtl_network)
data(bike_accidents)

# then plotting the data
tm_shape(mtl_network) + 
  tm_lines("black") + 
  tm_shape(bike_accidents) + 
  tm_dots("red", size = 0.2)

## ----include=FALSE------------------------------------------------------------
load(system.file("extdata", "results_vignette_network_build.rda",
                           package = "spNetwork", mustWork = TRUE))

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  # calculating the density values
#  densities <- nkde(mtl_network,
#                    events = bike_accidents,
#                    w = rep(1,nrow(bike_accidents)),
#                    samples = bike_accidents,
#                    kernel_name = "quartic",
#                    bw = 300, div= "bw",
#                    method = "discontinuous", digits = 2, tol = 0.5,
#                    grid_shape = c(1,1), max_depth = 8,
#                    agg = 5,
#                    sparse = TRUE,
#                    verbose = FALSE)
#  
#  bike_accidents$density <- densities * 1000

## ----message=FALSE, warning=FALSE---------------------------------------------
# mapping the density values
tm_shape(mtl_network) + 
  tm_lines(col = "black") +
  tm_shape(bike_accidents) + 
  tm_dots(col = "density", style = "kmeans",
          n = 6, size = 0.1, palette = "-RdYlBu")+
  tm_layout(legend.outside = TRUE) 

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  bike_accidents$weight <- 1
#  agg_points <- aggregate_points(bike_accidents, maxdist = 5)
#  
#  agg_points$OID <- 1:nrow(agg_points)
#  mtl_network$LineID <- 1:nrow(mtl_network)
#  
#  snapped_accidents <- snapPointsToLines2(agg_points,
#                                          mtl_network,
#                                          "LineID")

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  new_lines <- split_lines_at_vertex(mtl_network,
#                                     snapped_accidents,
#                                     snapped_accidents$nearest_line_id,
#                                     mindist = 0.1)

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  new_lines$OID <- 1:nrow(new_lines)
#  new_lines$length <- as.numeric(st_length(new_lines))
#  
#  graph_result <- build_graph(new_lines, 2, "length", attrs = TRUE)

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  btws <- igraph::betweenness(graph_result$graph, directed = FALSE,
#                              normalized = TRUE)
#  vertices <- graph_result$spvertices
#  vertices$btws <- btws

## ----message=FALSE, warning=FALSE---------------------------------------------
# mapping the betweenness
tm_shape(vertices) + 
  tm_dots(col = "btws", style = "kmeans",
          n = 6, size = 0.1, palette = "-RdYlBu")+
  tm_layout(legend.outside = TRUE) 


## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  # first: nn merging between snapped points and nodes
#  xy1 <- st_coordinates(snapped_accidents)
#  xy2 <- st_coordinates(vertices)
#  corr_nodes <- dbscan::kNN(x = xy2, query = xy1, k=1)$id
#  
#  snapped_accidents$btws <- vertices$btws[corr_nodes]
#  
#  # second: nn merging between original points and snapped points
#  xy1 <- st_coordinates(bike_accidents)
#  xy2 <- st_coordinates(snapped_accidents)
#  
#  corr_nodes <- dbscan::kNN(x = xy2, query = xy1, k=1)$id
#  bike_accidents$btws <- snapped_accidents$btws[corr_nodes]

## ----message=FALSE, warning=FALSE---------------------------------------------
# mapping the results
tm_shape(bike_accidents) + 
  tm_dots(col = "btws", style = "kmeans",
          n = 6, size = 0.1, palette = "-RdYlBu")+
  tm_layout(legend.outside = TRUE)

tm_shape(bike_accidents) + 
  tm_dots(col = "density", style = "kmeans",
          n = 6, size = 0.1, palette = "-RdYlBu")+
  tm_layout(legend.outside = TRUE) 

## ----message=FALSE, warning=FALSE---------------------------------------------
cor.test(bike_accidents$density, bike_accidents$btws)

