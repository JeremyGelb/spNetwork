## ----message=FALSE, warning=FALSE---------------------------------------------
# first load data and packages
library(sp)
library(maptools)
library(rgeos)
library(spNetwork)
library(tmap)
library(raster)

networkgpkg <- system.file("extdata", "networks.gpkg",
                           package = "spNetwork", mustWork = TRUE)

mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network",verbose = FALSE)

# calculating the length of each segment
mtl_network$length <- gLength(mtl_network, byid = TRUE)

# extracting the coordinates in one matrix
coords <- do.call(rbind, unlist(sp::coordinates(mtl_network), recursive = FALSE))
center <- matrix(apply(coords, 2, mean), nrow = 1)
df_center <- data.frame("OID" = 1)
sp::coordinates(df_center) <- center
raster::crs(df_center) <- raster::crs(mtl_network)

# then plotting the data
plot(mtl_network)
plot(df_center, add = TRUE, col = "red")

## ----message=FALSE, warning=FALSE---------------------------------------------
iso_results <- calc_isochrones(lines = mtl_network,
                               start_points = df_center,
                               dists = c(500,1000,2000),
                               weight = "length"
                               )


## ----message=FALSE, warning=FALSE---------------------------------------------
# creating a factor and changing order for better visualisation
iso_results$fac_dist <- as.factor(iso_results$distance)
iso_results <- iso_results[order(-1*iso_results$distance),]

tm_shape(mtl_network) + 
  tm_lines(col = "grey") +
  tm_shape(iso_results) + 
  tm_lines(col = "fac_dist",title.col = "distance (m)",
           palette = c("500"="#005f73", "1000"="#ca6702", "2000"="#9b2226"))+
  tm_layout(legend.outside = TRUE) + 
  tm_shape(df_center) + 
  tm_dots(col = "black", size = 0.1)


## ----message=FALSE, warning=FALSE---------------------------------------------
library(concaveman)

# identidying each isochrone
iso_results$iso_oid <- paste(iso_results$point_id,
                             iso_results$distance,
                             sep = "_")

# creating the polygons for each isochrone
polygons <- lapply(unique(iso_results$iso_oid), function(oid){
  
  # subsetting the required lines
  lines <- subset(iso_results, iso_results$iso_oid == oid)
  
  # extracting the coordinates of the lines
  coords <- do.call(rbind, unlist(sp::coordinates(lines),
                                  recursive = FALSE))
  poly_coords <- concaveman(points = coords, concavity = 3)
  poly <- sp::Polygons(list(sp::Polygon(poly_coords)), ID = oid)
  return(poly)
})

# creating a SpatialPolygonsDataFrame
df <- data.frame(
  iso_oid = unique(iso_results$iso_oid),
  distance = unique(iso_results$distance)
)

iso_sp <- sp::SpatialPolygonsDataFrame(sp::SpatialPolygons(polygons),
                                       df, match.ID = FALSE)
raster::crs(iso_sp) <- raster::crs(mtl_network)

## ----message=FALSE, warning=FALSE---------------------------------------------
# creating a factor and changing order for better visualisation
iso_sp$fac_dist <- as.factor(iso_sp$distance)
iso_sp <- iso_sp[order(-1*iso_sp$distance),]

tm_shape(iso_results) + 
  tm_lines(col = "black")+
  tm_shape(iso_sp) + 
  tm_polygons(col = "fac_dist",title.col = "distance (m)",
           palette = c("500"="#005f73", "1000"="#ca6702", "2000"="#9b2226"),
           alpha = 0.6, border.col = "white") +
  tm_shape(df_center) + 
  tm_dots(col = "black", size = 0.1)


## ----message=FALSE, warning=FALSE---------------------------------------------
library(smoothr)

simple_polygons <- gSimplify(iso_sp, tol = 50)
smoothed_polygons <- smooth(simple_polygons, method = "chaikin", 
                            refinements = 5)

smooth_iso <- sp::SpatialPolygonsDataFrame(smoothed_polygons,
                                           iso_sp@data,
                                           match.ID = FALSE)
tm_shape(iso_results) + 
  tm_lines(col = "black")+
  tm_shape(smooth_iso) +
  tm_polygons(col = "fac_dist",
           title.col = "distance (m)",
           palette = c("500"="#005f73", "1000"="#ca6702", "2000"="#9b2226"),
           border.col = "white",
           alpha = 0.6) +
  tm_shape(df_center) + 
  tm_dots(col = "black", size = 0.1) + 
  tm_layout(legend.outside = TRUE)

