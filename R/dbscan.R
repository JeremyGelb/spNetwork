
# for the dbscan algorithm

# pre-traitement :
# ajuster chacun des evenements sur le reseau
# ajouter ces evenements comme vertex sur le reseau
# generer un graph complet

networkgpkg <- system.file("extdata", "networks.gpkg",
                           package = "spNetwork", mustWork = TRUE)
eventsgpkg <- system.file("extdata", "events.gpkg",
                          package = "spNetwork", mustWork = TRUE)

lines <- rgdal::readOGR(networkgpkg,layer="mtl_network",verbose = FALSE)
events <- rgdal::readOGR(eventsgpkg,layer="bike_accidents", verbose = FALSE)

w <- rep(1,nrow(events))
digits = 2
tol = 0.1
agg = 5

# utiliser pour cela la lib https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/index.html

generate_network <- function(lines, events, w, digits, tol, agg){

  ## step1 cleaning the events
  events$weight <- w
  events <- clean_events(events,digits,agg)

  ## step2 defining the global IDS
  events$goid <- 1:nrow(events)

  ## step3 remove lines with no length
  lines$length <- gLength(lines,byid = TRUE)
  lines <- subset(lines, lines$length>0)
  lines$oid <- 1:nrow(lines)

  ## step4 snapping the points on the lines
  snapped_pts <- snapPointsToLines2(events, lines, "oid")
  new_lines <- add_vertices_lines(lines,snapped_pts,snapped_pts$nearest_line_id,tol)

  ## step5 make them simple lines
  new_lines <- simple_lines(new_lines)
  new_lines$length <- gLength(new_lines,byid = TRUE)
  new_lines <- subset(new_lines,new_lines$length>0)

  # step56 remove lines that are loops
  new_lines <- remove_loop_lines(new_lines,digits)
  new_lines$oid <- 1:nrow(new_lines)
  new_lines <- new_lines[c("length","oid")]

  ## step7 creating the graph
  graph_result <- build_graph(lines,digits = digits,line_weight = "length")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  edges <- graph_result$spedges

  ## step8 finding for each event, its corresponding node
  events$vertex_id <- closest_points(events, nodes)
  graph_result$events <- events
  return(graph_result)
}

dbscan_network <- function(graph_result, bw, minpts){
  edge_list <- graph_result$linelist

}

