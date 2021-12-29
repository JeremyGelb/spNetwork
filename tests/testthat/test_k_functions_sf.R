context("testing functions for k and g function analysis")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE SIMPLE K AND G FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the simple k function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))
  event <- st_as_sf(event, coords = c("x","y"))

  # calculating the observed values
  observed <- kfunctions(all_lines, event, 0, 6, 0.5, 2, 5, conf_int = 0.05, digits = 2, tol = 0.1, resolution = NULL, agg = NULL, verbose = TRUE)

  # after checking on a paper with a pen, the observed k and g values at distance 3 must be :
  expected_vals <- c(0.9, 1.5)
  diff <- observed$values[7,c("obs_k","obs_g")] - expected_vals
  diff <- round(sum(abs(diff)),10)
  expect_equal(diff, 0)

})



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE cross K AND G FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the cross k function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))

  event <- st_as_sf(event, coords = c("x","y"))

  As <- event[c(1,2),]
  Bs <- event[c(3,4),]

  # calculating the observed values

  observed <- cross_kfunctions(all_lines, As, Bs, 0, 6, 0.5, 2, 5, conf_int = 0.05, digits = 2, tol = 0.1, resolution = NULL, agg = NULL, verbose = TRUE)

  # after checking on a paper with a pen, the observed k and g values at distance 3 must be :
  expected_vals <- c(0.2, 0.4)
  diff <- observed$values[7,c("obs_k","obs_g")] - expected_vals
  diff <- round(sum(abs(diff)),10)
  expect_equal(diff, 0)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE FIRST RANDOMIZATION MATRIX FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the first randomization function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))
  event <- st_as_sf(event, coords = c("x","y"))

  # creating the graph on wich we will do the randomization
  lines <- all_lines
  points <- event
  agg <- NULL
  digits <- 2
  tol <- 0.1

  ## step0 : clean the points
  n <- nrow(points)
  points$goid <- seq_len(nrow(points))
  points$weight <- rep(1,nrow(points))
  points <- clean_events(points,digits,agg)

  probs <- NULL

  ## step1 : clean the lines
  if(is.null(probs)){
    lines$probs <- 1
  }else{
    lines$probs <- probs
  }

  lines$length <- as.numeric(st_length(lines))
  lines <- subset(lines, lines$length>0)
  lines$oid <- seq_len(nrow(lines))


  ## step2 : adding the points to the lines
  snapped_events <- snapPointsToLines2(points,lines,idField = "oid")
  new_lines <- split_lines_at_vertex(lines, snapped_events,
                                     snapped_events$nearest_line_id, tol)

  ## step3 : splitting the lines
  new_lines$length <- as.numeric(st_length(new_lines))
  new_lines <- subset(new_lines,new_lines$length>0)
  new_lines <- remove_loop_lines(new_lines,digits)
  new_lines$oid <- seq_len(nrow(new_lines))
  new_lines <- new_lines[c("length","oid","probs")]
  Lt <- sum(as.numeric(st_length(new_lines)))

  ## step4 : building the graph for the real case
  graph_result <- build_graph(new_lines,digits = digits,
                              line_weight = "length",
                              attrs = TRUE)
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
                                                           name = "probs")
  snapped_events$vertex_id <- closest_points(snapped_events, nodes)

  ## now time to test !
  set.seed(123)
  observed <- randomize_distmatrix(graph_result$graph,
                                   st_drop_geometry(graph_result$spedges),
                                   n = nrow(event))

  # NOTE : this result is obtained after checking for the function
  # by hand. It is hard to test it other wise because
  # we select the location of the points randomly
  expected <- rbind(
    c(0.0000000, 2.968420, 0.5314272, 1.1556292),
    c(2.9684197, 0.000000, 3.3807815, 4.1240489),
    c(0.5314272, 3.380782, 0.0000000, 0.7432674),
    c(1.1556292, 4.124049, 0.7432674, 0.0000000)
  )

  result <- sum(round(abs(observed-expected),6))
  expect_equal(result,0)

})




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE SECOND RANDOMIZATION MATRIX FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the second randomization function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))
  event <- st_as_sf(event, coords = c("x","y"))

  # creating the graph on wich we will do the randomization
  lines <- all_lines
  points <- event
  agg <- NULL
  digits <- 2
  tol <- 0.1

  ## step0 : clean the points
  n <- nrow(points)
  points$goid <- seq_len(nrow(points))
  points$weight <- rep(1,nrow(points))
  points <- clean_events(points,digits,agg)

  probs <- NULL

  ## step1 : clean the lines
  if(is.null(probs)){
    lines$probs <- 1
  }else{
    lines$probs <- probs
  }

  lines$length <- as.numeric(st_length(lines))
  lines <- subset(lines, lines$length>0)
  lines$oid <- seq_len(nrow(lines))


  ## step2 : adding the points to the lines
  snapped_events <- snapPointsToLines2(points,lines,idField = "oid")
  new_lines <- split_lines_at_vertex(lines, snapped_events,
                                     snapped_events$nearest_line_id, tol)

  ## step3 : splitting the lines
  new_lines$length <- as.numeric(st_length(new_lines))
  new_lines <- subset(new_lines,new_lines$length>0)
  new_lines <- remove_loop_lines(new_lines,digits)
  new_lines$oid <- seq_len(nrow(new_lines))
  new_lines <- new_lines[c("length","oid","probs")]
  Lt <- sum( as.numeric(st_length(new_lines)))
  #new_lines$weight <- gLength(new_lines, byid = TRUE)

  ## step4 : building the graph for the real case
  graph_result <- build_graph(new_lines,digits = digits,
                              line_weight = "length",
                              attrs = TRUE)
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  graph_result$spedges$probs <- igraph::get.edge.attribute(graph_result$graph,
                                                           name = "probs")
  snapped_events$vertex_id <- closest_points(snapped_events, nodes)

  ## now time to test !
  observed <- randomize_distmatrix2(graph_result$graph,
                                   st_drop_geometry(graph_result$spedges),
                                   n = nrow(event),
                                   resolution = 0.5,
                                   nsim = 5)

  ## first of all, the matrix should have only 0s in the diag
  test1 <- sapply(observed, function(mat){
    sum(diag(mat)) == 0
  })

  ## second, all the matrices must be symetric
  test2 <- sapply(observed, function(mat){
    isSymmetric.matrix(mat)
  })

  ## third, the distance between two points must be at max : 10
  test3 <- sapply(observed, function(mat){
    max(mat) <= 10
  })

  total_test <- any(c(!test1,!test2,!test3))

  expect_false(total_test)

})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR THE MULTICORE FUNCTIONS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the multicore simple k function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))
  event <- st_as_sf(event, coords = c("x","y"))

  # calculating the observed values
  future::plan(future::multisession(workers=1))
  observed <- kfunctions.mc(all_lines, event, 0, 6, 0.5, 2, 50, conf_int = 0.05, digits = 2, tol = 0.1, resolution = NULL, agg = NULL, verbose = TRUE)

  # after checking on a paper with a pen, the observed k and g values at distance 3 must be :
  expected_vals <- c(0.9, 1.5)
  diff <- observed$values[7,c("obs_k","obs_g")] - expected_vals
  diff <- round(sum(abs(diff)),10)

  expect_equal(diff, 0)

})


test_that("Testing the multicore cross k function", {

  # defining a simple situation
  wkt_lines <- c(
    "LINESTRING (0 5, 0 0)",
    "LINESTRING (-5 0, 0 0)",
    "LINESTRING (0 -5, 0 0)",
    "LINESTRING (5 0, 0 0)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of three events
  event <- data.frame(x=c(0,3,1,0),
                      y=c(3,0,0,1),
                      id = c(1,2,3,4))
  event <- st_as_sf(event, coords = c("x","y"))

  As <- event[c(1,2),]
  Bs <- event[c(3,4),]

  # calculating the observed values
  future::plan(future::multisession(workers=1))
  observed <- cross_kfunctions.mc(all_lines, As, Bs, 0, 6, 0.5, 2, 5, conf_int = 0.05, digits = 2, tol = 0.1, resolution = NULL, agg = NULL, verbose = TRUE)

  # after checking on a paper with a pen, the observed k and g values at distance 3 must be :
  expected_vals <- c(0.2, 0.4)
  diff <- observed$values[7,c("obs_k","obs_g")] - expected_vals
  diff <- round(sum(abs(diff)),10)
  expect_equal(diff, 0)

})


