
generate_random_points <- function(n,lines,nsim){
  #first, selecting the lines on which we will have the new points
  lines$length <- gLength(lines,byid=T)
  lines$prob <- lines$length/sum(lines$length)
  lines$oid <- 1:nrow(lines)
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  random_events <- lapply(1:nsim,function(i){
    setTxtProgressBar(pb, i)
    ids <- sample(x = lines$oid,
                  size = n,
                  replace=T,
                  prob=lines$prob)
    sel_lines <- lines[ids,]
    #then select for each line the distance from the first point
    nline <- nrow(sel_lines)
    dists <- runif(nline,min=rep(0,nline),max=sel_lines$length)
    #and generate the points
    pts <- lapply(1:nline,function(i){
      line <- sel_lines[i,]
      pt <- gInterpolate(line,dists[[i]])
      return(coordinates(pt))
    })
    pts <- data.frame(do.call(rbind,pts))
    pts$lines_goid <- sel_lines$goid
    coordinates(pts) <- pts[c("x","y")]
    crs(pts) <- crs(lines)
    return(pts)
  })
  return(random_events)
}


nkfunction <- function(events,w,lines,start,end,step,nsim,grid_shape,digits,tol,ci=0.95,test = 'lower.tail'){
  #step 0 : cleaning the events
  events$weight <- w
  events <- clean_events(events,digits,tol)
  lines$length <- gLength(lines,byid=TRUE)
  if(start==0){
    start <- start+step
  }

  #step 1 : generating the random set of points
  n <- nrow(events)
  lines$goid <- 1:nrow(lines)
  print("generating the random points sets...")
  random_events <- generate_random_points(n,lines,nsim)

  #step 2 calculate the k function fo each random set of points
  print("calculation of the random permutations...")
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  i<-0
  k_values_sim <- lapply(random_events,function(pt_set){
    setTxtProgressBar(pb, i)
    i<<- i+1
    pt_set$weight <- 1
    k_val <- kfunction_worker(pt_set,lines,start,end,step,grid_shape, lines_goid = pt_set$lines_goid)
  })
  k_val_table <- do.call(rbind,k_values_sim)
  #step 3 calculate the true k values
  k_val <- kfunction_worker(events, lines,start,end,step,grid_shape)

  # #step 4 : calculate pseudo p-values
  # if(test == 'lower.tail'){
  #   test_mat <- k_val <= k_val_table
  # }
  # if(test == 'upper.tail'){
  #   test_mat <- k_val >= k_val_table
  # }
  # test_count <- colSums(test_mat)
  # p_values <- test_count / nsim
  #
  # #step5 : generate chart
  # df <- data.table(
  #   "k" = k_val,
  #   "dists" = breaks <- seq(start,end,step),
  #   "lower" = apply(k_val_table,2,quantile,probs=1-ci),
  #   "upper" = apply(k_val_table,2,quantile,probs=ci)
  # )
  # ggplot(data = df)+
  #   geom_ribbon(aes(x=dists,ymin=lower,ymax = upper),fill=rgb(0.2,0.2,0.2,0.2))+
  #   geom_path(aes(x=dists,y=k),color="blue")

}



kfunction_worker <- function(events,lines,start,end,step,grid_shape, lines_goid = NULL){

  Lt <- sum(lines$length)
  n <- sum(events$weight)
  events$line_goid <- lines_goid

  #step1 : selecting the lines and the events in each quadra
  events$goid <- 1:nrow(events)
  grid <- build_grid(grid_shape,list(events,lines))
  tree_events <- build_quadtree(events)
  tree_lines <- build_quadtree(lines)

  selections <- lapply(1:length(grid),function(i){
    quadra <- grid[i,]
    base_events <- spatial_request(quadra,tree_events,events)
    if(nrow(base_events)==0){
      return(NULL)
    }
    buff <- gBuffer(base_events,width = end)
    sel_lines <-  spatial_request(buff,tree_lines,lines)
    sel_events <- spatial_request(buff,tree_events,events)
    return(list("events" = base_events,
                "lines" = sel_lines,
                "destinations" = sel_events))
  })
  selections <- selections[lengths(selections) != 0]

  #iterating over the selections
  results <- lapply(selections,function(sel){

    #extracting elements
    sel_events <- sel$events
    w1 <- sel_events$weight
    sel_lines <- sel$lines
    sel_destinations <- sel$destinations
    w2 <- sel_destinations$weight

    #splitting the lines at the event points
    sel_lines$oid <- 1:nrow(sel_lines)
    if(is.null(lines_goid)){
      snapped_events <- snapPointsToLines(sel_destinations,sel_lines,idField = "oid")
      sel_destinations$nearest_line_id <- snapped_events$nearest_line_id
    }else{
      sel_destinations$nearest_line_id <- left_join(sel_destinations@data,
                                                    sel_lines@data[c("oid","goid")],
                                                    by = c("lines_goid"="goid")
                                                    )$oid
    }
    new_lines <- add_vertices_lines(sel_lines,sel_destinations,sel_destinations$nearest_line_id,tol)
    new_lines <- simple_lines(new_lines)
    new_lines$length <- gLength(new_lines,byid = T)
    new_lines <- subset(new_lines,new_lines$length>0)

    #generating the graph
    graph_result <- build_graph(new_lines,digits = digits,line_weight = "length")
    graph <- graph_result$graph
    nodes <- graph_result$spvertices
    edges <- graph_result$spedges

    #finding for each event, its corresponding node
    sel_destinations$vertex_id <- closest_points(sel_destinations, nodes)

    sel_events$vertex_id <- left_join(sel_destinations@data[c("goid","vertex_id")],
                                      sel_events@data,by=c("goid" = "goid"))$vertex_id
    #calculating the distance matrix
    dist_mat <- igraph::distances(graph,sel_events$vertex_id,sel_destinations$vertex_id,mode = "out")

    #calculating the vector value for each event
    k_mat <- kfunction(dist_mat,start,end,step, w1, w2)
    k_tot <- colSums(k_mat)
    return(k_tot)

  })

  #combining the results
  k_values <- do.call(rbind,results)
  k_values <- colSums(k_values)

  #calculating the final values (a vector ok k_hat, one for each distance)
  k_values <- (Lt/(n*(n-1))) * (k_values/n)
  return(k_values)

}


kfunction <- function(dist_mat,start,end,step, w1, w2){
  breaks <- seq(start,end,step)
  cols <- lapply(breaks,function(dist){
    int_mat <- ifelse(dist_mat<=dist,1,0) * w2
    tot <- (rowSums(int_mat) - 1)* w1  #removing the origin point and multiplying by weight
    return(tot)
  })
  k_mat <- do.call(cbind,cols)
  return(k_mat)
}


