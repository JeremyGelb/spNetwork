#' @title Adaptive bw in one dimension
#'
#' @description Calculate adaptive bandwidths in one dimension
#'
#' @param events A numeric vector representing the moments of occurrence of events
#' @param w The weight of the events
#' @param samples A numeric vector representing the moments to sample
#' @param bw A float, the bandiwdth to use
#' @param kernel_name The name of the kernel to use
#' @keywords internal
adaptive_bw_1d <- function(events, w, bw, kernel_name){

  kern_func <- select_kernel(kernel_name)
  dens <- sapply(events, function(s){
    dists <- abs(s-events)
    return(sum(kern_func(dists,bw)*w))
  }) / bw

  lambdaf <- exp(sum(log(1/sqrt(dens)))/length(events))
  bws <- bw * (1/sqrt(dens)) * (1/lambdaf)
  return(bws)
}



#' @title Temporal Kernel density estimate
#'
#' @description Calculate the Temporal kernel density estimate based on sampling points in
#' time and events
#' @param events A numeric vector representing the moments of occurrence of events
#' @param w The weight of the events
#' @param samples A numeric vector representing the moments to sample
#' @param bw A float, the bandwidth to use
#' @param kernel_name The name of the kernel to use
#' @param adaptive Boolean
#' @export
#' @examples
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' bike_accidents <- sf::st_read(eventsgpkg,layer="bike_accidents")
#' bike_accidents$Date <- as.POSIXct(bike_accidents$Date, format = "%Y/%m/%d")
#' start <- min(bike_accidents$Date)
#' diff <- as.integer(difftime(bike_accidents$Date , start, units = "days"))
#' density <- tkde(diff, rep(1,length(diff)), seq(0,max(diff),1), 2, "quartic")
tkde <- function(events, w, samples, bw, kernel_name, adaptive = FALSE){

  ## selecting the kernel function
  kern_func <- select_kernel(kernel_name)

  ## calculating the adaptive bw if required
  if (adaptive){

    bws <- adaptive_bw_1d(events, bw, kern_func)

  }else{
    bws <- rep(bw, length(events))
  }

  ## calculating the kernel values with the bw
  kernel_values <- sapply(samples, function(s){
    dists <- abs(s-events)
    return(sum(kern_func(dists,bws)*w))
  })

  return((1/bw)*kernel_values)

}


#' @title Bandwidth selection for Temporal Kernel density estimate by likelihood cross validation
#'
#' @description Calculate the likelihood cross validation score for several bandwidths for the
#' Temporal Kernel density
#' @param events A numeric vector representing the moments of occurrence of events
#' @param w The weight of the events
#' @param bws A numeric vector, the bandiwdths to use
#' @param kernel_name The name of the kernel to use
#' @export
#' @examples
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' bike_accidents <- sf::st_read(eventsgpkg,layer="bike_accidents")
#' bike_accidents$Date <- as.POSIXct(bike_accidents$Date, format = "%Y/%m/%d")
#' start <- min(bike_accidents$Date)
#' diff <- as.integer(difftime(bike_accidents$Date , start, units = "days"))
#' w <- rep(1,length(diff))
#' scores <- bw_cv_likelihood_calc_tkde(diff, w, seq(10,60,10), "quartic")
bw_cv_likelihood_calc_tkde <- function(events, w, bws, kernel_name){

  kern_func <- select_kernel(kernel_name)
  n <- length(events)

  scores <- sapply(bws, function(bw){
    vals <- sapply(1:length(events), function(i){
      ei <- events[[i]]
      wi <- w[[i]]
      dists <- abs(ei-events)
      dens <- kern_func(dists,bw)*w
      tot_dens <- ((sum(dens) - (kern_func(0,bw)*wi)))
      return(tot_dens)
    })
    score <- ((sum(log(vals)))/n) - log((n-1)*bw)
    return(score)
  })
  return(scores)
}
