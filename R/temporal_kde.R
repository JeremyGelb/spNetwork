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
#' @param bw A float, the bandiwdth to use
#' @param kernel_name The name of the kernel to use
#' @param adaptive Boolean
#' @examples
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' bike_accidents <- rgdal::readOGR(eventsgpkg,layer="bike_accidents", verbose=FALSE)
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

  return(kernel_values)

}
