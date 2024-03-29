#' @param bws_net An ordered numeric vector with all the network bandwidths
#' @param bws_time An ordered numeric vector with all the time bandwidths
#' @param adaptive A boolean indicating if local bandwidths must be calculated
#' @param arr_bws_net An array with all the local netowrk bandwidths precalculated (for each event, and at each possible combinaison of network and temporal bandwidths). The dimensions must be c(length(net_bws), length(time_bws), nrow(events)))
#' @param arr_bws_time An array with all the local time bandwidths precalculated (for each event, and at each possible combinaison of network and temporal bandwidths). The dimensions must be c(length(net_bws), length(time_bws), nrow(events)))
#' @param trim_net_bws A numeric vector with the maximum local network bandwidth. If local bandwidths have higher values, they will be replaced by the corresponding value in this vector.
#' @param trim_time_bws A numeric vector with the maximum local time bandwidth. If local bandwidths have higher values, they will be replaced by the corresponding value in this vector.
