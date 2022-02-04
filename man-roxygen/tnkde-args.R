#' @param time_field The name of the field in events indicating when the events
#'   occurred. It must be a numeric field
#' @param bw_net The network kernel bandwidth (using the scale of the lines),
#'   can be a single float or a numeric vector if a different bandwidth must be
#'   used for each event.
#' @param bw_time The time kernel bandwidth, can be a single float or a numeric
#'   vector if a different bandwidth must be used for each event.
#' @param adaptive A Boolean, indicating if an adaptive bandwidth must be used.
#'   Both spatial and temporal bandwidths are adapted but separately.
#' @param trim_bw_net A float, indicating the maximum value for the adaptive
#'   network bandwidth
#' @param trim_bw_time A float, indicating the maximum value for the adaptive
#'   time bandwidth
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#'   events within the radius around each sampling point), "bw" (the bandwith)
#'   "none" (the simple sum).
