#' @param adaptive A boolean indicating if an adaptive bandwidth must be used.
#' If adaptive = TRUE, the local bandwidth are derived from the global bandwidths
#' calculated from bw_range and bw_step.
#' @param trim_bws A vector indicating the maximum value an adaptive bandwidth can
#' reach. Higher values will be trimmed. It must have the same length as
#' seq(bw_range[[1]],bw_range[[2]], bw_step).
#' @param mat_bws A matrix giving the bandwidths for each observation and for each global bandwidth.
#' This is usefull when the user want to use a different method from Abramson's smoothing regimen.
