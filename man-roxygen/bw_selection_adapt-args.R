#' @param adaptive A boolean indicating if an adaptive bandwidth must be used.
#' If adaptive = TRUE, the local bandwidth are derived from the global bandwidths (bws)
#' @param trim_bws A vector indicating the maximum value an adaptive bandwidth can
#' reach. Higher values will be trimmed. It must have the same length as bws.
#' @param mat_bws A matrix giving the bandwidths for each observation and for each global bandwidth.
#' This is usefull when the user want to use a different method from Abramson's smoothing regimen.