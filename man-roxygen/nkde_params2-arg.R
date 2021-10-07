#' @param bw The kernel bandwidth (using the scale of the lines)
#' @param adaptive A Boolean, indicating if an adaptive bandwidth must be
#' used
#' @param trim_bw A float, indicating the maximum value for the adaptive
#' bandwidth
#' @param div The divisor to use for the kernel. Must be "n" (the number of
#' events within the radius around each sampling point), "bw" (the bandwith)
#' "none" (the simple sum).
