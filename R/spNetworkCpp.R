#' @useDynLib spNetwork
#' @importFrom Rcpp sourceCpp
#' @import RcppProgress
#' @import RcppArmadillo
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("spNetwork", libpath)
}
