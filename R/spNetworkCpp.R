#' @useDynLib spNetwork
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("spNetwork", libpath)
}
