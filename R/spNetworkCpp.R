#' @useDynLib spNetwork
#' @importFrom Rcpp sourceCpp
#' @import methods Rcpp
"_PACKAGE"

.onUnload <- function (libpath) {
  library.dynam.unload("spNetwork", libpath)
}
