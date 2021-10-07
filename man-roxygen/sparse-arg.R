#' @param sparse A Boolean indicating if sparse or regular matrices should be
#'   used by the Rcpp functions. These matrices are used to store edge indices
#'   between two nodes in a graph. Regular matrices are faster, but require more
#'   memory, in particular with multiprocessing. Sparse matrices are slower (a
#'   bit), but require much less memory.
