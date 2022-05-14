// Include our definition of the student file (e.g. "")
#include "spatial_indexing.h"

//******************************************************
// ***** RTREE CLASS EXPORT VIA RCPP **** //
//******************************************************

// How to do it properly : https://github.com/r-pkg-examples/rcpp-modules-student

RCPP_EXPOSED_CLASS(Double)
RCPP_MODULE(spatial_index_cpp) {

  Rcpp::class_<spatial_index>("spatial_index")
    .constructor<Rcpp::NumericMatrix>()
    .method("tree_request", &spatial_index::tree_request)
  ;
}
