//  spatial_index.h

// We define here the spatial_index class

#include "spNetwork.h"
#include "boost_importer.h"

#ifndef spatial_index_H
#define spatial_index_H



//******************************************************
// ***** RTREE CLASS **** //
//******************************************************

// ---- Defining the class ----
class spatial_index {

  public:

    // Member variables
    NumericMatrix original_coords;
    MyRtree spindex;

    // constructor
    spatial_index (NumericMatrix x);

    // spatial request method
    IntegerVector tree_request(NumericVector reqBbox);


};

#endif /* spatial_index_H */
