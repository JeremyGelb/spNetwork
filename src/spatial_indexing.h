//  spatial_index.h

// We define here the spatial_index class

#ifndef spatial_index_H
#define spatial_index_H

#include "spNetwork.h"
#include "boost_importer.h"


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
