//  spatial_index.h

// We define here the spatial_index class

#ifndef spatial_index_H
#define spatial_index_H

#include "spNetwork.h"

// some boost libraries used to building an rtree
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry/index/predicates.hpp>

// **** some namespaces related to boost ****
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// ***** some global types used ****

// basic geometries
typedef bg::model::d2::point_xy<double> point_t;
typedef bg::model::linestring<point_t> linestring_t;
typedef bg::model::box<point_t> bbox;

// a simple vector of boost lines
typedef std::vector<linestring_t> lines_vector;

// some simple types for rtree
typedef std::pair<bbox, int> rtree_element;
typedef bgi::rtree< rtree_element, bgi::quadratic<16>> MyRtree;
typedef std::vector<rtree_element> vector_rtree_element ;

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
