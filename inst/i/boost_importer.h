#include <Rcpp.h>
using namespace Rcpp;

// some boost libraries used to building an rtree
#include <b/geometry.hpp>
#include <b/geometry/geometries/geometries.hpp>
#include <b/geometry/geometries/point.hpp>
#include <b/geometry/geometries/box.hpp>
#include <b/geometry/index/rtree.hpp>
#include <b/geometry/algorithms/distance.hpp>
#include <b/foreach.hpp>
#include <b/geometry/index/predicates.hpp>

// **** some namespaces related to boost ****
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// ***** some global types used ****

// basic geometries
typedef bg::model::d2::point_xy<double> point_t;
//typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::linestring<point_t> linestring_t;
typedef bg::model::box<point_t> bbox;

// a simple vector of boost lines
typedef std::vector<linestring_t> lines_vector;

// some simple types for rtree
typedef std::pair<bbox, int> rtree_element;
typedef bgi::rtree< rtree_element, bgi::quadratic<16>> MyRtree;
typedef std::vector<rtree_element> vector_rtree_element ;