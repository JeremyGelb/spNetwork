// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2015 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2015 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2015 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2014-2021.
// Modifications copyright (c) 2014-2021 Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle
// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_GEOMETRY_HPP
#define BOOST_GEOMETRY_GEOMETRY_HPP

#include <b/config.hpp>

#if defined(BOOST_NO_CXX14_CONSTEXPR)
#error "Use C++14 or higher to compile Boost.Geometry, or use Boost 1.72 or lower."
#endif

// Shortcut to include all header files

#include <b/geometry/core/closure.hpp>
#include <b/geometry/core/coordinate_dimension.hpp>
#include <b/geometry/core/coordinate_system.hpp>
#include <b/geometry/core/coordinate_type.hpp>
#include <b/geometry/core/cs.hpp>
#include <b/geometry/core/geometry_types.hpp>
#include <b/geometry/core/interior_type.hpp>
#include <b/geometry/core/point_order.hpp>
#include <b/geometry/core/point_type.hpp>
#include <b/geometry/core/ring_type.hpp>
#include <b/geometry/core/tag.hpp>
#include <b/geometry/core/tag_cast.hpp>
#include <b/geometry/core/tags.hpp>
#include <b/geometry/core/visit.hpp>

// Core algorithms
#include <b/geometry/core/access.hpp>
#include <b/geometry/core/exterior_ring.hpp>
#include <b/geometry/core/interior_rings.hpp>
#include <b/geometry/core/radian_access.hpp>
#include <b/geometry/core/radius.hpp>
#include <b/geometry/core/topological_dimension.hpp>

#include <b/geometry/arithmetic/arithmetic.hpp>
#include <b/geometry/arithmetic/dot_product.hpp>

#include <b/geometry/strategies/strategies.hpp>

#include <b/geometry/algorithms/append.hpp>
#include <b/geometry/algorithms/area.hpp>
#include <b/geometry/algorithms/assign.hpp>
#include <b/geometry/algorithms/azimuth.hpp>
#include <b/geometry/algorithms/buffer.hpp>
#include <b/geometry/algorithms/centroid.hpp>
#include <b/geometry/algorithms/clear.hpp>
#include <b/geometry/algorithms/closest_points.hpp>
#include <b/geometry/algorithms/comparable_distance.hpp>
#include <b/geometry/algorithms/convert.hpp>
#include <b/geometry/algorithms/convex_hull.hpp>
#include <b/geometry/algorithms/correct.hpp>
#include <b/geometry/algorithms/covered_by.hpp>
#include <b/geometry/algorithms/crosses.hpp>
#include <b/geometry/algorithms/densify.hpp>
#include <b/geometry/algorithms/difference.hpp>
#include <b/geometry/algorithms/discrete_frechet_distance.hpp>
#include <b/geometry/algorithms/discrete_hausdorff_distance.hpp>
#include <b/geometry/algorithms/disjoint.hpp>
#include <b/geometry/algorithms/distance.hpp>
#include <b/geometry/algorithms/envelope.hpp>
#include <b/geometry/algorithms/equals.hpp>
#include <b/geometry/algorithms/expand.hpp>
#include <b/geometry/algorithms/for_each.hpp>
#include <b/geometry/algorithms/intersection.hpp>
#include <b/geometry/algorithms/intersects.hpp>
#include <b/geometry/algorithms/is_convex.hpp>
#include <b/geometry/algorithms/is_empty.hpp>
#include <b/geometry/algorithms/is_simple.hpp>
#include <b/geometry/algorithms/is_valid.hpp>
#include <b/geometry/algorithms/length.hpp>
#include <b/geometry/algorithms/line_interpolate.hpp>
#include <b/geometry/algorithms/make.hpp>
#include <b/geometry/algorithms/num_geometries.hpp>
#include <b/geometry/algorithms/num_interior_rings.hpp>
#include <b/geometry/algorithms/num_points.hpp>
#include <b/geometry/algorithms/num_segments.hpp>
#include <b/geometry/algorithms/overlaps.hpp>
#include <b/geometry/algorithms/perimeter.hpp>
#include <b/geometry/algorithms/point_on_surface.hpp>
#include <b/geometry/algorithms/relate.hpp>
#include <b/geometry/algorithms/relation.hpp>
#include <b/geometry/algorithms/remove_spikes.hpp>
#include <b/geometry/algorithms/reverse.hpp>
#include <b/geometry/algorithms/simplify.hpp>
#include <b/geometry/algorithms/sym_difference.hpp>
#include <b/geometry/algorithms/touches.hpp>
#include <b/geometry/algorithms/transform.hpp>
#include <b/geometry/algorithms/union.hpp>
#include <b/geometry/algorithms/unique.hpp>
#include <b/geometry/algorithms/within.hpp>

// check includes all concepts
#include <b/geometry/geometries/concepts/check.hpp>

#include <b/geometry/srs/srs.hpp>

#include <b/geometry/util/for_each_coordinate.hpp>
#include <b/geometry/util/math.hpp>
#include <b/geometry/util/select_coordinate_type.hpp>
#include <b/geometry/util/select_most_precise.hpp>

#include <b/geometry/views/box_view.hpp>
#include <b/geometry/views/closeable_view.hpp>
#include <b/geometry/views/identity_view.hpp>
#include <b/geometry/views/reversible_view.hpp>
#include <b/geometry/views/segment_view.hpp>

#include <b/geometry/io/io.hpp>
#include <b/geometry/io/dsv/write.hpp>
#include <b/geometry/io/svg/svg_mapper.hpp>
#include <b/geometry/io/svg/write.hpp>
#include <b/geometry/io/wkt/read.hpp>
#include <b/geometry/io/wkt/write.hpp>

#include <b/geometry/algorithms/is_convex.hpp>
#include <b/geometry/algorithms/point_on_surface.hpp>

#endif // BOOST_GEOMETRY_GEOMETRY_HPP
