// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2017 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2014-2020.
// Modifications copyright (c) 2014-2020 Oracle and/or its affiliates.

// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_STRATEGIES_HPP
#define BOOST_GEOMETRY_STRATEGIES_STRATEGIES_HPP


#include <b/geometry/strategies/tags.hpp>

#include <b/geometry/strategies/azimuth.hpp>
#include <b/geometry/strategies/buffer.hpp>
#include <b/geometry/strategies/centroid.hpp>
#include <b/geometry/strategies/compare.hpp>
#include <b/geometry/strategies/convex_hull.hpp>
#include <b/geometry/strategies/covered_by.hpp>
#include <b/geometry/strategies/densify.hpp>
#include <b/geometry/strategies/disjoint.hpp>
#include <b/geometry/strategies/distance.hpp>
#include <b/geometry/strategies/intersection.hpp>
#include <b/geometry/strategies/intersection_strategies.hpp> // for backward compatibility
#include <b/geometry/strategies/side.hpp>
#include <b/geometry/strategies/transform.hpp>
#include <b/geometry/strategies/within.hpp>

#include <b/geometry/strategies/cartesian/azimuth.hpp>
#include <b/geometry/strategies/cartesian/box_in_box.hpp>
#include <b/geometry/strategies/cartesian/buffer_end_flat.hpp>
#include <b/geometry/strategies/cartesian/buffer_end_round.hpp>
#include <b/geometry/strategies/cartesian/buffer_join_miter.hpp>
#include <b/geometry/strategies/cartesian/buffer_join_round.hpp>
#include <b/geometry/strategies/cartesian/buffer_join_round_by_divide.hpp>
#include <b/geometry/strategies/cartesian/buffer_point_circle.hpp>
#include <b/geometry/strategies/cartesian/buffer_point_square.hpp>
#include <b/geometry/strategies/cartesian/buffer_side_straight.hpp>
#include <b/geometry/strategies/cartesian/centroid_average.hpp>
#include <b/geometry/strategies/cartesian/centroid_bashein_detmer.hpp>
#include <b/geometry/strategies/cartesian/centroid_weighted_length.hpp>
#include <b/geometry/strategies/cartesian/densify.hpp>
#include <b/geometry/strategies/cartesian/disjoint_segment_box.hpp>
#include <b/geometry/strategies/cartesian/distance_pythagoras.hpp>
#include <b/geometry/strategies/cartesian/distance_pythagoras_point_box.hpp>
#include <b/geometry/strategies/cartesian/distance_pythagoras_box_box.hpp>
#include <b/geometry/strategies/cartesian/distance_projected_point.hpp>
#include <b/geometry/strategies/cartesian/distance_projected_point_ax.hpp>
#include <b/geometry/strategies/cartesian/distance_segment_box.hpp>
#include <b/geometry/strategies/cartesian/intersection.hpp>
#include <b/geometry/strategies/cartesian/point_in_box.hpp>
#include <b/geometry/strategies/cartesian/point_in_point.hpp>
#include <b/geometry/strategies/cartesian/point_in_poly_franklin.hpp>
#include <b/geometry/strategies/cartesian/point_in_poly_crossings_multiply.hpp>
#include <b/geometry/strategies/cartesian/point_in_poly_winding.hpp>
#include <b/geometry/strategies/cartesian/line_interpolate.hpp>

#include <b/geometry/strategies/spherical/azimuth.hpp>
#include <b/geometry/strategies/spherical/densify.hpp>
#include <b/geometry/strategies/spherical/disjoint_segment_box.hpp>
#include <b/geometry/strategies/spherical/distance_haversine.hpp>
#include <b/geometry/strategies/spherical/distance_cross_track.hpp>
#include <b/geometry/strategies/spherical/distance_cross_track_box_box.hpp>
#include <b/geometry/strategies/spherical/distance_cross_track_point_box.hpp>
#include <b/geometry/strategies/spherical/distance_segment_box.hpp>
#include <b/geometry/strategies/spherical/compare.hpp>
#include <b/geometry/strategies/spherical/intersection.hpp>
#include <b/geometry/strategies/spherical/point_in_point.hpp>
#include <b/geometry/strategies/spherical/point_in_poly_winding.hpp>
#include <b/geometry/strategies/spherical/line_interpolate.hpp>
#include <b/geometry/strategies/spherical/ssf.hpp>

#include <b/geometry/strategies/geographic/azimuth.hpp>
#include <b/geometry/strategies/geographic/buffer_end_round.hpp>
#include <b/geometry/strategies/geographic/buffer_join_miter.hpp>
#include <b/geometry/strategies/geographic/buffer_join_round.hpp>
#include <b/geometry/strategies/geographic/buffer_point_circle.hpp>
#include <b/geometry/strategies/geographic/buffer_side_straight.hpp>
#include <b/geometry/strategies/geographic/densify.hpp>
#include <b/geometry/strategies/geographic/disjoint_segment_box.hpp>
#include <b/geometry/strategies/geographic/distance.hpp>
#include <b/geometry/strategies/geographic/distance_andoyer.hpp>
#include <b/geometry/strategies/geographic/distance_cross_track.hpp>
#include <b/geometry/strategies/geographic/distance_cross_track_box_box.hpp>
#include <b/geometry/strategies/geographic/distance_cross_track_point_box.hpp>
#include <b/geometry/strategies/geographic/distance_segment_box.hpp>
#include <b/geometry/strategies/geographic/distance_thomas.hpp>
#include <b/geometry/strategies/geographic/distance_vincenty.hpp>
#include <b/geometry/strategies/geographic/distance_karney.hpp>
#include <b/geometry/strategies/geographic/intersection.hpp>
//#include <b/geometry/strategies/geographic/intersection_elliptic.hpp>
#include <b/geometry/strategies/geographic/point_in_poly_winding.hpp>
#include <b/geometry/strategies/geographic/line_interpolate.hpp>
#include <b/geometry/strategies/geographic/side.hpp>
#include <b/geometry/strategies/geographic/side_andoyer.hpp>
#include <b/geometry/strategies/geographic/side_thomas.hpp>
#include <b/geometry/strategies/geographic/side_vincenty.hpp>

#include <b/geometry/strategies/agnostic/buffer_distance_symmetric.hpp>
#include <b/geometry/strategies/agnostic/buffer_distance_asymmetric.hpp>
#include <b/geometry/strategies/agnostic/point_in_box_by_side.hpp>
#include <b/geometry/strategies/agnostic/point_in_point.hpp>
#include <b/geometry/strategies/agnostic/point_in_poly_winding.hpp>
#include <b/geometry/strategies/agnostic/simplify_douglas_peucker.hpp>

#include <b/geometry/strategies/strategy_transform.hpp>

#include <b/geometry/strategies/transform/matrix_transformers.hpp>
#include <b/geometry/strategies/transform/map_transformer.hpp>
#include <b/geometry/strategies/transform/inverse_transformer.hpp>

// TEMP

#include <b/geometry/strategy/area.hpp>
#include <b/geometry/strategy/envelope.hpp>
#include <b/geometry/strategy/expand.hpp>
#include <b/geometry/strategy/relate.hpp>

#include <b/geometry/strategy/cartesian/area.hpp>
#include <b/geometry/strategy/cartesian/envelope.hpp>
#include <b/geometry/strategy/cartesian/envelope_box.hpp>
#include <b/geometry/strategy/cartesian/envelope_multipoint.hpp>
#include <b/geometry/strategy/cartesian/envelope_point.hpp>
#include <b/geometry/strategy/cartesian/envelope_segment.hpp>
#include <b/geometry/strategy/cartesian/expand_box.hpp>
#include <b/geometry/strategy/cartesian/expand_point.hpp>
#include <b/geometry/strategy/cartesian/expand_segment.hpp>
#include <b/geometry/strategy/cartesian/side_by_triangle.hpp>
#include <b/geometry/strategy/cartesian/side_robust.hpp>

#include <b/geometry/strategy/geographic/area.hpp>
#include <b/geometry/strategy/geographic/envelope.hpp>
#include <b/geometry/strategy/geographic/envelope_segment.hpp>
#include <b/geometry/strategy/geographic/expand_segment.hpp>

#include <b/geometry/strategy/spherical/area.hpp>
#include <b/geometry/strategy/spherical/envelope.hpp>
#include <b/geometry/strategy/spherical/envelope_box.hpp>
#include <b/geometry/strategy/spherical/envelope_multipoint.hpp>
#include <b/geometry/strategy/spherical/envelope_point.hpp>
#include <b/geometry/strategy/spherical/envelope_segment.hpp>
#include <b/geometry/strategy/spherical/expand_box.hpp>
#include <b/geometry/strategy/spherical/expand_point.hpp>
#include <b/geometry/strategy/spherical/expand_segment.hpp>

#include <b/geometry/strategies/cartesian.hpp>
#include <b/geometry/strategies/geographic.hpp>
#include <b/geometry/strategies/spherical.hpp>


#endif // BOOST_GEOMETRY_STRATEGIES_STRATEGIES_HPP
