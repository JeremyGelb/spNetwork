// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2014-2021, Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_IS_VALID_IMPLEMENTATION_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_IS_VALID_IMPLEMENTATION_HPP

#include <b/geometry/algorithms/detail/is_valid/pointlike.hpp>
#include <b/geometry/algorithms/detail/is_valid/linear.hpp>
#include <b/geometry/algorithms/detail/is_valid/polygon.hpp>
#include <b/geometry/algorithms/detail/is_valid/multipolygon.hpp>
#include <b/geometry/algorithms/detail/is_valid/ring.hpp>
#include <b/geometry/algorithms/detail/is_valid/segment.hpp>
#include <b/geometry/algorithms/detail/is_valid/box.hpp>

#include <b/geometry/strategies/relate/cartesian.hpp>
#include <b/geometry/strategies/relate/geographic.hpp>
#include <b/geometry/strategies/relate/spherical.hpp>

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_IS_VALID_IMPLEMENTATION_HPP
