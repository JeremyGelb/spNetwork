// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2015 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2015 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2015 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2015-2021.
// Modifications copyright (c) 2015-2021, Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_ENVELOPE_IMPLEMENTATION_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_ENVELOPE_IMPLEMENTATION_HPP

#include <b/geometry/core/exterior_ring.hpp>
#include <b/geometry/core/interior_rings.hpp>
#include <b/geometry/core/tags.hpp>

#include <b/geometry/algorithms/is_empty.hpp>

#include <b/geometry/algorithms/detail/envelope/areal.hpp>
#include <b/geometry/algorithms/detail/envelope/box.hpp>
#include <b/geometry/algorithms/detail/envelope/geometry_collection.hpp>
#include <b/geometry/algorithms/detail/envelope/linear.hpp>
#include <b/geometry/algorithms/detail/envelope/multipoint.hpp>
#include <b/geometry/algorithms/detail/envelope/point.hpp>
#include <b/geometry/algorithms/detail/envelope/range.hpp>
#include <b/geometry/algorithms/detail/envelope/segment.hpp>

#include <b/geometry/algorithms/dispatch/envelope.hpp>

#include <b/geometry/strategies/envelope/cartesian.hpp>
#include <b/geometry/strategies/envelope/geographic.hpp>
#include <b/geometry/strategies/envelope/spherical.hpp>

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_ENVELOPE_IMPLEMENTATION_HPP
