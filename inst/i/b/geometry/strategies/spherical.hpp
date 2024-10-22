// Boost.Geometry

// Copyright (c) 2020-2021, Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_STRATEGIES_SPHERICAL_HPP
#define BOOST_GEOMETRY_STRATEGIES_SPHERICAL_HPP


#include <b/geometry/strategies/area/spherical.hpp>
#include <b/geometry/strategies/azimuth/spherical.hpp>
#include <b/geometry/strategies/buffer/spherical.hpp>
#include <b/geometry/strategies/convex_hull/spherical.hpp>
#include <b/geometry/strategies/distance/spherical.hpp>
#include <b/geometry/strategies/envelope/spherical.hpp>
#include <b/geometry/strategies/expand/spherical.hpp>
#include <b/geometry/strategies/io/spherical.hpp>
#include <b/geometry/strategies/index/spherical.hpp>
#include <b/geometry/strategies/is_convex/spherical.hpp>
#include <b/geometry/strategies/relate/spherical.hpp>
#include <b/geometry/strategies/simplify/spherical.hpp>


namespace boost { namespace geometry
{


namespace strategies
{


template
<
    typename RadiusTypeOrSphere = double,
    typename CalculationType = void
>
class spherical
    // derived from the umbrella strategy defining the most strategies
    : public strategies::index::detail::spherical<RadiusTypeOrSphere, CalculationType>
{
    using base_t = strategies::index::detail::spherical<RadiusTypeOrSphere, CalculationType>;

public:
    spherical() = default;

    template <typename RadiusOrSphere>
    explicit spherical(RadiusOrSphere const& radius_or_sphere)
        : base_t(radius_or_sphere)
    {}

    static auto azimuth()
    {
        return strategy::azimuth::spherical<CalculationType>();
    }

    static auto point_order()
    {
        return strategy::point_order::spherical<CalculationType>();
    }
};


} // namespace strategies


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_STRATEGIES_SPHERICAL_HPP
