// Boost.Geometry

// Copyright (c) 2020-2021, Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_STRATEGIES_GEOGRAPHIC_HPP
#define BOOST_GEOMETRY_STRATEGIES_GEOGRAPHIC_HPP


#include <b/geometry/strategies/area/geographic.hpp>
#include <b/geometry/strategies/azimuth/geographic.hpp>
#include <b/geometry/strategies/buffer/geographic.hpp>
#include <b/geometry/strategies/convex_hull/geographic.hpp>
#include <b/geometry/strategies/distance/geographic.hpp>
#include <b/geometry/strategies/envelope/geographic.hpp>
#include <b/geometry/strategies/expand/geographic.hpp>
#include <b/geometry/strategies/io/geographic.hpp>
#include <b/geometry/strategies/index/geographic.hpp>
#include <b/geometry/strategies/is_convex/geographic.hpp>
#include <b/geometry/strategies/relate/geographic.hpp>
#include <b/geometry/strategies/simplify/geographic.hpp>


namespace boost { namespace geometry
{


namespace strategies
{


template
<
    typename FormulaPolicy = strategy::andoyer,
    typename Spheroid = srs::spheroid<double>,
    typename CalculationType = void
>
class geographic
    // derived from the umbrella strategy defining the most strategies
    : public index::geographic<FormulaPolicy, Spheroid, CalculationType>
{
    using base_t = index::geographic<FormulaPolicy, Spheroid, CalculationType>;

public:
    geographic() = default;

    explicit geographic(Spheroid const& spheroid)
        : base_t(spheroid)
    {}

    auto azimuth() const
    {
        return strategy::azimuth::geographic
            <
                FormulaPolicy, Spheroid, CalculationType
            >(base_t::m_spheroid);
    }

    auto point_order() const
    {
        return strategy::point_order::geographic
            <
                FormulaPolicy, Spheroid, CalculationType
            >(base_t::m_spheroid);
    }
};


} // namespace strategies


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_STRATEGIES_GEOGRAPHIC_HPP
