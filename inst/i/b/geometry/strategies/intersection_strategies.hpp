// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2016-2020.
// Modifications copyright (c) 2016-2020, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_INTERSECTION_HPP
#define BOOST_GEOMETRY_STRATEGIES_INTERSECTION_HPP


//#include <type_traits>

#include <b/geometry/policies/relate/intersection_policy.hpp>
#include <b/geometry/policies/robustness/segment_ratio_type.hpp>

#include <b/geometry/strategies/intersection.hpp>
#include <b/geometry/strategies/intersection_result.hpp>
#include <b/geometry/strategies/side.hpp>
#include <b/geometry/strategies/cartesian/intersection.hpp>
#include <b/geometry/strategies/spherical/intersection.hpp>
#include <b/geometry/strategies/spherical/ssf.hpp>


namespace boost { namespace geometry
{


/*!
\brief "compound strategy", containing a segment-intersection-strategy
       and a side-strategy
 */
template
<
    typename Tag,
    typename Geometry1,
    typename Geometry2,
    typename IntersectionPoint,
    typename RobustPolicy,
    typename CalculationType = void
>
struct intersection_strategies
{
private :
    // for development BOOST_STATIC_ASSERT((! std::is_same<RobustPolicy, void>::type::value));

    typedef segment_intersection_points
    <
        IntersectionPoint,
        typename detail::segment_ratio_type
        <
            IntersectionPoint, RobustPolicy
        >::type
    > ip_type;

public:
    typedef policies::relate::segments_intersection_policy
        <
            ip_type
        > intersection_policy_type;

    typedef typename strategy::intersection::services::default_strategy
            <
                Tag,
                CalculationType
            >::type segment_intersection_strategy_type;

    typedef typename strategy::side::services::default_strategy
        <
            Tag,
            CalculationType
        >::type side_strategy_type;

    typedef RobustPolicy rescale_policy_type;
};


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_STRATEGIES_INTERSECTION_HPP
