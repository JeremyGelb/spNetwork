// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2011-2015 Akira Takahashi
// Copyright (c) 2011-2015 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2015-2020.
// Modifications copyright (c) 2015-2020, Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_GEOMETRIES_ADAPTED_FUSION_HPP
#define BOOST_GEOMETRY_GEOMETRIES_ADAPTED_FUSION_HPP


#include <cstddef>
#include <type_traits>

#include <b/fusion/include/is_sequence.hpp>
#include <b/fusion/include/size.hpp>
#include <b/fusion/include/tag_of.hpp>
#include <b/fusion/include/front.hpp>
#include <b/fusion/include/at.hpp>
#include <b/fusion/mpl.hpp>

#include <b/mpl/and.hpp>
#include <b/mpl/count_if.hpp>
#include <b/mpl/front.hpp>
#include <b/mpl/placeholders.hpp>
#include <b/mpl/pop_front.hpp>
#include <b/mpl/size.hpp>

#include <b/type_traits/is_same.hpp>

#include <b/geometry/core/access.hpp>
#include <b/geometry/core/coordinate_dimension.hpp>
#include <b/geometry/core/coordinate_system.hpp>
#include <b/geometry/core/coordinate_type.hpp>
#include <b/geometry/core/point_type.hpp>
#include <b/geometry/core/tags.hpp>


namespace boost { namespace geometry
{

namespace fusion_adapt_detail
{

template <class Sequence>
struct all_same :
    std::integral_constant
        <
            bool,
            boost::mpl::count_if<
                Sequence,
                boost::is_same<
                    typename boost::mpl::front<Sequence>::type,
                    boost::mpl::_
                >
            >::value == boost::mpl::size<Sequence>::value
        >
{};

template <class Sequence>
struct is_coordinate_size
    : std::integral_constant
        <
            bool,
            (boost::fusion::result_of::size<Sequence>::value == 2
          || boost::fusion::result_of::size<Sequence>::value == 3)
        >
{};

template
<
    typename Sequence,
    bool IsSequence = boost::fusion::traits::is_sequence<Sequence>::value
>
struct is_fusion_sequence
    : std::integral_constant
        <
            bool,
            (fusion_adapt_detail::is_coordinate_size<Sequence>::value
          && fusion_adapt_detail::all_same<Sequence>::value)
        >
{};

template<typename Sequence>
struct is_fusion_sequence<Sequence, false>
    : std::false_type
{};


} // namespace fusion_adapt_detail


#ifndef DOXYGEN_NO_TRAITS_SPECIALIZATIONS
namespace traits
{

// Boost Fusion Sequence, 2D or 3D
template <typename Sequence>
struct coordinate_type
    <
        Sequence,
        std::enable_if_t
            <
                fusion_adapt_detail::is_fusion_sequence<Sequence>::value
            >
    >
{
    typedef typename boost::mpl::front<Sequence>::type type;
};


template <typename Sequence>
struct dimension
    <
        Sequence,
        std::enable_if_t
            <
                fusion_adapt_detail::is_fusion_sequence<Sequence>::value
            >
    > : boost::mpl::size<Sequence>
{};


template <typename Sequence, std::size_t Dimension>
struct access
    <
        Sequence,
        Dimension,
        std::enable_if_t
            <
                fusion_adapt_detail::is_fusion_sequence<Sequence>::value
            >
    >
{
    typedef typename coordinate_type<Sequence>::type ctype;

    static inline ctype get(Sequence const& point)
    {
        return boost::fusion::at_c<Dimension>(point);
    }

    template <class CoordinateType>
    static inline void set(Sequence& point, CoordinateType const& value)
    {
        boost::fusion::at_c<Dimension>(point) = value;
    }
};


template <typename Sequence>
struct tag
    <
        Sequence,
        std::enable_if_t
            <
                fusion_adapt_detail::is_fusion_sequence<Sequence>::value
            >
    >
{
    typedef point_tag type;
};


} // namespace traits

#endif // DOXYGEN_NO_TRAITS_SPECIALIZATIONS


}} // namespace boost::geometry


// Convenience registration macro to bind a Fusion sequence to a CS
#define BOOST_GEOMETRY_REGISTER_BOOST_FUSION_CS(CoordinateSystem) \
    namespace boost { namespace geometry { namespace traits { \
    template <typename Sequence> \
    struct coordinate_system \
               < \
                   Sequence, \
                   std::enable_if_t \
                       < \
                           fusion_adapt_detail::is_fusion_sequence<Sequence>::value \
                       > \
               > \
    { typedef CoordinateSystem type; }; \
    }}}


#endif // BOOST_GEOMETRY_GEOMETRIES_ADAPTED_FUSION_HPP
