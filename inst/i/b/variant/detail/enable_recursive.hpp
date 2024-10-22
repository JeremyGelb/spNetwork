//-----------------------------------------------------------------------------
// boost variant/detail/enable_recursive.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2003
// Eric Friedman
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_VARIANT_DETAIL_ENABLE_RECURSIVE_HPP
#define BOOST_VARIANT_DETAIL_ENABLE_RECURSIVE_HPP

#include <b/variant/detail/enable_recursive_fwd.hpp>
#include <b/variant/variant_fwd.hpp>

#if !defined(BOOST_VARIANT_NO_FULL_RECURSIVE_VARIANT_SUPPORT)
#   include <b/mpl/apply.hpp>
#   include <b/mpl/eval_if.hpp>
#   include <b/mpl/lambda.hpp>
#endif

#include <b/variant/detail/substitute.hpp>
#include <b/mpl/aux_/config/ctps.hpp>
#include <b/mpl/bool_fwd.hpp>
#include <b/mpl/if.hpp>
#include <b/mpl/or.hpp>
#include <b/type_traits/is_pointer.hpp>
#include <b/type_traits/is_reference.hpp>
#include <b/type_traits/is_same.hpp>

#include <b/variant/recursive_wrapper.hpp>

namespace boost {
namespace detail { namespace variant {

#   define BOOST_VARIANT_AUX_ENABLE_RECURSIVE_IMPL(T,Dest,Source) \
    substitute< T , Dest , Source > \
    /**/

///////////////////////////////////////////////////////////////////////////////
// (detail) metafunction enable_recursive
//
// See b/variant/detail/enable_recursive_fwd.hpp for more information.
//


template <typename T, typename RecursiveVariant, typename NoWrapper>
struct enable_recursive
    : BOOST_VARIANT_AUX_ENABLE_RECURSIVE_IMPL(
          T, RecursiveVariant, ::boost::recursive_variant_
        )
{
};

template <typename T, typename RecursiveVariant>
struct enable_recursive< T,RecursiveVariant,mpl::false_ >
{
private: // helpers, for metafunction result (below)

    typedef typename BOOST_VARIANT_AUX_ENABLE_RECURSIVE_IMPL(
          T, RecursiveVariant, ::boost::recursive_variant_
        )::type t_;

public: // metafunction result

    // [Wrap with recursive_wrapper only if rebind really changed something:]
    typedef typename mpl::if_<
          mpl::or_<
              is_same< t_,T >
            , is_reference<t_>
            , is_pointer<t_>
            >
        , t_
        , boost::recursive_wrapper<t_>
        >::type type;

};


///////////////////////////////////////////////////////////////////////////////
// (detail) metafunction class quoted_enable_recursive
//
// Same behavior as enable_recursive metafunction (see above).
//
template <typename RecursiveVariant, typename NoWrapper>
struct quoted_enable_recursive
{
    template <typename T>
    struct apply
        : enable_recursive<T, RecursiveVariant, NoWrapper>
    {
    };
};

}} // namespace detail::variant
} // namespace boost

#endif // BOOST_VARIANT_DETAIL_ENABLE_RECURSIVE_HPP
