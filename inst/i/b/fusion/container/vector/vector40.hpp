/*=============================================================================
    Copyright (c) 2014 Kohei Takahashi

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef FUSION_VECTOR40_11052014_2316
#define FUSION_VECTOR40_11052014_2316

#include <b/config.hpp>
#include <b/fusion/support/config.hpp>
#include <b/fusion/container/vector/detail/config.hpp>

///////////////////////////////////////////////////////////////////////////////
// Without variadics, we will use the PP version
///////////////////////////////////////////////////////////////////////////////
#if !defined(BOOST_FUSION_HAS_VARIADIC_VECTOR)
#   include <b/fusion/container/vector/detail/cpp03/vector40.hpp>
#else

///////////////////////////////////////////////////////////////////////////////
// C++11 interface
///////////////////////////////////////////////////////////////////////////////
#include <b/fusion/container/vector/vector_fwd.hpp>
#include <b/fusion/container/vector/vector.hpp>

#endif
#endif

