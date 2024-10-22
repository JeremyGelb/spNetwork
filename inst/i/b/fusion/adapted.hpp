/*=============================================================================
    Copyright (c) 2001-2011 Joel de Guzman
    Copyright (c) 2005-2006 Dan Marsden

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_ADAPTED_30122005_1420)
#define BOOST_FUSION_ADAPTED_30122005_1420

#include <b/fusion/support/config.hpp>
#include <b/fusion/adapted/adt.hpp>
#include <b/fusion/adapted/array.hpp>
#include <b/fusion/adapted/boost_array.hpp>
#include <b/fusion/adapted/boost_tuple.hpp>
#include <b/fusion/adapted/mpl.hpp>
#include <b/fusion/adapted/std_pair.hpp>
#include <b/fusion/adapted/struct.hpp>

// The std_tuple_iterator adaptor only supports implementations
// using variadic templates
#if !defined(BOOST_NO_CXX11_VARIADIC_TEMPLATES)
#include <b/fusion/adapted/std_tuple.hpp>
#endif

#endif
