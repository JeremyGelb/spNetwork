/*=============================================================================
    Copyright (c) 2001-2011 Joel de Guzman
    Copyright (c) 2006 Dan Marsden

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_EQUAL_TO_IMPL_20060128_1423)
#define FUSION_EQUAL_TO_IMPL_20060128_1423

#include <b/fusion/support/config.hpp>
#include <b/fusion/mpl.hpp>

#include <b/mpl/lambda.hpp>
#include <b/mpl/and.hpp>
#include <b/mpl/transform_view.hpp>
#include <b/mpl/zip_view.hpp>
#include <b/mpl/vector.hpp>
#include <b/mpl/unpack_args.hpp>
#include <b/mpl/find_if.hpp>
#include <b/mpl/end.hpp>
#include <b/mpl/bool.hpp>
#include <b/mpl/equal_to.hpp>

#include <b/type_traits/is_same.hpp>

#include <b/fusion/iterator/equal_to.hpp>

namespace boost { namespace fusion {

    struct zip_view_iterator_tag;

    namespace detail
    {
        template<typename It1, typename It2>
        struct zip_iterators_equal
        {
            typedef mpl::zip_view<mpl::vector2<typename It1::iterators, typename It2::iterators> > zipped;
            typedef mpl::transform_view<zipped, mpl::unpack_args<result_of::equal_to<mpl::_,mpl::_> > > transformed;

            typedef typename mpl::find_if<transformed, mpl::equal_to<mpl::_, mpl::false_> >::type found;

            typedef typename is_same<typename mpl::end<transformed>::type, found>::type type;
        };
    }

    namespace extension
    {
        template<typename Tag>
        struct equal_to_impl;

        template<>
        struct equal_to_impl<zip_view_iterator_tag>
        {
            template<typename It1, typename It2>
            struct apply
                : detail::zip_iterators_equal<It1, It2>::type
            {};
        };
    }
}}

#endif
