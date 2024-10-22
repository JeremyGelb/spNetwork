/*=============================================================================
    Copyright (c) 2005-2012 Joel de Guzman
    Copyright (c) 2005-2006 Dan Marsden

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_DEQUE_AT_IMPL_09122006_2017)
#define BOOST_FUSION_DEQUE_AT_IMPL_09122006_2017

#include <b/fusion/support/config.hpp>
#include <b/fusion/container/deque/detail/keyed_element.hpp>

#include <b/mpl/eval_if.hpp>
#include <b/mpl/equal_to.hpp>
#include <b/mpl/assert.hpp>
#include <b/mpl/identity.hpp>

#include <b/type_traits/is_const.hpp>
#include <b/type_traits/add_const.hpp>
#include <b/type_traits/add_reference.hpp>

namespace boost { namespace fusion
{
    struct deque_tag;

    namespace extension
    {
        template<typename T>
        struct at_impl;

        template<>
        struct at_impl<deque_tag>
        {
            template<typename Sequence, typename N>
            struct apply
            {
                typedef typename Sequence::next_up next_up;
                typedef typename Sequence::next_down next_down;
                BOOST_MPL_ASSERT_RELATION(next_down::value, !=, next_up::value);

                static int const offset = next_down::value + 1;
                typedef mpl::int_<(N::value + offset)> adjusted_index;
                typedef typename
                    detail::keyed_element_value_at<Sequence, adjusted_index>::type
                element_type;

                typedef typename
                    add_reference<
                      typename mpl::eval_if<
                      is_const<Sequence>,
                      add_const<element_type>,
                      mpl::identity<element_type> >::type
                    >::type
                type;

                BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
                static type call(Sequence& seq)
                {
                    return seq.get(adjusted_index());
                }
            };
        };
    }
}}

#endif
