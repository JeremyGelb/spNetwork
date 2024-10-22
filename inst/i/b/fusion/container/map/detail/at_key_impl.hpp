/*=============================================================================
    Copyright (c) 2001-2013 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_MAP_DETAIL_AT_KEY_IMPL_02042013_0821)
#define BOOST_FUSION_MAP_DETAIL_AT_KEY_IMPL_02042013_0821

#include <b/fusion/support/config.hpp>
#include <b/fusion/support/detail/access.hpp>
#include <b/type_traits/is_const.hpp>
#include <b/mpl/at.hpp>
#include <b/mpl/identity.hpp>
#include <b/utility/declval.hpp>

namespace boost { namespace fusion
{
    struct map_tag;

    namespace extension
    {
        template <typename Tag>
        struct at_key_impl;

        template <>
        struct at_key_impl<map_tag>
        {
            template <typename Sequence, typename Key>
            struct apply
            {
                typedef
                    decltype(boost::declval<Sequence>().get(mpl::identity<Key>()))
                type;

                BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
                static type
                call(Sequence& m)
                {
                    return m.get(mpl::identity<Key>());
                }
            };

            template <typename Sequence, typename Key>
            struct apply<Sequence const, Key>
            {
                typedef
                    decltype(boost::declval<Sequence const>().get(mpl::identity<Key>()))
                type;

                BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
                static type
                call(Sequence const& m)
                {
                    return m.get(mpl::identity<Key>());
                }
            };
        };
    }
}}

#endif
