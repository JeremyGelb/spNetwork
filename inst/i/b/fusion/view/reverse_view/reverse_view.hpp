/*=============================================================================
    Copyright (c) 2001-2011 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_REVERSE_VIEW_07202005_0836)
#define FUSION_REVERSE_VIEW_07202005_0836

#include <b/fusion/support/config.hpp>
#include <b/fusion/support/detail/access.hpp>
#include <b/fusion/support/is_view.hpp>
#include <b/fusion/support/category_of.hpp>
#include <b/fusion/view/reverse_view/reverse_view_iterator.hpp>
#include <b/fusion/view/reverse_view/detail/begin_impl.hpp>
#include <b/fusion/view/reverse_view/detail/end_impl.hpp>
#include <b/fusion/view/reverse_view/detail/at_impl.hpp>
#include <b/fusion/view/reverse_view/detail/value_at_impl.hpp>
#include <b/fusion/support/sequence_base.hpp>
#include <b/fusion/sequence/intrinsic/begin.hpp>
#include <b/fusion/sequence/intrinsic/end.hpp>
#include <b/fusion/sequence/intrinsic/size.hpp>
#include <b/type_traits/is_base_of.hpp>
#include <b/static_assert.hpp>
#include <b/mpl/bool.hpp>
#include <b/mpl/eval_if.hpp>
#include <b/mpl/inherit.hpp>
#include <b/mpl/identity.hpp>

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#endif

namespace boost { namespace fusion
{
    struct reverse_view_tag;
    struct fusion_sequence_tag;

    template <typename Sequence>
    struct reverse_view : sequence_base<reverse_view<Sequence> >
    {
        typedef reverse_view_tag fusion_tag;
        typedef fusion_sequence_tag tag; // this gets picked up by MPL
        typedef mpl::true_ is_view;

        typedef Sequence seq_type;
        typedef typename traits::category_of<Sequence>::type category;
        typedef typename result_of::begin<Sequence>::type first_type;
        typedef typename result_of::end<Sequence>::type last_type;
        typedef typename result_of::size<Sequence>::type size;

        BOOST_STATIC_ASSERT((
            is_base_of<
                bidirectional_traversal_tag
              , typename traits::category_of<first_type>::type>::value));

        BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
        reverse_view(Sequence& in_seq)
            : seq(in_seq)
        {}

        BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
        first_type first() const { return fusion::begin(seq); }
        BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
        last_type last() const { return fusion::end(seq); }
        typename mpl::if_<traits::is_view<Sequence>, Sequence, Sequence&>::type seq;
    };
}}

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#endif


