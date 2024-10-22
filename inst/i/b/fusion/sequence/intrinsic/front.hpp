/*=============================================================================
    Copyright (c) 2001-2011 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_FRONT_09162005_0343)
#define FUSION_FRONT_09162005_0343

#include <b/fusion/support/config.hpp>
#include <b/fusion/sequence/intrinsic_fwd.hpp>
#include <b/fusion/sequence/intrinsic/begin.hpp>
#include <b/fusion/iterator/deref.hpp>
#include <b/mpl/bool.hpp>

namespace boost { namespace fusion
{
    struct fusion_sequence_tag;

    namespace result_of
    {
        template <typename Sequence>
        struct front
            : result_of::deref<typename result_of::begin<Sequence>::type>
        {};
    }
    
    template <typename Sequence>
    BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
    inline typename result_of::front<Sequence>::type
    front(Sequence& seq)
    {
        return *fusion::begin(seq);
    }

    template <typename Sequence>
    BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
    inline typename result_of::front<Sequence const>::type
    front(Sequence const& seq)
    {
        return *fusion::begin(seq);
    }
}}

#endif
