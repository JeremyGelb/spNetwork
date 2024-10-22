/*=============================================================================
    Copyright (c) 2001-2011 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_VECTOR10_05042005_0257)
#define FUSION_VECTOR10_05042005_0257

#include <b/fusion/support/config.hpp>
#include <b/fusion/container/vector/detail/cpp03/vector10_fwd.hpp>
#include <b/fusion/support/sequence_base.hpp>
#include <b/fusion/support/is_sequence.hpp>
#include <b/fusion/support/detail/access.hpp>
#include <b/fusion/iterator/next.hpp>
#include <b/fusion/iterator/deref.hpp>
#include <b/fusion/sequence/intrinsic/begin.hpp>
#include <b/fusion/container/vector/detail/at_impl.hpp>
#include <b/fusion/container/vector/detail/value_at_impl.hpp>
#include <b/fusion/container/vector/detail/begin_impl.hpp>
#include <b/fusion/container/vector/detail/end_impl.hpp>

#include <b/mpl/void.hpp>
#include <b/mpl/int.hpp>
#include <b/mpl/bool.hpp>
#include <b/mpl/at.hpp>
#include <b/mpl/vector/vector10.hpp>
#include <b/type_traits/is_convertible.hpp>
#include <b/utility/enable_if.hpp>

#include <b/preprocessor/dec.hpp>
#include <b/preprocessor/iteration/iterate.hpp>
#include <b/preprocessor/repetition/enum.hpp>
#include <b/preprocessor/repetition/enum_shifted.hpp>
#include <b/preprocessor/repetition/enum_params.hpp>
#include <b/preprocessor/repetition/enum_binary_params.hpp>
#include <b/preprocessor/repetition/repeat_from_to.hpp>

namespace boost { namespace fusion
{
    struct vector_tag;
    struct fusion_sequence_tag;
    struct random_access_traversal_tag;

    template <typename Dummy>
    struct vector0 : sequence_base<vector0<Dummy> >
    {
        typedef mpl::vector0<> types;
        typedef vector_tag fusion_tag;
        typedef fusion_sequence_tag tag; // this gets picked up by MPL
        typedef mpl::false_ is_view;
        typedef random_access_traversal_tag category;
        typedef mpl::int_<0> size;

        BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
        vector0() BOOST_NOEXCEPT {}

        template<typename Sequence>
        BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED
        vector0(Sequence const& /*seq*/) BOOST_NOEXCEPT
        {}
    };
}}

#if !defined(BOOST_FUSION_DONT_USE_PREPROCESSED_FILES)
#include <b/fusion/container/vector/detail/cpp03/preprocessed/vector10.hpp>
#else
#if defined(__WAVE__) && defined(BOOST_FUSION_CREATE_PREPROCESSED_FILES)
#pragma wave option(preserve: 2, line: 0, output: "preprocessed/vector10.hpp")
#endif

/*=============================================================================
    Copyright (c) 2001-2011 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    This is an auto-generated file. Do not edit!
==============================================================================*/

#if defined(__WAVE__) && defined(BOOST_FUSION_CREATE_PREPROCESSED_FILES)
#pragma wave option(preserve: 1)
#endif

namespace boost { namespace fusion
{
    struct vector_tag;
    struct fusion_sequence_tag;
    struct random_access_traversal_tag;

#define FUSION_HASH #
// expand vector1 to vector10
#define BOOST_PP_FILENAME_1 <b/fusion/container/vector/detail/cpp03/vector_n.hpp>
#define BOOST_PP_ITERATION_LIMITS (1, 10)
#include BOOST_PP_ITERATE()
#undef FUSION_HASH
}}

#if defined(__WAVE__) && defined(BOOST_FUSION_CREATE_PREPROCESSED_FILES)
#pragma wave option(output: null)
#endif

#endif // BOOST_FUSION_DONT_USE_PREPROCESSED_FILES

#endif
