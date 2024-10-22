
#if !defined(BOOST_PP_IS_ITERATING)

///// header body

#ifndef BOOST_MPL_AUX_TEMPLATE_ARITY_HPP_INCLUDED
#define BOOST_MPL_AUX_TEMPLATE_ARITY_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2001-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id$
// $Date$
// $Revision$

#include <b/mpl/aux_/config/ttp.hpp>
#include <b/mpl/aux_/config/lambda.hpp>

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include <b/mpl/aux_/template_arity_fwd.hpp>
#   include <b/mpl/int.hpp>
#   if !defined(BOOST_MPL_CFG_NO_FULL_LAMBDA_SUPPORT)
#   if defined(BOOST_MPL_CFG_EXTENDED_TEMPLATE_PARAMETERS_MATCHING)
#       include <b/mpl/aux_/type_wrapper.hpp>
#   endif
#   else
#       include <b/mpl/aux_/has_rebind.hpp>
#   endif
#endif

#include <b/mpl/aux_/config/static_constant.hpp>
#include <b/mpl/aux_/config/use_preprocessed.hpp>

#if !defined(BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER template_arity.hpp
#   include <b/mpl/aux_/include_preprocessed.hpp>

#else

#   if !defined(BOOST_MPL_CFG_NO_FULL_LAMBDA_SUPPORT)
#   if defined(BOOST_MPL_CFG_EXTENDED_TEMPLATE_PARAMETERS_MATCHING)

#   include <b/mpl/limits/arity.hpp>
#   include <b/mpl/aux_/preprocessor/range.hpp>
#   include <b/mpl/aux_/preprocessor/repeat.hpp>
#   include <b/mpl/aux_/preprocessor/params.hpp>
#   include <b/mpl/aux_/nttp_decl.hpp>

#   include <b/preprocessor/seq/fold_left.hpp>
#   include <b/preprocessor/comma_if.hpp>
#   include <b/preprocessor/iterate.hpp>
#   include <b/preprocessor/inc.hpp>
#   include <b/preprocessor/cat.hpp>

#   define AUX778076_ARITY BOOST_PP_INC(BOOST_MPL_LIMIT_METAFUNCTION_ARITY)

namespace boost { namespace mpl { namespace aux {

template< BOOST_MPL_AUX_NTTP_DECL(int, N) > struct arity_tag
{
    typedef char (&type)[(unsigned)N + 1];
};

#   define AUX778076_MAX_ARITY_OP(unused, state, i_) \
    ( BOOST_PP_CAT(C,i_) > 0 ? BOOST_PP_CAT(C,i_) : state ) \
/**/

template<
      BOOST_MPL_PP_PARAMS(AUX778076_ARITY, BOOST_MPL_AUX_NTTP_DECL(int, C))
    >
struct max_arity
{
    BOOST_STATIC_CONSTANT(int, value = 
          BOOST_PP_SEQ_FOLD_LEFT(
              AUX778076_MAX_ARITY_OP
            , -1
            , BOOST_MPL_PP_RANGE(1, AUX778076_ARITY)
            )
        );
};

#   undef AUX778076_MAX_ARITY_OP

arity_tag<0>::type arity_helper(...);

#   define BOOST_PP_ITERATION_LIMITS (1, AUX778076_ARITY)
#   define BOOST_PP_FILENAME_1 <b/mpl/aux_/template_arity.hpp>
#   include BOOST_PP_ITERATE()

template< typename F, BOOST_MPL_AUX_NTTP_DECL(int, N) >
struct template_arity_impl
{
    BOOST_STATIC_CONSTANT(int, value = 
          sizeof(::boost::mpl::aux::arity_helper(type_wrapper<F>(),arity_tag<N>())) - 1
        );
};

#   define AUX778076_TEMPLATE_ARITY_IMPL_INVOCATION(unused, i_, F) \
    BOOST_PP_COMMA_IF(i_) template_arity_impl<F,BOOST_PP_INC(i_)>::value \
/**/

template< typename F >
struct template_arity
{
    BOOST_STATIC_CONSTANT(int, value = (
          max_arity< BOOST_MPL_PP_REPEAT(
              AUX778076_ARITY
            , AUX778076_TEMPLATE_ARITY_IMPL_INVOCATION
            , F
            ) >::value
        ));
        
    typedef mpl::int_<value> type;
};

#   undef AUX778076_TEMPLATE_ARITY_IMPL_INVOCATION

#   undef AUX778076_ARITY

}}}

#   endif // BOOST_MPL_CFG_EXTENDED_TEMPLATE_PARAMETERS_MATCHING
#   else // BOOST_MPL_CFG_NO_FULL_LAMBDA_SUPPORT

#   include <b/mpl/aux_/config/eti.hpp>

namespace boost { namespace mpl { namespace aux {

template< bool >
struct template_arity_impl
{
    template< typename F > struct result_
        : mpl::int_<-1>
    {
    };
};

template<>
struct template_arity_impl<true>
{
    template< typename F > struct result_
        : F::arity
    {
    };
};

template< typename F >
struct template_arity
    : template_arity_impl< ::boost::mpl::aux::has_rebind<F>::value >
        ::template result_<F>
{
};

#if defined(BOOST_MPL_CFG_MSVC_ETI_BUG)
template<>
struct template_arity<int>
    : mpl::int_<-1>
{
};
#endif

}}}

#   endif // BOOST_MPL_CFG_NO_FULL_LAMBDA_SUPPORT

#endif // BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#endif // BOOST_MPL_AUX_TEMPLATE_ARITY_HPP_INCLUDED

///// iteration

#else
#define i_ BOOST_PP_FRAME_ITERATION(1)

template<
      template< BOOST_MPL_PP_PARAMS(i_, typename P) > class F
    , BOOST_MPL_PP_PARAMS(i_, typename T)
    >
typename arity_tag<i_>::type
arity_helper(type_wrapper< F<BOOST_MPL_PP_PARAMS(i_, T)> >, arity_tag<i_>);

#undef i_
#endif // BOOST_PP_IS_ITERATING
