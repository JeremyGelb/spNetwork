/*
 [auto_generated]
 b/numeric/odeint/util/n_ary_helper.hpp

 Macros to generate scale_sumN and for_eachN functors.

 Copyright 2013 Karsten Ahnert
 Copyright 2013 Mario Mulansky
 Copyright 2013 Pascal Germroth

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_N_ARY_HELPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_UTIL_N_ARY_HELPER_HPP_INCLUDED

#include <b/preprocessor/repetition.hpp>

// like BOOST_PP_ENUM_SHIFTED but with a comma in front like _TRAILING
#define BOOST_ODEINT_ENUM_TRAILING_SHIFTED_PARAMS(count, param) \
    BOOST_PP_COMMA_IF(BOOST_PP_DEC(count)) \
    BOOST_PP_ENUM_SHIFTED_PARAMS(count, param)

#define BOOST_ODEINT_ENUM_TRAILING_SHIFTED_BINARY_PARAMS(count, p1, p2) \
    BOOST_PP_COMMA_IF(BOOST_PP_DEC(count)) \
    BOOST_PP_ENUM_SHIFTED_BINARY_PARAMS(count, p1, p2)

// like BOOST_PP_ENUM_SHIFTED_BINARY_PARAMS(n, p1, p2) but p2 is shifted left.
// generate "p1 ## 0 = p2, p1 ## 1 = p3 ## 0, p1 ## 2 = p3 ## 1"
#define BOOST_ODEINT_ENUM_LSHIFTED_BINARY_PARAMS(count, p1, p2, p3) \
    BOOST_PP_ENUM(count, BOOST_ODEINT_ENUM_LSHIFTED_BINARY_PARAMS_, (p1, p2, p3))
#define BOOST_ODEINT_ENUM_LSHIFTED_BINARY_PARAMS_(z, n, data) \
    BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(3, 0, data), n) \
    BOOST_PP_IF(n, \
        BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(3, 2, data), BOOST_PP_DEC(n)), \
        BOOST_PP_TUPLE_ELEM(3, 1, data))

// like BOOST_PP_ENUM_BINARY_PARAMS(n, p1, p2) but with statements.
// "p1 ## 0 p2 ## 0 ; p1 ## 1 p2 ## 1 ; ..."
#define BOOST_ODEINT_ENUM_BINARY_STATEMENTS(count, p1, p2) \
    BOOST_PP_REPEAT(count, BOOST_ODEINT_ENUM_BINARY_STATEMENTS_, (p1, p2))
#define BOOST_ODEINT_ENUM_BINARY_STATEMENTS_(z, n, data) \
    BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2, 0, data), n) \
    BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2, 1, data), n) ;

// like BOOST_PP_ENUM_BINARY_PARAMS(n, p1, p2) but p2 is in parens.
// "p1 ## 0 (p2 ## 0) , p1 ## 1 (p2 ## 1) , ..."
#define BOOST_ODEINT_ENUM_UNARY_CALLS(count, p1, p2) \
    BOOST_PP_ENUM(count, BOOST_ODEINT_ENUM_UNARY_CALLS_, (p1, p2))
#define BOOST_ODEINT_ENUM_SHIFTED_UNARY_CALLS(count, p1, p2) \
    BOOST_PP_ENUM_SHIFTED(count, BOOST_ODEINT_ENUM_UNARY_CALLS_, (p1, p2))
#define BOOST_ODEINT_ENUM_TRAILING_SHIFTED_UNARY_CALLS(count, p1, p2) \
    BOOST_PP_COMMA_IF(BOOST_PP_DEC(count)) \
    BOOST_PP_ENUM_SHIFTED(count, BOOST_ODEINT_ENUM_UNARY_CALLS_, (p1, p2))
#define BOOST_ODEINT_ENUM_UNARY_CALLS_(z, n, data) \
      BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2, 0, data), n) \
    ( BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2, 1, data), n) )


// maximum arity + 1 for scale_sum and for_each
#define BOOST_ODEINT_N_ARY_MAX 16


// generate scale_sum1 to scale_sumN, operator body generated by macro(N)
#define BOOST_ODEINT_GEN_SCALE_SUM(macro) \
    BOOST_PP_REPEAT_FROM_TO(1, BOOST_ODEINT_N_ARY_MAX, BOOST_ODEINT_GEN_SCALE_SUM_, macro)
#define BOOST_ODEINT_GEN_SCALE_SUM_(z, n, macro) \
    template< BOOST_ODEINT_ENUM_LSHIFTED_BINARY_PARAMS(n, class Fac, = double, = Fac) > \
    struct BOOST_PP_CAT(scale_sum, n) \
    { \
        BOOST_ODEINT_ENUM_BINARY_STATEMENTS(n, const Fac, m_alpha) \
        \
        BOOST_PP_CAT(scale_sum, n) \
        ( BOOST_PP_ENUM_BINARY_PARAMS(n, Fac, alpha) ) \
        : BOOST_ODEINT_ENUM_UNARY_CALLS(n, m_alpha, alpha) {} \
        \
        template< BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), class T) > \
        void operator()( T0 &t0 \
            BOOST_ODEINT_ENUM_TRAILING_SHIFTED_BINARY_PARAMS(BOOST_PP_INC(n), const T, &t) \
        ) const \
        { macro(n) } \
        typedef void result_type; \
    };

// generate for_each1 to for_eachN, body generated by macro(N)
#define BOOST_ODEINT_GEN_FOR_EACH(macro) \
    BOOST_PP_REPEAT_FROM_TO(1, BOOST_ODEINT_N_ARY_MAX, BOOST_ODEINT_GEN_FOR_EACH_, macro)
#define BOOST_ODEINT_GEN_FOR_EACH_(z, n, macro) \
    template< BOOST_PP_ENUM_PARAMS(n, class S) , class Op > \
    static void for_each##n ( BOOST_PP_ENUM_BINARY_PARAMS(n, S, &s) , Op op ) \
    { macro(n) }


#endif
