# /* **************************************************************************
#  *                                                                          *
#  *     (C) Copyright Paul Mensonides 2002.
#  *     Distributed under the Boost Software License, Version 1.0. (See
#  *     accompanying file LICENSE_1_0.txt or copy at
#  *     http://www.boost.org/LICENSE_1_0.txt)
#  *                                                                          *
#  ************************************************************************** */
#
# /* Revised by Edward Diener (2020) */
#
# /* See http://www.boost.org for most recent version. */
#
# ifndef BOOST_PREPROCESSOR_REPETITION_DEDUCE_R_HPP
# define BOOST_PREPROCESSOR_REPETITION_DEDUCE_R_HPP
#
# include <b/preprocessor/config/config.hpp>
#
# if ~BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_STRICT()
#
# include <b/preprocessor/detail/auto_rec.hpp>
# include <b/preprocessor/repetition/for.hpp>
#
# /* BOOST_PP_DEDUCE_R */
#
# define BOOST_PP_DEDUCE_R() BOOST_PP_AUTO_REC(BOOST_PP_FOR_P, 256)
#
# else
#
# /* BOOST_PP_DEDUCE_R */
#
# include <b/preprocessor/arithmetic/dec.hpp>
# include <b/preprocessor/detail/auto_rec.hpp>
# include <b/preprocessor/repetition/for.hpp>
# include <b/preprocessor/config/limits.hpp>
#
# if BOOST_PP_LIMIT_FOR == 256
# define BOOST_PP_DEDUCE_R() BOOST_PP_DEC(BOOST_PP_AUTO_REC(BOOST_PP_FOR_P, 256))
# elif BOOST_PP_LIMIT_FOR == 512
# define BOOST_PP_DEDUCE_R() BOOST_PP_DEC(BOOST_PP_AUTO_REC(BOOST_PP_FOR_P, 512))
# elif BOOST_PP_LIMIT_FOR == 1024
# define BOOST_PP_DEDUCE_R() BOOST_PP_DEC(BOOST_PP_AUTO_REC(BOOST_PP_FOR_P, 1024))
# else
# error Incorrect value for the BOOST_PP_LIMIT_FOR limit
# endif
#
# endif
#
# endif
