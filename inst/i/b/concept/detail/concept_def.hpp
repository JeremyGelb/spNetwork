// Copyright David Abrahams 2006. Distributed under the Boost
// Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_CONCEPT_DETAIL_CONCEPT_DEF_DWA200651_HPP
# define BOOST_CONCEPT_DETAIL_CONCEPT_DEF_DWA200651_HPP
# include <b/preprocessor/seq/for_each_i.hpp>
# include <b/preprocessor/seq/enum.hpp>
# include <b/preprocessor/comma_if.hpp>
# include <b/preprocessor/cat.hpp>
#endif // BOOST_CONCEPT_DETAIL_CONCEPT_DEF_DWA200651_HPP

// BOOST_concept(SomeName, (p1)(p2)...(pN))
//
// Expands to "template <class p1, class p2, ...class pN> struct SomeName"
//
// Also defines an equivalent SomeNameConcept for backward compatibility.
// Maybe in the next release we can kill off the "Concept" suffix for good.
# define BOOST_concept(name, params)                                            \
    template < BOOST_PP_SEQ_FOR_EACH_I(BOOST_CONCEPT_typename,~,params) >       \
    struct name; /* forward declaration */                                      \
                                                                                \
    template < BOOST_PP_SEQ_FOR_EACH_I(BOOST_CONCEPT_typename,~,params) >       \
    struct BOOST_PP_CAT(name,Concept)                                           \
      : name< BOOST_PP_SEQ_ENUM(params) >                                       \
    {                                                                           \
    };                                                                          \
                                                                                \
    template < BOOST_PP_SEQ_FOR_EACH_I(BOOST_CONCEPT_typename,~,params) >       \
    struct name                                                                
    
// Helper for BOOST_concept, above.
# define BOOST_CONCEPT_typename(r, ignored, index, t) \
    BOOST_PP_COMMA_IF(index) typename t

