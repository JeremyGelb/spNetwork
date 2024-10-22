/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_COMPILER_MICROTEC_H
#define BOOST_PREDEF_COMPILER_MICROTEC_H

#include <b/predef/version_number.h>
#include <b/predef/make.h>

/* tag::reference[]
= `BOOST_COMP_MRI`

http://www.mentor.com/microtec/[Microtec C/{CPP}] compiler.

[options="header"]
|===
| {predef_symbol} | {predef_version}

| `+_MRI+` | {predef_detection}
|===
*/ // end::reference[]

#define BOOST_COMP_MRI BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(_MRI)
#   define BOOST_COMP_MRI_DETECTION BOOST_VERSION_NUMBER_AVAILABLE
#endif

#ifdef BOOST_COMP_MRI_DETECTION
#   if defined(BOOST_PREDEF_DETAIL_COMP_DETECTED)
#       define BOOST_COMP_MRI_EMULATED BOOST_COMP_MRI_DETECTION
#   else
#       undef BOOST_COMP_MRI
#       define BOOST_COMP_MRI BOOST_COMP_MRI_DETECTION
#   endif
#   define BOOST_COMP_MRI_AVAILABLE
#   include <b/predef/detail/comp_detected.h>
#endif

#define BOOST_COMP_MRI_NAME "Microtec C/C++"

#endif

#include <b/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_COMP_MRI,BOOST_COMP_MRI_NAME)

#ifdef BOOST_COMP_MRI_EMULATED
#include <b/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_COMP_MRI_EMULATED,BOOST_COMP_MRI_NAME)
#endif
