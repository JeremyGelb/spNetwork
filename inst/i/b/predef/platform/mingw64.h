/*
Copyright Rene Rivera 2008-2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_PLAT_MINGW64_H
#define BOOST_PREDEF_PLAT_MINGW64_H

#include <b/predef/version_number.h>
#include <b/predef/make.h>

/* tag::reference[]
= `BOOST_PLAT_MINGW64`

https://mingw-w64.org/[MinGW-w64] platform.
Version number available as major, minor, and patch.

[options="header"]
|===
| {predef_symbol} | {predef_version}

| `+__MINGW64__+` | {predef_detection}

| `+__MINGW64_VERSION_MAJOR+`, `+__MINGW64_VERSION_MINOR+` | V.R.0
|===
*/ // end::reference[]

#define BOOST_PLAT_MINGW64 BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__MINGW64__)
#   include <_mingw.h>
#   if !defined(BOOST_PLAT_MINGW64_DETECTION) && (defined(__MINGW64_VERSION_MAJOR) && defined(__MINGW64_VERSION_MINOR))
#       define BOOST_PLAT_MINGW64_DETECTION \
            BOOST_VERSION_NUMBER(__MINGW64_VERSION_MAJOR,__MINGW64_VERSION_MINOR,0)
#   endif
#   if !defined(BOOST_PLAT_MINGW64_DETECTION)
#       define BOOST_PLAT_MINGW64_DETECTION BOOST_VERSION_NUMBER_AVAILABLE
#   endif
#endif

#ifdef BOOST_PLAT_MINGW64_DETECTION
#   define BOOST_PLAT_MINGW64_AVAILABLE
#   if defined(BOOST_PREDEF_DETAIL_PLAT_DETECTED)
#       define BOOST_PLAT_MINGW64_EMULATED BOOST_PLAT_MINGW64_DETECTION
#   else
#       undef BOOST_PLAT_MINGW64
#       define BOOST_PLAT_MINGW64 BOOST_PLAT_MINGW64_DETECTION
#   endif
#   include <b/predef/detail/platform_detected.h>
#endif

#define BOOST_PLAT_MINGW64_NAME "MinGW-w64"

#endif

#include <b/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_PLAT_MINGW64,BOOST_PLAT_MINGW64_NAME)

#ifdef BOOST_PLAT_MINGW64_EMULATED
#include <b/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_PLAT_MINGW64_EMULATED,BOOST_PLAT_MINGW64_NAME)
#endif
