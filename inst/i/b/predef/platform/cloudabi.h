/*
  Copyright 2017 James E. King, III
  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_PLAT_CLOUDABI_H
#define BOOST_PREDEF_PLAT_CLOUDABI_H

#include <b/predef/version_number.h>
#include <b/predef/make.h>

/* tag::reference[]
= `BOOST_PLAT_CLOUDABI`

https://github.com/NuxiNL/cloudabi[CloudABI] platform.

[options="header"]
|===
| {predef_symbol} | {predef_version}

| `+__CloudABI__+` | {predef_detection}
|===
*/ // end::reference[]

#define BOOST_PLAT_CLOUDABI BOOST_VERSION_NUMBER_NOT_AVAILABLE

#if defined(__CloudABI__)
#   undef BOOST_PLAT_CLOUDABI
#   define BOOST_PLAT_CLOUDABI BOOST_VERSION_NUMBER_AVAILABLE
#endif

#if BOOST_PLAT_CLOUDABI
#   define BOOST_PLAT_CLOUDABI_AVAILABLE
#   include <b/predef/detail/platform_detected.h>
#endif

#define BOOST_PLAT_CLOUDABI_NAME "CloudABI"

#endif

#include <b/predef/detail/test.h>
BOOST_PREDEF_DECLARE_TEST(BOOST_PLAT_CLOUDABI,BOOST_PLAT_CLOUDABI_NAME)
