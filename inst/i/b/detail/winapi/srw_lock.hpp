/*
 * Copyright 2017 Andrey Semashev
 *
 * Distributed under the Boost Software License, Version 1.0.
 * See http://www.boost.org/LICENSE_1_0.txt
 *
 * This header is deprecated, use b/winapi/srw_lock.hpp instead.
 */

#ifndef BOOST_DETAIL_WINAPI_SRW_LOCK_HPP
#define BOOST_DETAIL_WINAPI_SRW_LOCK_HPP

#include <b/config/header_deprecated.hpp>

BOOST_HEADER_DEPRECATED("<b/winapi/srw_lock.hpp>")

#include <b/winapi/srw_lock.hpp>
#include <b/detail/winapi/detail/deprecated_namespace.hpp>

#ifdef BOOST_HAS_PRAGMA_ONCE
#pragma once
#endif

// Deprecated
#define BOOST_DETAIL_WINAPI_SRWLOCK_INIT BOOST_WINAPI_SRWLOCK_INIT

#endif // BOOST_DETAIL_WINAPI_SRW_LOCK_HPP
