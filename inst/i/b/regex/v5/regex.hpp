/*
 *
 * Copyright (c) 1998-2002
 * John Maddock
 *
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         regex.cpp
  *   VERSION      see <b/version.hpp>
  *   DESCRIPTION: Declares boost::basic_regex<> and associated
  *                functions and classes. This header is the main
  *                entry point for the template regex code.
  */

#ifndef BOOST_RE_REGEX_HPP_INCLUDED
#define BOOST_RE_REGEX_HPP_INCLUDED

#ifdef __cplusplus

// what follows is all C++ don't include in C builds!!

#include <b/regex/config.hpp>
#include <b/regex/v5/regex_workaround.hpp>
#include <b/regex_fwd.hpp>
#include <b/regex/regex_traits.hpp>
#include <b/regex/v5/error_type.hpp>
#include <b/regex/v5/match_flags.hpp>
#include <b/regex/v5/regex_raw_buffer.hpp>
#include <b/regex/pattern_except.hpp>
#include <b/regex/v5/char_regex_traits.hpp>
#include <b/regex/v5/states.hpp>
#include <b/regex/v5/regbase.hpp>
#include <b/regex/v5/basic_regex.hpp>
#include <b/regex/v5/basic_regex_creator.hpp>
#include <b/regex/v5/basic_regex_parser.hpp>
#include <b/regex/v5/sub_match.hpp>
#include <b/regex/v5/regex_format.hpp>
#include <b/regex/v5/match_results.hpp>
#include <b/regex/v5/perl_matcher.hpp>

namespace boost{
#ifdef BOOST_REGEX_NO_FWD
typedef basic_regex<char, regex_traits<char> > regex;
#ifndef BOOST_NO_WREGEX
typedef basic_regex<wchar_t, regex_traits<wchar_t> > wregex;
#endif
#endif

typedef match_results<const char*> cmatch;
typedef match_results<std::string::const_iterator> smatch;
#ifndef BOOST_NO_WREGEX
typedef match_results<const wchar_t*> wcmatch;
typedef match_results<std::wstring::const_iterator> wsmatch;
#endif

} // namespace boost

#include <b/regex/v5/regex_match.hpp>
#include <b/regex/v5/regex_search.hpp>
#include <b/regex/v5/regex_iterator.hpp>
#include <b/regex/v5/regex_token_iterator.hpp>
#include <b/regex/v5/regex_grep.hpp>
#include <b/regex/v5/regex_replace.hpp>
#include <b/regex/v5/regex_merge.hpp>
#include <b/regex/v5/regex_split.hpp>

#endif  // __cplusplus

#endif  // include































