///////////////////////////////////////////////////////////////////////////////
/// \file algorithm.hpp
///   Includes the range-based versions of the algorithms in the
///   C++ standard header file <algorithm>
//
/////////////////////////////////////////////////////////////////////////////

// Copyright 2009 Neil Groves.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// Acknowledgements:
// This code uses combinations of ideas, techniques and code snippets
// from: Thorsten Ottosen, Eric Niebler, Jeremy Siek,
// and Vladimir Prus'
//
// The original mutating algorithms that served as the first version
// were originally written by Vladimir Prus'
// <ghost@cs.msu.su> code from Boost Wiki

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef BOOST_RANGE_ALGORITHM_HPP_INCLUDED_01012009
#define BOOST_RANGE_ALGORITHM_HPP_INCLUDED_01012009

#include <b/range/concepts.hpp>
#include <b/range/iterator_range.hpp>
#include <b/range/difference_type.hpp>
#include <b/range/detail/range_return.hpp>
#include <b/iterator/iterator_traits.hpp>
#include <b/next_prior.hpp>
#include <algorithm>

// Non-mutating algorithms
#include <b/range/algorithm/adjacent_find.hpp>
#include <b/range/algorithm/count.hpp>
#include <b/range/algorithm/count_if.hpp>
#include <b/range/algorithm/equal.hpp>
#include <b/range/algorithm/for_each.hpp>
#include <b/range/algorithm/find.hpp>
#include <b/range/algorithm/find_end.hpp>
#include <b/range/algorithm/find_first_of.hpp>
#include <b/range/algorithm/find_if.hpp>
#include <b/range/algorithm/lexicographical_compare.hpp>
#include <b/range/algorithm/mismatch.hpp>
#include <b/range/algorithm/search.hpp>
#include <b/range/algorithm/search_n.hpp>

// Mutating algorithms
#include <b/range/algorithm/copy.hpp>
#include <b/range/algorithm/copy_backward.hpp>
#include <b/range/algorithm/fill.hpp>
#include <b/range/algorithm/fill_n.hpp>
#include <b/range/algorithm/generate.hpp>
#include <b/range/algorithm/inplace_merge.hpp>
#include <b/range/algorithm/merge.hpp>
#include <b/range/algorithm/nth_element.hpp>
#include <b/range/algorithm/partial_sort.hpp>
#include <b/range/algorithm/partial_sort_copy.hpp>
#include <b/range/algorithm/partition.hpp>
#include <b/range/algorithm/random_shuffle.hpp>
#include <b/range/algorithm/remove.hpp>
#include <b/range/algorithm/remove_copy.hpp>
#include <b/range/algorithm/remove_copy_if.hpp>
#include <b/range/algorithm/remove_if.hpp>
#include <b/range/algorithm/replace.hpp>
#include <b/range/algorithm/replace_copy.hpp>
#include <b/range/algorithm/replace_copy_if.hpp>
#include <b/range/algorithm/replace_if.hpp>
#include <b/range/algorithm/reverse.hpp>
#include <b/range/algorithm/reverse_copy.hpp>
#include <b/range/algorithm/rotate.hpp>
#include <b/range/algorithm/rotate_copy.hpp>
#include <b/range/algorithm/sort.hpp>
#include <b/range/algorithm/stable_partition.hpp>
#include <b/range/algorithm/stable_sort.hpp>
#include <b/range/algorithm/transform.hpp>
#include <b/range/algorithm/unique.hpp>
#include <b/range/algorithm/unique_copy.hpp>

// Binary search
#include <b/range/algorithm/binary_search.hpp>
#include <b/range/algorithm/equal_range.hpp>
#include <b/range/algorithm/lower_bound.hpp>
#include <b/range/algorithm/upper_bound.hpp>

// Set operations of sorted ranges
#include <b/range/algorithm/set_algorithm.hpp>

// Heap operations
#include <b/range/algorithm/heap_algorithm.hpp>

// Minimum and Maximum
#include <b/range/algorithm/max_element.hpp>
#include <b/range/algorithm/min_element.hpp>

// Permutations
#include <b/range/algorithm/permutation.hpp>

#endif // include guard

