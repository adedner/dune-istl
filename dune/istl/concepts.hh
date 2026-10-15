// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_CONCEPTS_HH
#define DUNE_ISTL_CONCEPTS_HH

#include <concepts>
#include <iterator>

#include <dune/common/indices.hh>

/**
  * Some concepts related to data structures from dune-istl
  */
namespace Dune::ISTL::Concept
{
  // A vector is a dense 1d container
  template<class V>
  concept Vector = requires(V vector)
  {
    { vector.size() } -> std::convertible_to<std::size_t>;
    vector[Dune::Indices::_0];
  };

  // A matrix is a dense 2d container
  template<class M>
  concept Matrix = requires(M matrix)
  {
    { matrix.N() } -> std::convertible_to<std::size_t>;
    { matrix.M() } -> std::convertible_to<std::size_t>;
    matrix[Dune::Indices::_0][Dune::Indices::_0];
  };

  // An iterator is something we can increment, compare, and dereference
  template<class It>
  concept Iterator = requires(It it)
  {
    { ++it } -> std::convertible_to<It&>;
    { it++ } -> std::convertible_to<It>;
    { it == it } -> std::convertible_to<bool>;
    { it != it } -> std::convertible_to<bool>;
    { *it } -> std::convertible_to<typename std::iterator_traits<It>::reference>;
  };

  // An indexed iterator is an iterator with a method .index()
  template<class It>
  concept IndexedIterator = Iterator<It> && requires(It it)
  {
    { it.index() } -> std::convertible_to<std::size_t>;
  };

  // A range is something we can iterator over
  template<class R>
  concept Range = requires(R range)
  {
    { range.begin() } -> Iterator;
    { range.end() } -> Iterator;
  };

  // A nested range is a Range whose elements are ranges too
  template<class M>
  concept NestedRange = Range<M> && requires(M matrix)
  {
    { *matrix.begin() } -> Range;
  };

  template<class V>
  concept IndexedRange = Range<V> && requires(V vector)
  {
    { vector.size() } -> std::convertible_to<std::size_t>;
    { vector.begin() } -> IndexedIterator;
  };

  template<class M>
  concept NestedIndexedRange = NestedRange<M> && requires(M matrix)
  {
    { matrix.N() } -> std::convertible_to<std::size_t>;
    { matrix.M() } -> std::convertible_to<std::size_t>;
    { matrix.begin() } -> IndexedIterator;
    { matrix.begin()->begin() } -> IndexedIterator;
  };

  // A scalar is neither of the above
  template<class S>
  concept Scalar = not (Vector<S> || Matrix<S> || Range<S>);

} // end Dune::ISTL::Concept

#endif // DUNE_ISTL_CONCEPTS_HH