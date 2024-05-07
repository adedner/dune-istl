// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_SIZEUTILITIES_HH
#define DUNE_ISTL_SIZEUTILITIES_HH

#include <type_traits>

#include <dune/common/rangeutilities.hh>
#include <dune/common/typeutilities.hh>

namespace Dune {
namespace Hybrid {
namespace Impl {

template <class M>
using DynamicMatrixIndexAccessConcept = decltype((std::declval<M>()[0][0], true));

template <class V>
using DynamicVectorIndexAccessConcept = decltype((std::declval<V>()[0], true));

// if matrix has dynamic element access and provides method N()
template <class Matrix,
  DynamicMatrixIndexAccessConcept<Matrix> = true>
auto numRows (const Matrix& matrix, Dune::PriorityTag<2>)
  -> decltype(std::size_t(matrix.N()))
{
  return matrix.N();
}

// if matrix provides static constexpr method N()
template <class Matrix>
auto numRows (const Matrix& matrix, Dune::PriorityTag<1>)
  -> std::integral_constant<std::size_t, std::size_t(Matrix::N())>
{
  return {};
}

// if matrix has dynamic element access and provides method M()
template <class Matrix,
  DynamicMatrixIndexAccessConcept<Matrix> = true>
auto numCols (const Matrix& matrix, Dune::PriorityTag<2>)
  -> decltype(std::size_t(matrix.M()))
{
  return matrix.M();
}

// if matrix provides static constexpr method M()
template <class Matrix>
auto numCols (const Matrix& matrix, Dune::PriorityTag<1>)
  -> std::integral_constant<std::size_t, std::size_t(Matrix::M())>
{
  return {};
}

// if vector has dynamic element access and provides method size()
template <class Vector,
  DynamicVectorIndexAccessConcept<Vector> = true>
auto numEntries (const Vector& vector, Dune::PriorityTag<2>)
  -> decltype(std::size_t(vector.size()))
{
  return vector.size();
}

// if vector provides static constexpr method size()
template <class Vector>
auto numEntries (const Vector& vector, Dune::PriorityTag<1>)
  -> std::integral_constant<std::size_t, std::size_t(Vector::size())>
{
  return {};
}

} // end namespace Impl


/// \brief Return the number of matrix rows
/**
 * The number of rows is returned either as integer or integral-constant
 * depending on whether the matrix has dynamic or static size and can be iterated
 * dynamically or statically only.
 *
 * Dynamic size information is preferred over static size information, if element
 * access allows dynamic indices.
 **/
template <class Matrix>
auto numRows (const Matrix& matrix)
{
  return Impl::numRows(matrix, Dune::PriorityTag<42>{});
}

/// \brief Return the number of matrix columns
/**
 * The number of columns is returned either as integer or integral-constant
 * depending on whether the matrix has dynamic or static size and can be iterated
 * dynamically or statically only.
 *
 * Dynamic size information is preferred over static size information, if element
 * access allows dynamic indices.
 **/
template <class Matrix>
auto numCols (const Matrix& matrix)
{
  return Impl::numCols(matrix, Dune::PriorityTag<42>{});
}

/// \brief Return the number of vector entries
/**
 * The number of entries is returned either as integer or integral-constant
 * depending on whether the vector has dynamic or static size and can be iterated
 * dynamically or statically only.
 *
 * Dynamic size information is preferred over static size information, if element
 * access allows dynamic indices.
 **/
template <class Vector>
auto numEntries (const Vector& vector)
{
  return Impl::numEntries(vector, Dune::PriorityTag<42>{});
}

/// \brief Return a dynamic or static index-range over the matrix rows
/**
 * \relates numRows
 **/
template <class Matrix>
auto rows (const Matrix& matrix)
{
  return Dune::range(numRows(matrix));
}

/// \brief Return a dynamic or static index-range over the matrix columns
/**
 * \relates numCols
 **/
template <class Matrix>
auto cols (const Matrix& matrix)
{
  return Dune::range(numCols(matrix));
}

/// \brief Return a dynamic or static index-range over the vector entries
/**
 * \relates numEntries
 **/
template <class Vector>
auto entries (const Vector& vector)
{
  return Dune::range(numEntries(vector));
}

}} // end namespace Dune::Hybrid

#endif // DUNE_ISTL_SIZEUTILITIES_HH
