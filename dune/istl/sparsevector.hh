// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_SPARSEVECTOR_HH
#define DUNE_ISTL_SPARSEVECTOR_HH

#include <iterator>
#include <span>
#include <utility>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/std/extents.hh>
#include <dune/common/std/mdarray.hh>
#include <dune/common/std/no_unique_address.hh>
#include <dune/istl/bvector.hh>

namespace Dune {
namespace Impl {
template <class It>
struct SparseVectorIterator
  : public Dune::IteratorFacadeForTraits<SparseVectorIterator<It>,
      Dune::DefaultIteratorTraits<std::forward_iterator_tag, decltype(std::declval<It>()->second)>>
{
  SparseVectorIterator (It it) : it_(it) {}
  auto& operator*() const { return it_->second; }
  SparseVectorIterator& operator++() { ++it_; return *this; }
  SparseVectorIterator operator++(int) { SparseVectorIterator tmp(*this); ++(*this); return tmp; }
  bool operator==(const SparseVectorIterator& other) const { return it_ == other.it_; };
  std::size_t index() const { return it_->first; }

  It it_;
};

} // end namespace Impl


template <class T, std::size_t S = std::dynamic_extent, std::size_t C = S>
class SparseReservedVector
{
  using self_type = SparseReservedVector;

public:
  using block_type = T;
  using field_type = typename BlockType<block_type>::field_type;
  using value_type = T;
  using reference = T&;
  using const_reference = T const&;

public:
  constexpr SparseReservedVector ()
    : size_{S == std::dynamic_extent ? 0 : S}
    , data_{C == std::dynamic_extent ? size_.extent(0) : C}
  {}

  explicit constexpr SparseReservedVector (std::size_t size)
        requires (S == std::dynamic_extent)
    : size_{size}
    , data_{C == std::dynamic_extent ? size_.extent(0) : C}
  {}

  explicit constexpr SparseReservedVector (std::size_t size, std::size_t capacity)
        requires (S == std::dynamic_extent && C == std::dynamic_extent)
    : size_{size}
    , data_{capacity}
  {}

public:
  /// \name Iterators
  /// @{

  constexpr auto begin() noexcept { return Impl::SparseVectorIterator{data_.container_data()}; }
  constexpr auto end() noexcept { return Impl::SparseVectorIterator{data_.container_data() + nnz_}; }
  constexpr auto begin() const noexcept { return Impl::SparseVectorIterator{data_.container_data()}; }
  constexpr auto end() const noexcept { return Impl::SparseVectorIterator{data_.container_data() + nnz_}; }
  constexpr auto cbegin() const noexcept { return Impl::SparseVectorIterator{data_.container_data()}; }
  constexpr auto cend() const noexcept { return Impl::SparseVectorIterator{data_.container_data() + nnz_}; }

  /// @}


  /// \name Capacity
  /// @{

  constexpr std::size_t size() const { return size_.extent(0); }
  constexpr std::size_t capacity() const { return data_.extent(0); }
  constexpr std::size_t nnz() const { return nnz_; }

  /// @}


  /// \name Modifiers
  /// @{

  void insert(std::size_t pos, T value)
  {
    assert(nnz_ < capacity());
    if (nnz_ < capacity())
      data_[nnz_++] = std::pair<std::size_t,T>{pos,std::move(value)};
  }

  /// @}


  /// \name Vector-space operations
  /// @{

  self_type& operator= (const field_type& scalar)
  {
    for (std::size_t i = 0; i < nnz_; ++i)
      data_[i].second = scalar;
  }

  self_type& operator+= (const field_type& scalar)
  {
    for (std::size_t i = 0; i < nnz_; ++i)
      data_[i].second += scalar;
  }

  self_type& operator-= (const field_type& scalar)
  {
    for (std::size_t i = 0; i < nnz_; ++i)
      data_[i].second -= scalar;
  }

  self_type& operator*= (const field_type& scalar)
  {
    for (std::size_t i = 0; i < nnz_; ++i)
      data_[i].second *= scalar;
  }

  self_type& operator/= (const field_type& scalar)
  {
    for (std::size_t i = 0; i < nnz_; ++i)
      data_[i].second /= scalar;
  }

  //! Two norm sqrt(sum over squared values of entries)
  typename FieldTraits<field_type>::real_type two_norm () const
  {
    using std::sqrt;
    return sqrt(two_norm2());
  }

  //! Square of the two-norm (the sum over the squared values of the entries)
  typename FieldTraits<field_type>::real_type two_norm2 () const
  {
    typename FieldTraits<field_type>::real_type sum = 0;
    for (auto it = begin(); it != end(); ++it)
      sum += Impl::asVector(it->second).two_norm2();
    return sum;
  }

  //! Infinity norm (maximum of absolute values of entries)
  typename FieldTraits<field_type>::real_type infinity_norm () const
        requires (!HasNaN<field_type>::value)
  {
    using real_type = typename FieldTraits<field_type>::real_type;
    using std::max;

    real_type norm = 0;
    for (auto it = begin(); it != end(); ++it) {
      real_type const a = Impl::asVector(it->second).infinity_norm();
      norm = max(a, norm);
    }
    return norm;
  }

  /// @}

private:
  DUNE_NO_UNIQUE_ADDRESS Std::extents<std::size_t,S> size_;
  Std::mdarray<std::pair<std::size_t,T>, Std::extents<std::size_t,C>> data_;
  std::size_t nnz_ = 0;
};

} // end namespace Dune

#endif // DUNE_ISTL_SPARSEVECTOR_HH
