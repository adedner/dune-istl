// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_ISTL_ACCESS_HH
#define DUNE_ISTL_ACCESS_HH

#include <type_traits>
#include <dune/common/hybridutilities.hh>

// these are needed to specialize
#include <array>
#include <dune/common/reservedvector.hh>
#include <dune/common/fvector.hh>

/** \file File providing helper functions to work in ISTL containers using multi-indices
 *
 */

namespace Dune
{

    namespace ISTL
    {

        // Template aliases for using detection idiom.
        template<class C>
        using dynamicIndexAccess_t = decltype(std::declval<C>()[0]);

        template<class C>
        using staticIndexAccess_t = decltype(std::declval<C>()[Dune::Indices::_0]);

        template<class C>
        using resizeMethod_t = decltype(std::declval<C>().resize(0));

        // Short cuts for feature detection
        template<class C>
        using hasDynamicIndexAccess = Dune::Std::is_detected<dynamicIndexAccess_t, std::remove_reference_t<C>>;

        template<class C>
        using hasStaticIndexAccess = Dune::Std::is_detected<staticIndexAccess_t, std::remove_reference_t<C>>;

        template<class C>
        using hasResizeMethod = Dune::Std::is_detected<resizeMethod_t, std::remove_reference_t<C>>;

        template<class C>
        using isDynamicVector = Dune::Std::is_detected<dynamicIndexAccess_t, std::remove_reference_t<C>>;

        template<class C>
        using isStaticVector = Dune::Std::bool_constant<
            Dune::Std::is_detected_v<staticIndexAccess_t, std::remove_reference_t<C>>
            and not Dune::Std::is_detected_v<dynamicIndexAccess_t, std::remove_reference_t<C>>>;

        template<class C>
        using isScalar = Dune::Std::bool_constant<not Dune::Std::is_detected_v<staticIndexAccess_t, std::remove_reference_t<C>>>;

        template<class C>
        using isVector = Dune::Std::bool_constant<Dune::Std::is_detected_v<staticIndexAccess_t, std::remove_reference_t<C>>>;

        /** \interface MultiIndex
            \brief Interface requirements for a MultiIndex

            To access hierarchic containers in ISTL, we provide herlpe
            functions, which give access to an entry or entries, given
            a multi-index. The multi-index should fulfill the
            following requirements:

            * `pop()` pop from out to inner index
            * `size()`

            We further provide specializations for flat indices (just
            a single integer). If you are using a specific
            implementation of a multi-index that does not fulfill the
            above requirements, a simple view-only wrapper should be
            sufficient to use our helpers.
         */

        namespace Impl {
            // helper to get the max size and work around issues with non-const max_size
            template<typename C>
            struct maxSize : public std::integral_constant<unsigned int, 99> {}; // in an case we have to make sure that the recursion ends

            // specialization for std::array, ReservedVector and FieldVector
            template<typename T, std::size_t N>
            struct maxSize<std::array<T,N>> : public std::integral_constant<unsigned int, N> {};

            template<typename T, int n>
            struct maxSize<Dune::ReservedVector<T,n>> : public std::integral_constant<unsigned int, n> {};

            template<typename T, int n>
            struct maxSize<Dune::FieldVector<T,n>> : public std::integral_constant<unsigned int, n> {};

            template<class C>
            inline constexpr unsigned int maxSize_v = maxSize<C>::value;

            /**
               \brief recursive helper function to call a functor f for a list of entries in an ISTL vector
            */
            template <class F, class MultiIndex, std::size_t i, class... Blocks>
            void applyAtIndex (F&& f, MultiIndex const& mi, index_constant<i>, Blocks&&... blocks)
            {
                if constexpr(std::is_integral<MultiIndex>{})
                {
                    f(blocks..., mi);
                }
                else {
                    if constexpr(i < maxSize_v<MultiIndex>) {
                        if (i < mi.size()) {
                            if constexpr(
                                std::conjunction_v<hasDynamicIndexAccess<Blocks>...>
                                ) {
                                // access block directly
                                applyAtIndex(f,mi,index_constant<i+1>{},blocks[mi[i]]...);
                            }
                            else if constexpr(
                                std::conjunction_v<hasStaticIndexAccess<Blocks>...>
                                ) {
                                // access block by searching for the static index equal to mi[i]
                                using Size = decltype(Hybrid::size(std::get<0>(std::make_tuple(blocks...))));
                                return Hybrid::switchCases(std::make_index_sequence<Size::value>{},mi[i],
                                    [&](auto ii) { applyAtIndex(f,mi,index_constant<i+1>{},blocks[ii]...); });
                            }
                            else {
                                // leaf entry in the container
                                f(blocks..., mi);
                            }
                        }
                        else {
                            // end of multi-index
                            f(blocks..., mi);
                        }
                    }
                    else {
                        // end of multi-index (static condition)
                        f(blocks..., mi);
                    }
                }
            }
        } // end namespace Impl

        /**
           \brief call a functor f for a list of entries in an ISTL vector or in a set of ISTL vectors
        */
        template <class F, class MultiIndex, class... Vectors>
        void applyAtIndex (F&& f, MultiIndex const& mi, Vectors&&... vectors)
        {
            constexpr index_constant<0> level;
            Impl::applyAtIndex(f,mi,level,std::forward<Vectors>(vectors)...);
        }

        template <class F, class Indices, class... Vectors>
        void forEachIndex (F&& f, Indices const& indices, Vectors&&... vectors)
        {
            for (auto const& index : indices)
                applyAtIndex(f,index,std::forward<Vectors>(vectors)...);
        }

    } // end namespace ISTL

} // end namespace Dune

#endif // DUNE_ISTL_ACCESS_HH
