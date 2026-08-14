// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_ISTL_ACCESS_HH
#define DUNE_ISTL_ACCESS_HH

#include <type_traits>
#include <dune/common/hybridutilities.hh>

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

        /**
           \brief all a functor f for a list of entries in an ISTL vector
         */
        template <class MultiIndex, class Block, class F, std::size_t i = 0>
        void applyAtIndex (MultiIndex const& mi, Block&& block, F&& f, index_constant<i> = {})
        {
            if constexpr(std::is_integral<MultiIndex>{})
            {
                f(block, mi);
            }
            else {
                if constexpr(i < mi.max_size()) {
                    if (i < mi.size()) {
                        if constexpr(hasDynamicIndexAccess<Block>{}) {
                            // access block directly
                            applyAtIndex(mi,block[mi[i]],f,index_constant<i+1>{});
                        }
                        else if constexpr(hasStaticIndexAccess<Block>{}) {
                            // access block by searching for the static index equal to mi[i]
                            using Size = decltype(Hybrid::size(block));
                            return Hybrid::switchCases(std::make_index_sequence<Size::value>{},mi[i],
                                [&](auto ii) { applyAtIndex(mi,block[ii],f,index_constant<i+1>{}); });
                        }
                        else {
                            // leaf entry in the container
                            f(block, mi);
                        }
                    }
                    else {
                        // end of multi-index
                        f(block, mi);
                    }
                }
                else {
                    // end of multi-index (static condition)
                    f(block, mi);
                }
            }
        }

        template <class Indices, class Container, class F>
        void forEachIndex (Indices const& indices, Container&& container, F&& f)
        {
            for (auto const& index : indices)
                applyAtIndex(index,container,f);
        }

        // template <class Indices, class Container, class F>
        // void hierarchicForEach (Container&& container, F&& f)
        // {
        //     for (auto const& index : indices)
        //         applyAtIndex(index,container,f);
        // }

    } // end namespace ISTL

} // end namespace Dune

#endif // DUNE_ISTL_ACCESS_HH
