/*
  we want to check different options about how to evaluate
  weighted/parallel scalar products.

  The main thing is, that the SP evaluate only a subset of entries.

  We see different algorithmic approaches:

  == (A) (nested) bool vector

  Given a nested bool vector that matches the structure of the ISTL
  vector, we simply mask the scalar product with the bool entries.

  == (B) skip list

  We define a list `std::vector<MultiIndex>` that lists which entries
  should not be considered in the scalar product.

  == (C) "reverse" skip list

  We use the same skip list as before.

  1. compute full scalar product
  2. evaluate a "sparse" scalar product only on the entries of the skip list
  3. subtract (2) from (1)

  == weighted scalar product with sparse diagonal matrix

  Conceptually the nicest approach is to write these as weighted inner
  product. This allows a slightly more general setting. We could
  define a very special diagonal matrix format that is by default 1
  and lists only those entries that differ.

 */

#include <type_traits>
#include <vector>
#include <iostream>
#include <cassert>

#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/access.hh>
#include <dune/istl/foreach.hh>

namespace Dune {
namespace ScalarProductHelper {

   /** \brief Traverse a blocked vector and call a functor at each scalar entry
    *
    *  The functor `f` is assumed to have the signature
    *
    *    void(auto&& entry, std::size_t offset)
    *
    *  taking a scalar entry and the current flat offset (index)
    *  of this position.
    *
    *  It returns the total number of scalar entries. Similar to `dimension()` for
    *  some DUNE vector types.
    */
    template <class F, class... Vectors>
    void applyToVector(F&& f, Vectors&&... vectors)
    {
        // are we at the leaf?
        static const bool isScalar =
            std::conjunction_v<Impl::IsScalar<std::decay_t<Vectors>>...>;
        // check consistency between all vectors
        static_assert(
            std::conjunction_v<
                std::conditional_t<
                    Impl::IsScalar<std::decay_t<Vectors>>::value == isScalar,
                    std::true_type, std::false_type>...>
            );
        if constexpr(isScalar)
        {
            f(vectors...);
        }
        else
        {
            // get size
            auto sz = std::make_tuple(ForEach::size(vectors)...);

            // check consistency of size
            // TODO
            Hybrid::forEach(Dune::range(std::get<0>(sz)), [&](auto i) {
                applyToVector(f, vectors[i]...);
            });
        }
    }

    template<typename V, typename B>
    auto dot(const V& a, const V& b, const B& boolVector)
    {
        double sp(0);
        applyToVector(
            [&sp](double a, double b, bool w)
            {
                // only evaluate if w == true
                sp += w?a*b:0.0;
            },
            a, b, boolVector);
        return sp;
    }

    template<typename V, typename B>
    auto dot(const V& a, const V& b)
    {
        double sp(0);
        applyToVector(
            [&sp](double a, double b)
            {
                sp += a*b;
            },
            a, b);
        return sp;
    }

    template<typename MI1, typename MI2, std::size_t L>
    bool cmp(const MI1 & idx1, const MI2 & idx2,
        index_constant<L> level)
    {
        bool same = true;
        for (std::size_t i=0; i<L; i++)
            same &= (idx1[i] == idx2[i]);
        return same;
    }

    template<typename MI1, typename MI2>
    bool cmp(const MI1 & idx1, const MI2 & idx2)
    {
        bool same = true;
        std::size_t level = std::min(idx1.size(), idx2.size());
        for (std::size_t i=0; i<level; i++)
            same &= (idx1[i] == idx2[i]);
        return same;
    }

    template<typename V, typename It, std::size_t max_level, std::size_t L>
    double SP_skip (const V & x, const V & y,
        It& w, const It wend,
        std::array<int, max_level> & indices,
        bool critical,
        index_constant<L> level)
    {
        // TODO support static loop for MultiType vectors

        // are we at the leaf?
        static const bool isScalar =
            Impl::IsScalar<std::decay_t<V>>::value;

        if constexpr (isScalar)
        {
            return x*y;
        }
        else
        {
            // std::cout << "l" << level << "\n";
            double sp = 0;
            for (indices[level]=0; indices[level]<x.size(); indices[level]++)
            {
                bool doCheck = (w!=wend && critical && (*w)[level] == indices[level]);
                if (doCheck && (w->size() == L+1) && cmp(*w,indices,level)) // check AND end of skip-index AND indices match
                {
                    w++;
                }
                else {
                    if constexpr(ISTL::hasDynamicIndexAccess<V>{}) {
                        // access block directly
                        int i = indices[level];
                        sp += SP_skip(x[i], y[i], w, wend, indices, doCheck, index_constant<L+1>{});
                    }
                    else if constexpr(ISTL::hasStaticIndexAccess<V>{}) {
                        // access block by searching for the static index equal to mi[i]
                        using Size = decltype(Hybrid::size(x));
                        Hybrid::switchCases(std::make_index_sequence<Size::value>{},indices[level],
                            [&](auto i) {
                                sp += SP_skip(x[i], y[i], w, wend, indices, doCheck, index_constant<L+1>{}); });
                    }
                }
            }
            return sp;
        }
    }

}} // end namespace Dune::ScalarProductHelper

using namespace Dune;

/*
  Calculate norm of
  ((1,2),(2,3),(3,4),(4,5),(5,6))

  with an inner weight

  ((1,1),(1,1),(0,0),(0,1),(1,0))

  i.e. skipping entries

  (2,0),(2,1),(3,0),(4,1)

  or with a dynamic index: (2),(3,0),(4,1)

  The result of the weighted sp x * w * x = 68
 */

std::vector<std::array<int,2>> skipEntries {{2,0},{2,1},{3,0},{4,1}};
std::vector<std::array<double,2>> x {{1,2},{2,3},{3,4},{4,5},{5,6}};

template<typename V, typename B>
double sp_A (const V & x, const V & y, const B & boolVec) {
    return Dune::ScalarProductHelper::dot(x, y, boolVec);
}

template<typename V, typename S>
auto sp_B(const V & x, const V & y, const S & skipEntries)
{
    static const unsigned int max_level = 3;
    std::array<int, max_level> indices {0,0,0};
    auto it = skipEntries.begin();
    auto itend = skipEntries.end();
    return Dune::ScalarProductHelper::SP_skip(x, y, it, itend, indices, true,
        index_constant<0>{});
}

template<typename V, typename S>
auto sp_C(const V & x, const V & y, const S & skipEntries)
{
    double sp = x.dot(y); // dot(x,y);
    double skip = 0;
    ISTL::forEachIndex (skipEntries, y,
        [&skip](auto & v, auto & mi) { skip += dot(v,v); }
        );

    return sp-skip;
}

/*
  Calculate norm of
  ((1,2),(2,3),(3,4),(4,5),(5,6))

  with an inner weight

  ((1,1),(1,1),(0,0),(0,1),(1,0))

  i.e. skipping entries

  (2,0),(2,1),(3,0),(4,1)

  or with a dynamic index: (2),(3,0),(4,1)

  The result of the weighted sp x * w * x = 68
 */
auto generate_basic(bool unbalanced)
{
    std::vector<ReservedVector<int,2>> skipIdx {{2,0},{2,1},{3,0},{4,1}};
    std::vector<ReservedVector<int,2>> skipIdx2 {{2},{3,0},{4,1}};
    BlockVector<FieldVector<double,2>> x {{1,2},{2,3},{3,4},{4,5},{5,6}};
    BlockVector<FieldVector<bool,2>> useBool {{1,1},{1,1},{0,0},{0,1},{1,0}};
    if (unbalanced)
        return std::make_tuple(x, useBool, skipIdx2);
    else
        return std::make_tuple(x, useBool, skipIdx);
}

template<std::size_t B, std::size_t maxSize = 2>
auto generate_nested(int N,
    double propability,
    double propability2
    )
{
    // generate large vector
    BlockVector<FieldVector<double,B>> x(N);
    for (unsigned int i=0; i<N; i++)
        for (unsigned int j=0; j<B; j++)
            x[i][j] = drand48();

    // generate bool vector
    // and skip entries from bool vector
    BlockVector<FieldVector<bool,B>> useBool(N);
    using MultiIndex = Dune::ReservedVector<unsigned int,maxSize>;
    std::vector<MultiIndex> skipIdx;
    skipIdx.reserve(2*propability*N*B);
    for (unsigned int i=0; i<N; i++)
    {
        bool might_skip_block = (drand48()<propability)?true:false;
        bool skip_block = true;
        for (unsigned int j=0; j<B; j++) {
            bool skip_entry = drand48()>(1-propability2)?true:false;
            useBool[i][j] = not(might_skip_block && skip_entry);
            skip_block = skip_block && (not useBool[i][j]);
        }
        if (skip_block)
        {
            MultiIndex idx({i});
            skipIdx.push_back(idx);
        }
        else
        {
            for (unsigned int j=0; j<B; j++) {
                if (useBool[i][j] == false) {
                    MultiIndex idx({i,j});
                    skipIdx.push_back(idx);
                }
            }
        }
    }

    // print skip and bool vector
    if constexpr (0) {
        for (auto && _B : useBool)
        {
            std::cout << ": ";
            for (auto && _b : _B)
                std::cout << _b << "\t";
            std::cout << "\n";
        }

        for (auto && idx : skipIdx)
            std::cout << "- " << idx << std::endl;
    }

    return std::make_tuple(x, useBool, skipIdx);
}

template<std::size_t maxSize = 1>
auto generate_flat(int N,
    double propability)
{
    // generate large vector
    BlockVector<double> x(N);
    for (unsigned int i=0; i<N; i++)
            x[i] = drand48();

    // generate bool vector
    // and skip entries from bool vector
    BlockVector<unsigned char> useBool(N);
    using MultiIndex = Dune::ReservedVector<unsigned int,maxSize>;
    std::vector<MultiIndex> skipIdx;
    skipIdx.reserve(2*propability*N);
    for (unsigned int i=0; i<N; i++)
    {
        bool skip_entry = (drand48()<propability)?true:false;
        useBool[i] = not(skip_entry);
        if (skip_entry)
        {
            MultiIndex idx({i});
            skipIdx.push_back(idx);
        }
    }

    return std::make_tuple(x, useBool, skipIdx);
}

auto generate_multitype(int N1, int N2,
    double propability = 0.01
    )
{
    // generate a tailor-hood like structure
    auto [x1, useBool1, skipIdx1] = generate_flat(N1,propability);
    auto [x2, useBool2, skipIdx2] = generate_nested<3>(N2,propability,1);

    // combine blocks
    MultiTypeBlockVector<
        BlockVector<double>,
        BlockVector<FieldVector<double,3>>
        //decltype(x1),decltype(x2)
        > x(x1,x2);
    MultiTypeBlockVector<decltype(useBool1),decltype(useBool2)> useBool(useBool1,useBool2);

    // merge index lists
    using MultiIndex = Dune::ReservedVector<unsigned int,3>;
    std::vector<MultiIndex> skipIdx;
    skipIdx.reserve(skipIdx1.size()+skipIdx2.size());
    for (auto i : skipIdx1)
    {
        MultiIndex idx({0});
        for (auto j : i)
            idx.push_back(j);
        skipIdx.push_back(idx);
    }
    for (auto i : skipIdx2)
    {
        MultiIndex idx({1});
        for (auto j : i)
            idx.push_back(j);
        skipIdx.push_back(idx);
    }

    return std::make_tuple(x, useBool, skipIdx);
}

template<typename F, typename... Args>
void do_test(std::string name, F && f, Args&&... args)
{
    auto [x, useBool, skipIdx] = f(args...);

    // statistics
    std::size_t skipped = 0;
    std::size_t sz =
        flatVectorForEach(useBool,
            [&skipped](auto&& entry, std::size_t offset) { skipped += (entry==false)?1:0; });

    // run performance tests
    Dune::Timer t;
    t.start();
    auto xA = sp_A(x,x,useBool);
    double tA = t.elapsed();
    t.reset();
    auto xB = sp_B(x,x,skipIdx);
    double tB = t.elapsed();
    t.reset();
    auto xC = sp_C(x,x,skipIdx);
    double tC = t.elapsed();
    t.stop();

    // report
    std::cout << "----------------------------------------------\n" << name << std::endl;
    std::cout << "Stats: skip " << skipped << " of " << sz << " entries "
              << "(" << 100 * skipped / sz << "%), "
              << " marked " << skipIdx.size() << " Indices\n";
    std::cout << "\talgo A\t" << xA << "\t(" << tA << " s)" << std::endl;
    std::cout << "\talgo B\t" << xB << "\t(" << tB << " s)"
              << " -> diff " << std::abs(xA-xB) << std::endl;
    std::cout << "\talgo C\t" << xC << "\t(" << tC << " s)"
              << " -> diff " << std::abs(xA-xC) << std::endl;

    // check result
    if ( std::abs(xA-xB)/xA > 1e-12 )
        std::cerr << "ERROR: result of algoB differs from reference"
                  << " (rel error " << std::abs(xA-xB)/xA << ")" << std::endl;
    if ( std::abs(xA-xC)/xA > 1e-12 )
        std::cerr << "ERROR: result of algoC differs from reference"
                  << " (rel error " << std::abs(xA-xC)/xA << ")" << std::endl;
}

#define test(F,INFO,...) do_test(std::string(#F) + INFO, generate_ ## F, __VA_ARGS__)

int main()
{
    using namespace Dune::Indices;
    test(basic," [balanced]",false);
    test(basic," [unbalanced]",true);
    test(nested<5> ," [simple]", 10, 0.1, 0.9);
    test(nested<5> ," [larger]", 10000000, 0.01, 0.95);
    test(nested<5> ," [many skipped]", 10000000, 0.2, 0.5);
    test(nested<10>," [few skipped]", 10000000,0.01,0.5);
    test(nested<10>," [only blocks]", 10000000,0.01,1.0);
    test(nested<10>," [20Mx10]", 20000000,0.01,1.0);
    test(nested<2> ," [100Mx2]", 100000000,0.01,1.0);
    test(nested<1> ," [200Mx1]", 200000000,0.01,1.0);
    test(nested<10> ," [long MI]", 20000000,0.01,0.5);
    test(multitype," [2/3 levels]", 2000000,8000000,0.01);
    std::cout << "----------------------------------------------\n";
}
