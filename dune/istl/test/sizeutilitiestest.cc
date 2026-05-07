// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/sizeutilities.hh>

using namespace Dune;

// Import the static constants _0, _1, etc
using namespace Indices;


int main(int argc, char** argv) try
{
  Dune::TestSuite suite;

  // test vectors...

  using Vector1 = Dune::MultiTypeBlockVector<Dune::FieldVector<double,1>, Dune::FieldVector<double,1>>;
  Vector1 vector1;
  auto nEntries1 = Dune::Hybrid::numEntries(vector1);
  static_assert(nEntries1 == Vector1::size());
  Dune::Hybrid::forEach(Dune::Hybrid::entries(vector1), [&](auto i) {
    vector1[i] = 1.0;
  });

  using Vector2 = Dune::FieldVector<double,4>;
  Vector2 vector2;
  auto nEntries2 = Dune::Hybrid::numEntries(vector2);
  static_assert(std::is_same_v<decltype(nEntries2), std::size_t>);
  suite.check(nEntries2 == vector2.size());
  Dune::Hybrid::forEach(Dune::Hybrid::entries(vector2), [&](auto i) {
    vector2[i] = 1.0;
  });

  using Vector3 = Dune::BlockVector<Dune::FieldVector<double,1>>;
  Vector3 vector3(7);
  auto nEntries3 = Dune::Hybrid::numEntries(vector3);
  static_assert(std::is_same_v<decltype(nEntries3), std::size_t>);
  suite.check(nEntries3 == vector3.size());
  Dune::Hybrid::forEach(Dune::Hybrid::entries(vector3), [&](auto i) {
    vector3[i] = 1.0;
  });

  // test matrices...

  using Matrix1 = Dune::MultiTypeBlockMatrix<
    Dune::MultiTypeBlockVector<Dune::FieldMatrix<double,1,1>, Dune::FieldMatrix<double,1,1>>,
    Dune::MultiTypeBlockVector<Dune::FieldMatrix<double,1,1>, Dune::FieldMatrix<double,1,1>> >;
  Matrix1 matrix1;
  auto nRows1 = Dune::Hybrid::numRows(matrix1);
  auto nCols1 = Dune::Hybrid::numCols(matrix1);
  static_assert(nRows1 == Matrix1::N());
  static_assert(nCols1 == Matrix1::M());
  Dune::Hybrid::forEach(Dune::Hybrid::rows(matrix1), [&](auto i) {
    Dune::Hybrid::forEach(Dune::Hybrid::cols(matrix1), [&](auto j) {
      matrix1[i][j] = 1.0;
    });
  });

  using Matrix2 = Dune::FieldMatrix<double,3,3>;
  Matrix2 matrix2;
  auto nRows2 = Dune::Hybrid::numRows(matrix2);
  auto nCols2 = Dune::Hybrid::numCols(matrix2);
  static_assert(std::is_same_v<decltype(nRows2), std::size_t>);
  static_assert(std::is_same_v<decltype(nCols2), std::size_t>);
  suite.check(nRows2 == matrix2.N());
  suite.check(nCols2 == matrix2.M());
  Dune::Hybrid::forEach(Dune::Hybrid::rows(matrix2), [&](auto i) {
    Dune::Hybrid::forEach(Dune::Hybrid::cols(matrix2), [&](auto j) {
      matrix2[i][j] = 1.0;
    });
  });

  using Matrix3 = Dune::Matrix<Dune::FieldMatrix<double,1,1>>;
  Matrix3 matrix3(9,9);
  auto nRows3 = Dune::Hybrid::numRows(matrix3);
  auto nCols3 = Dune::Hybrid::numCols(matrix3);
  static_assert(std::is_same_v<decltype(nRows3), std::size_t>);
  static_assert(std::is_same_v<decltype(nCols3), std::size_t>);
  suite.check(nRows3 == matrix3.N());
  suite.check(nCols3 == matrix3.M());
  Dune::Hybrid::forEach(Dune::Hybrid::rows(matrix3), [&](auto i) {
    Dune::Hybrid::forEach(Dune::Hybrid::cols(matrix3), [&](auto j) {
      matrix3[i][j] = 1.0;
    });
  });

  return suite.exit();
}
catch (Dune::Exception& e)
{
  std::cerr << "DUNE reported an exception: " << e << std::endl;
  return 1;
}
catch (std::exception& e)
{
  std::cerr << "C++ reported an exception: " << e.what() << std::endl;
  return 2;
} catch (...)
{
  std::cerr << "Unknown exception encountered!" << std::endl;
  return 3;
}
