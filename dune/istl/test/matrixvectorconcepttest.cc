#if HAVE_CONFIG_H
#include <config.h>
#endif

#if __has_include(<concepts>)

#include <array>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/concepts/matrix.hh>
#include <dune/common/concepts/linearmap.hh>
#include <dune/common/concepts/vector.hh>
#include <dune/common/concepts/vectorspace.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/scaledidmatrix.hh>

int main(int argc, char** argv)
{
  using M = Dune::FieldMatrix<double,2,2>;
  using Mat = Dune::BCRSMatrix<M>;
  static_assert(Dune::Concept::Matrix<Mat>);
  static_assert(Dune::Concept::MutableMatrix<Mat>);
  static_assert(Dune::Concept::TraversableMatrix<Mat>);
  static_assert(Dune::Concept::ResizeableMatrix<Mat>);
  static_assert(Dune::Concept::VectorSpace<Mat>);

  using V = Dune::FieldVector<double,2>;
  using Vec = Dune::BlockVector<V>;
  static_assert(Dune::Concept::Vector<Vec>);
  static_assert(Dune::Concept::VectorSpace<Vec>);
  static_assert(Dune::Concept::LinearMap<Mat, Vec, Vec>);
  static_assert(Dune::Concept::TransposableLinearMap<Mat, Vec, Vec>);
  static_assert(Dune::Concept::HermitianLinearMap<Mat, Vec, Vec>);

  using Mat2 = Dune::ScaledIdentityMatrix<double,2>;
  static_assert(Dune::Concept::Matrix<Mat2>);
  static_assert(Dune::Concept::MutableMatrix<Mat2>);
  static_assert(Dune::Concept::TraversableMatrix<Mat2>);
  static_assert(Dune::Concept::ResizeableMatrix<Mat2>);
  static_assert(Dune::Concept::VectorSpace<Mat2>);
  static_assert(Dune::Concept::LinearMap<Mat2, V, V>);
  static_assert(Dune::Concept::TransposableLinearMap<Mat2, V, V>);
  static_assert(Dune::Concept::HermitianLinearMap<Mat2, V, V>);

  using Mat3 = Dune::Matrix<M>;
  static_assert(Dune::Concept::Matrix<Mat3>);
  static_assert(Dune::Concept::MutableMatrix<Mat3>);
  static_assert(Dune::Concept::TraversableMatrix<Mat3>);
  static_assert(Dune::Concept::ResizeableMatrix<Mat3>);
  static_assert(Dune::Concept::VectorSpace<Mat3>);
  static_assert(Dune::Concept::LinearMap<Mat3, Vec, Vec>);
  static_assert(Dune::Concept::TransposableLinearMap<Mat3, Vec, Vec>);
  static_assert(Dune::Concept::HermitianLinearMap<Mat3, Vec, Vec>);
}

#else
int main()
{
  return 77;
}
#endif
