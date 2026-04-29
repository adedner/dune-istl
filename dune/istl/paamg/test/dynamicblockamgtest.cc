// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Template for dense (DynamicMatrix) blocks
template<typename MatrixBlockType, typename VectorBlockType>
struct BlockTraits;

// Specialization for DynamicMatrix dense blocks
template<>
struct BlockTraits<Dune::DynamicMatrix<double>, Dune::DynamicVector<double>>
{
  static constexpr const char* name() { return "DynamicMatrix"; }

  static void setRectangularDiagonal(Dune::DynamicMatrix<double>& block, std::size_t rows, std::size_t cols, double value)
  {
    block.resize(rows, cols);
    block = 0.0;

    const std::size_t diagSize = std::min(rows, cols);
    for (std::size_t k = 0; k < diagSize; ++k)
      block[k][k] = value;
  }
};

// Specialization for BCRSMatrix sparse blocks
template<>
struct BlockTraits<Dune::BCRSMatrix<double>, Dune::BlockVector<double>>
{
  static constexpr const char* name() { return "BCRSMatrix"; }

  static void setRectangularDiagonal(Dune::BCRSMatrix<double>& block, std::size_t rows, std::size_t cols, double value)
  {
    const std::size_t diagSize = std::min(rows, cols);
    block.setBuildMode(Dune::BCRSMatrix<double>::row_wise);
    block.setSize(rows, cols, diagSize*3);

    auto iter = block.createbegin();
    for (std::size_t i = 0; i < rows; ++i, ++iter) {
      if (i < diagSize)
        iter.insert(i);
      if (i + 2 < cols)
        iter.insert(i + 2);
      if (i >= 2 && i - 2 < cols)
        iter.insert(i - 2);
    }

    for (std::size_t k = 0; k < diagSize; ++k)
      block[k][k] = value;
  }
};

template<typename MatrixBlockType, typename VectorBlockType>
int nodeIndex(int x, int y, int N)
{
  return y * N + x;
}

template<typename MatrixBlockType, typename VectorBlockType>
bool isBoundary(int x, int y, int N)
{
  return (x == 0 || y == 0 || x == N - 1 || y == N - 1);
}

template<typename MatrixBlockType, typename VectorBlockType>
using Matrix_t = Dune::BCRSMatrix<MatrixBlockType>;

template<typename MatrixBlockType, typename VectorBlockType>
using Vector_t = Dune::BlockVector<VectorBlockType>;

template<typename MatrixBlockType, typename VectorBlockType>
Matrix_t<MatrixBlockType, VectorBlockType> setupDynamicBlockAnisotropic2d(int N, std::vector<std::size_t>& blockSizes, double eps)
{
  const int nodes = N * N;
  blockSizes.resize(nodes);

  for (int y = 0; y < N; ++y) {
    for (int x = 0; x < N; ++x) {
      const int row = nodeIndex<MatrixBlockType, VectorBlockType>(x, y, N);
      blockSizes[row] = 1 + static_cast<std::size_t>((x + 2 * y) % 5);
    }
  }

  Dune::MatrixIndexSet pattern(nodes, nodes);
  for (int y = 0; y < N; ++y) {
    for (int x = 0; x < N; ++x) {
      const int row = nodeIndex<MatrixBlockType, VectorBlockType>(x, y, N);

      pattern.add(row, row);
      if (isBoundary<MatrixBlockType, VectorBlockType>(x, y, N))
        continue;

      pattern.add(row, nodeIndex<MatrixBlockType, VectorBlockType>(x - 1, y, N));
      pattern.add(row, nodeIndex<MatrixBlockType, VectorBlockType>(x + 1, y, N));
      pattern.add(row, nodeIndex<MatrixBlockType, VectorBlockType>(x, y - 1, N));
      pattern.add(row, nodeIndex<MatrixBlockType, VectorBlockType>(x, y + 1, N));
    }
  }

  Matrix_t<MatrixBlockType, VectorBlockType> matrix;
  pattern.exportIdx(matrix);

  for (int y = 0; y < N; ++y) {
    for (int x = 0; x < N; ++x) {
      const int row = nodeIndex<MatrixBlockType, VectorBlockType>(x, y, N);
      const std::size_t rowSize = blockSizes[row];

      for (auto col = matrix[row].begin(); col != matrix[row].end(); ++col) {
        const int colIdx = col.index();
        const std::size_t colSize = blockSizes[colIdx];

        if (isBoundary<MatrixBlockType, VectorBlockType>(x, y, N)) {
          if (colIdx == row)
            BlockTraits<MatrixBlockType, VectorBlockType>::setRectangularDiagonal(*col, rowSize, colSize, 1.0);
          else
            BlockTraits<MatrixBlockType, VectorBlockType>::setRectangularDiagonal(*col, rowSize, colSize, 0.0);
          continue;
        }

        if (colIdx == row)
          BlockTraits<MatrixBlockType, VectorBlockType>::setRectangularDiagonal(*col, rowSize, colSize, 2.0 + 2.0 * eps);
        else if (colIdx == nodeIndex<MatrixBlockType, VectorBlockType>(x - 1, y, N) ||
                 colIdx == nodeIndex<MatrixBlockType, VectorBlockType>(x + 1, y, N))
          BlockTraits<MatrixBlockType, VectorBlockType>::setRectangularDiagonal(*col, rowSize, colSize, -eps);
        else
          BlockTraits<MatrixBlockType, VectorBlockType>::setRectangularDiagonal(*col, rowSize, colSize, -1.0);
      }
    }
  }

  return matrix;
}

template<typename MatrixBlockType, typename VectorBlockType>
void verifyRectangularOffDiagonalBlocks(const Matrix_t<MatrixBlockType, VectorBlockType>& matrix, const std::vector<std::size_t>& blockSizes)
{
  bool sawRectangularOffDiagonal = false;

  for (std::size_t row = 0; row < matrix.N(); ++row) {
    for (auto col = matrix[row].begin(); col != matrix[row].end(); ++col) {
      const std::size_t colIdx = col.index();
      const std::size_t expectedRows = blockSizes[row];
      const std::size_t expectedCols = blockSizes[colIdx];

      if (col->N() != expectedRows || col->M() != expectedCols)
        throw std::runtime_error("dynamic block has unexpected shape");

      if (colIdx != row && expectedRows != expectedCols)
        sawRectangularOffDiagonal = true;
    }
  }

  if (!sawRectangularOffDiagonal)
    throw std::runtime_error("did not encounter a rectangular off-diagonal block");
}

template<typename MatrixBlockType, typename VectorBlockType>
void resizeAndAssign(Vector_t<MatrixBlockType, VectorBlockType>& vector, const std::vector<std::size_t>& blockSizes, double value)
{
  vector.resize(blockSizes.size());
  for (std::size_t i = 0; i < blockSizes.size(); ++i) {
    vector[i].resize(blockSizes[i]);
    vector[i] = value;
  }
}

template<typename MatrixBlockType, typename VectorBlockType>
void testDynamicBlockAMG(int N, int coarsenTarget, int maxLevel)
{
  std::vector<std::size_t> blockSizes;
  Matrix_t<MatrixBlockType, VectorBlockType> matrix = setupDynamicBlockAnisotropic2d<MatrixBlockType, VectorBlockType>(N, blockSizes, 1.0);
  verifyRectangularOffDiagonalBlocks<MatrixBlockType, VectorBlockType>(matrix, blockSizes);

  Vector_t<MatrixBlockType, VectorBlockType> x;
  Vector_t<MatrixBlockType, VectorBlockType> b;
  Vector_t<MatrixBlockType, VectorBlockType> ones;
  resizeAndAssign<MatrixBlockType, VectorBlockType>(x, blockSizes, 0.0);
  resizeAndAssign<MatrixBlockType, VectorBlockType>(b, blockSizes, 0.0);
  resizeAndAssign<MatrixBlockType, VectorBlockType>(ones, blockSizes, 1.0);

  matrix.mv(ones, b);

  // Write matrix structure to SVG
  std::string svgFilename = std::string("dynamicblockamgtest_") + BlockTraits<MatrixBlockType, VectorBlockType>::name() + ".svg";
  std::ofstream svgFile(svgFilename);
  Dune::writeSVGMatrix(svgFile, matrix);
  svgFile.close();
  std::cout << "Matrix structure written to " << svgFilename << std::endl;

  using Matrix = Matrix_t<MatrixBlockType, VectorBlockType>;
  using Vector = Vector_t<MatrixBlockType, VectorBlockType>;
  using Operator = Dune::MatrixAdapter<Matrix, Vector, Vector>;

  Operator op(matrix);

  using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::UnSymmetricCriterion<Matrix, Dune::Amg::FrobeniusNorm>>;
  using Smoother = Dune::SeqJac<Matrix, Vector, Vector, 2>;
  using AMG = Dune::Amg::AMG<Operator, Vector, Smoother>;

  using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
  SmootherArgs smootherArgs;
  smootherArgs.iterations = 1;
  smootherArgs.relaxationFactor = 1.0;

  Criterion criterion(15, coarsenTarget);
  criterion.setDefaultValuesIsotropic(2);
  criterion.setAlpha(0.67);
  criterion.setBeta(1.0e-4);
  criterion.setMaxLevel(maxLevel);
  criterion.setSkipIsolated(false);
  criterion.setNoPreSmoothSteps(1);
  criterion.setNoPostSmoothSteps(1);

  AMG amg(op, criterion, smootherArgs);

  Dune::BiCGSTABSolver<Vector> solver(op, amg, 1e-8, 250, 2);
  Dune::InverseOperatorResult result;
  solver.apply(x, b, result);

  if (!result.converged)
    throw std::runtime_error("AMG with dynamic rectangular blocks did not converge");

  std::cout << "dynamicblockamgtest (" << BlockTraits<MatrixBlockType, VectorBlockType>::name()
            << ") converged in " << result.iterations
            << " iterations with reduction " << result.reduction << std::endl;
}

} // namespace

int main(int argc, char** argv)
try
{
  int N = 24;
  int coarsenTarget = 800;
  int maxLevel = 10;

  if (argc > 1)
    N = std::atoi(argv[1]);

  if (argc > 2)
    coarsenTarget = std::atoi(argv[2]);

  if (argc > 3)
    maxLevel = std::atoi(argv[3]);

  std::cout << "Testing AMG with DynamicMatrix dense blocks..." << std::endl;
  testDynamicBlockAMG<Dune::DynamicMatrix<double>, Dune::DynamicVector<double>>(N, coarsenTarget, maxLevel);

  std::cout << "Testing AMG with BCRSMatrix sparse blocks..." << std::endl;
  testDynamicBlockAMG<Dune::BCRSMatrix<double>, Dune::BlockVector<double>>(N, coarsenTarget, maxLevel);

  return 0;
}
catch (const std::exception& e)
{
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 1;
}
catch (...)
{
  std::cerr << "Dune reported an unknown error." << std::endl;
  return 1;
}
