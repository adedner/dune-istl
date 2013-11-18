// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#undef NDEBUG // make sure assert works

#include <dune/common/float_cmp.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/istl/bcrsmatrix.hh>

typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > ScalarMatrix;

void buildMatrix(ScalarMatrix& m)
{
  m.entry(0,0) = 1.0; m.entry(0,1) = 1.0; m.entry(0,2) = 1.0;  m.entry(0,3) = 1.0;
  m.entry(1,0) = 1.0; m.entry(1,1) = 1.0; m.entry(1,2) = 1.0;
  m.entry(2,1) = 1.0; m.entry(2,2) = 1.0; m.entry(2,3) = 1.0;
  m.entry(3,2) = 1.0; m.entry(3,3) = 1.0; m.entry(3,4) = 1.0;
  m.entry(4,3) = 1.0; m.entry(4,4) = 1.0; m.entry(4,5) = 1.0;
  m.entry(5,4) = 1.0; m.entry(5,5) = 1.0; m.entry(5,6) = 1.0;
  m.entry(6,5) = 1.0; m.entry(6,6) = 1.0; m.entry(6,7) = 1.0;
  m.entry(7,6) = 1.0; m.entry(7,7) = 1.0; m.entry(7,8) = 1.0;
  m.entry(8,7) = 1.0; m.entry(8,8) = 1.0; m.entry(8,9) = 1.0;
  m.entry(9,8) = 1.0; m.entry(9,9) = 1.0;
  // add some more entries in random order
  m.entry(7,3) = 1.0;
  m.entry(6,0) = 1.0;
  m.entry(3,8) = 1.0;
}

template<typename M>
void setMatrix(M& m)
{
  m[0][0] = 1.0; m[0][1] = 1.0; m[0][2] = 1.0;  m[0][3] = 1.0;
  m[1][0] = 1.0; m[1][1] = 1.0; m[1][2] = 1.0;
  m[2][1] = 1.0; m[2][2] = 1.0; m[2][3] = 1.0;
  m[3][2] = 1.0; m[3][3] = 1.0; m[3][4] = 1.0;
  m[4][3] = 1.0; m[4][4] = 1.0; m[4][5] = 1.0;
  m[5][4] = 1.0; m[5][5] = 1.0; m[5][6] = 1.0;
  m[6][5] = 1.0; m[6][6] = 1.0; m[6][7] = 1.0;
  m[7][6] = 1.0; m[7][7] = 1.0; m[7][8] = 1.0;
  m[8][7] = 1.0; m[8][8] = 1.0; m[8][9] = 1.0;
  m[9][8] = 1.0; m[9][9] = 1.0;
  // add some more entries in random order
  m[7][3] = 1.0;
  m[6][0] = 1.0;
  m[3][8] = 1.0;
}

void testImplicitBuild()
{
  ScalarMatrix m(10,10,3,0.1,ScalarMatrix::implicit);
  buildMatrix(m);
  ScalarMatrix::CompressionStatistics stats = m.compress();
  assert(Dune::FloatCmp::eq(stats.avg,33./10.));
  assert(stats.maximum == 4);
  assert(stats.overflow_total == 4);
  setMatrix(m);
}

void testImplicitBuildWithInsufficientOverflow()
{
  try {
    ScalarMatrix m(10,10,1,0,ScalarMatrix::implicit);
    // add diagonal entries + completely fill the first row with entries
    // with the current base buffer of 4 * avg, that should be enough to make
    // compress fail.
    for (int i = 0; i < 10; ++i)
      {
        m.entry(i,i) = 1.0;
        m.entry(0,i) = 1.0;
      }
    m.compress();
    assert(false && "compress() should have thrown an exception");
  } catch (Dune::ImplicitModeOverflowExhausted& e) {
    // test passed
  }
}

void testSetterInterface()
{
  ScalarMatrix m;
  m.setBuildMode(ScalarMatrix::implicit);
  m.setImplicitBuildModeParameters(3,0.1);
  m.setSize(10,10);
  buildMatrix(m);
  ScalarMatrix::CompressionStatistics stats = m.compress();
  assert(Dune::FloatCmp::eq(stats.avg,33.0/10.0));
  assert(stats.maximum == 4);
  assert(stats.overflow_total == 4);
}

void testDoubleSetSize()
{
  ScalarMatrix m;
  m.setBuildMode(ScalarMatrix::implicit);
  m.setImplicitBuildModeParameters(3,0.1);
  m.setSize(14,14);
  m.setSize(10,10);
  buildMatrix(m);
  ScalarMatrix::CompressionStatistics stats = m.compress();
  assert(Dune::FloatCmp::eq(stats.avg,33.0/10.0));
  assert(stats.maximum == 4);
  assert(stats.overflow_total == 4);
}

void testInvalidBuildModeConstructorCall()
{
  try {
    ScalarMatrix m(10,10,1,-1.0,ScalarMatrix::random);
    assert(false && "Constructor should have thrown an exception!");
  } catch (Dune::BCRSMatrixError& e) {
    // test passed
  }
}

void testNegativeOverflowConstructorCall()
{
  try {
    ScalarMatrix m(10,10,1,-1.0,ScalarMatrix::implicit);
    assert(false && "Constructor should have thrown an exception!");
  } catch (Dune::BCRSMatrixError& e) {
    // test passed
  }
}

void testInvalidSetImplicitBuildModeParameters()
{
  try {
    ScalarMatrix m;
    m.setBuildMode(ScalarMatrix::implicit);
    m.setImplicitBuildModeParameters(1,-1.0);
    assert(false && "setImplicitBuildModeParameters() should have thrown an exception!");
  } catch (Dune::BCRSMatrixError& e) {
    // test passed
  }
}

void testSetImplicitBuildModeParametersAfterSetSize()
{
  try {
    ScalarMatrix m;
    m.setBuildMode(ScalarMatrix::implicit);
    m.setImplicitBuildModeParameters(3,0.1);
    m.setSize(10,10);
    m.setImplicitBuildModeParameters(4,0.1);
    assert(false && "setImplicitBuildModeParameters() should have thrown an exception!");
  } catch (Dune::InvalidStateException& e) {
    // test passed
  }
}

void testSetSizeWithNonzeroes()
{
  try {
    ScalarMatrix m;
    m.setBuildMode(ScalarMatrix::implicit);
    m.setImplicitBuildModeParameters(3,0.1);
    m.setSize(10,10,300);
    assert(false && "setSize() should have thrown an exception!");
  } catch (Dune::BCRSMatrixError& e) {
    // test passed
  }
}

void testCopyConstructionAndAssignment()
{
  ScalarMatrix m(10,10,3,0.1,ScalarMatrix::implicit);
  buildMatrix(m);
  m.compress();
  ScalarMatrix m2(m);
  m2 = 3.0;
  ScalarMatrix m3(m);
  m3 = m2;
  ScalarMatrix m4;
  m4 = m;
}

void testInvalidCopyConstruction()
{
  try {
    ScalarMatrix m(10,10,3,0.1,ScalarMatrix::implicit);
    buildMatrix(m);
    ScalarMatrix m2(m);
    assert(false && "copy constructor should have thrown an exception!");
  } catch (Dune::InvalidStateException& e) {
    // test passed
  }
}

void testInvalidCopyAssignment()
{
  ScalarMatrix m(10,10,3,0.1,ScalarMatrix::implicit);
  buildMatrix(m);
  // copy incomplete matrix into empty one
  try {
    ScalarMatrix m2;
    assert(false && "operator=() should have thrown an exception!");
  } catch (Dune::InvalidStateException& e) {
    // test passed
  }
  // copy incomplete matrix into full one
  try {
    ScalarMatrix m2(10,10,3,0.1,ScalarMatrix::implicit);
    buildMatrix(m2);
    m2.compress();
    m2 = m;
    assert(false && "operator=() should have thrown an exception!");
  } catch (Dune::InvalidStateException& e) {
    // test passed
  }
  // copy fully build matrix into half-built one
  m.compress();
  try {
    ScalarMatrix m2(10,10,3,0.1,ScalarMatrix::implicit);
    buildMatrix(m2);
    m2 = m;
    assert(false && "operator=() should have thrown an exception!");
  } catch (Dune::InvalidStateException& e) {
    // test passed
  }
}

void testEntryConsistency()
{
  ScalarMatrix m(10,10,3,0.1,ScalarMatrix::implicit);
  assert(Dune::FloatCmp::eq(static_cast<const double&>(m.entry(0,3)),0.0));
  assert(Dune::FloatCmp::eq(static_cast<const double&>(m.entry(7,6)),0.0));
  buildMatrix(m);
  assert(Dune::FloatCmp::eq(static_cast<const double&>(m.entry(0,3)),1.0));
  assert(Dune::FloatCmp::eq(static_cast<const double&>(m.entry(7,6)),1.0));
  m.entry(4,4) += 3.0;
  assert(Dune::FloatCmp::eq(static_cast<const double&>(m.entry(4,4)),4.0));
  m.compress();
  assert(Dune::FloatCmp::eq(static_cast<const double&>(m[0][3]),1.0));
  assert(Dune::FloatCmp::eq(static_cast<const double&>(m[7][6]),1.0));
  assert(Dune::FloatCmp::eq(static_cast<const double&>(m[4][4]),4.0));
}

void testEntryAfterCompress()
{
  try {
    ScalarMatrix m(10,10,3,0.1,ScalarMatrix::implicit);
    buildMatrix(m);
    m.compress();
    m.entry(3,3);
    assert(false && "entry() should have thrown an exception!");
  } catch (Dune::BCRSMatrixError& e) {
    // test passed
  }
}

void testBracketOperatorBeforeCompress()
{
  try {
    ScalarMatrix m(10,10,3,0.1,ScalarMatrix::implicit);
    buildMatrix(m);
    m[3][3];
    assert(false && "operator[]() should have thrown an exception!");
  } catch (Dune::BCRSMatrixError& e) {
    // test passed
  }
}

void testConstBracketOperatorBeforeCompress()
{
  try {
    ScalarMatrix m(10,10,3,0.1,ScalarMatrix::implicit);
    buildMatrix(m);
    const_cast<const ScalarMatrix&>(m)[3][3];
    assert(false && "operator[]() should have thrown an exception!");
  } catch (Dune::BCRSMatrixError& e) {
    // test passed
  }
}

void testImplicitMatrixBuilder()
{
  ScalarMatrix m(10,10,3,0.1,ScalarMatrix::implicit);
  Dune::ImplicitMatrixBuilder<ScalarMatrix> b(m);
  setMatrix(b);
  m.compress();
  setMatrix(m);
}

void testImplicitMatrixBuilderExtendedConstructor()
{
  ScalarMatrix m;
  Dune::ImplicitMatrixBuilder<ScalarMatrix> b(m,10,10,3,0.1);
  setMatrix(b);
  m.compress();
  setMatrix(m);
}

int main()
{
  try{
    testImplicitBuild();
    testImplicitBuildWithInsufficientOverflow();
    testSetterInterface();
    testDoubleSetSize();
    testInvalidBuildModeConstructorCall();
    testNegativeOverflowConstructorCall();
    testInvalidSetImplicitBuildModeParameters();
    testSetImplicitBuildModeParametersAfterSetSize();
    testSetSizeWithNonzeroes();
    testCopyConstructionAndAssignment();
    testInvalidCopyConstruction();
    testEntryConsistency();
    testEntryAfterCompress();
    testBracketOperatorBeforeCompress();
    testConstBracketOperatorBeforeCompress();
    testImplicitMatrixBuilder();
    testImplicitMatrixBuilderExtendedConstructor();
  }catch(Dune::Exception& e) {
    std::cerr << e <<std::endl;
    return 1;
  }
}