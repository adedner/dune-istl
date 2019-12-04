// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <dune/istl/solverrepository.hh>

// include all solvers that should be included:
#include <dune/istl/umfpack.hh>
#include <dune/istl/superlu.hh>

#include <dune/istl/bcrsmatrix.hh>

// Create Instances of the factories to add all the registered solvers to the factory
using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;

struct UniqueTag {};

struct Initializer {
  Initializer(){
    addRegistryToFactory<Dune::TypeList<Matrix, Vector, Vector>>
      (Dune::DirectSolverFactory<Matrix,Vector,Vector>::instance(), DirectSolverTag{});
  }
};
inline Initializer init;