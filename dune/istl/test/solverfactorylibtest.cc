#include <config.h>

#include <iostream>

#include <dune/common/parametertreeparser.hh>

#include <dune/istl/solverfactory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/bcrsmatrix.hh>

using namespace Dune;

int main(int argc, char** argv){
  BCRSMatrix<double> mat;
  MatrixAdapter<BCRSMatrix<double>, BlockVector<double>, BlockVector<double>> op(mat);
  ParameterTree config;
  config["type"] = "cgsolver";
  config["verbose"] = "0";
  config["maxit"] = "10";
  config["reduction"] = "1e-5";
  config.sub("preconditioner")["type"] = "ssor";
  ParameterTreeParser::readOptions(argc, argv, config);
  auto solver = getSolverFromFactory(stackobject_to_shared_ptr(op), config);
  std::cout << Dune::className(*solver) << std::endl;
  return 0;
}
