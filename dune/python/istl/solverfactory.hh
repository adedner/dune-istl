#ifndef DUNE_PYTHON_ISTL_SOLVERFACTORY_HH
#define DUNE_PYTHON_ISTL_SOLVERFACTORY_HH

#include <dune/python/common/parametertree.hh>
#include <dune/istl/solverfactory.hh>
#include <dune/python/istl/solvers.hh>

namespace Dune {

  namespace Python {

    namespace {
      template<class SF>
      struct SFTraits;

      template<class Op>
      struct SFTraits<SolverFactory<Op>>{
        using Operator = Op;
        using InvOperator = InverseOperator<typename Op::domain_type, typename Op::range_type>;
      };

    }

    template<class SF, class... options>
    inline void registerSolverFactory( pybind11::module module,
                                       pybind11::class_< SF, options...> cls)
    {
      using Operator = typename SFTraits<SF>::Operator;
      initSolverFactories<Operator>();

      try{
        using InverseOp = typename SFTraits<SF>::InvOperator;
        pybind11::class_<InverseOp, std::shared_ptr<InverseOp>> invOpCls(module, "InverseOperator");
        registerInverseOperator(invOpCls);
      }catch(std::runtime_error& e){
        // InverseOperator is (probably) already registered
      }
      cls.def_static("get", [](std::shared_ptr<Operator> op, const ParameterTree& config){
        return SF::get(op, config);
      });
    }
  }
}

#endif
