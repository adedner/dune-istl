add_python_targets(istl
  __init__
)
dune_add_pybind11_module(NAME _istl)
target_link_libraries(_istl PUBLIC Dune::istl)
