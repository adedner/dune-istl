# .. cmake_module::
#
#    This modules content is executed whenever a module required or suggests dune-istl!
#

dune_find_package(METIS TARGET METIS::METIS)
dune_find_package(ParMETIS TARGET ParMETIS::ParMETIS)
include(AddParMETISFlags)

dune_find_package(SuperLU 5.0 TARGET SuperLU::SuperLU)
include(AddSuperLUFlags)

find_package(ARPACKPP)
include(AddARPACKPPFlags)

dune_find_package(SuiteSparse OPTIONAL_COMPONENTS CHOLMOD LDL SPQR UMFPACK
  TARGET SuiteSparse::SuiteSparse)
include(AddSuiteSparseFlags)

# enable / disable backwards compatibility w.r.t. category
set(DUNE_ISTL_SUPPORT_OLD_CATEGORY_INTERFACE 1
  CACHE BOOL "Enable/Disable the backwards compatibility of the category enum/method in dune-istl solvers, preconditioner, etc. '1'")
