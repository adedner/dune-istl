
add_subdirectory("test")

install(FILES
  blas.hh
  blockbicgstab.hh
  blockcg.hh
  blockgmres.hh
  blockinnerproduct.hh
  matrixalgebra.hh
  utils.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/istl/blockkrylov)
