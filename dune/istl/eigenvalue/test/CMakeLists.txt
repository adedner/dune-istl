add_executable(poweriterationtest cond2test.cc)
target_link_libraries(poweriterationtest PRIVATE dune-istl)
dune_add_test(
  TARGET poweriterationtest
  CMD_ARGS 40)

if(SuperLU_FOUND)
  add_executable(poweriterationsuperlutest cond2test.cc)
  add_dune_superlu_flags(poweriterationsuperlutest)
  dune_add_test(
    TARGET poweriterationsuperlutest
    CMD_ARGS 40)
  target_link_libraries(poweriterationsuperlutest PRIVATE dune-istl)
endif()

if(ARPACKPP_FOUND)
  add_executable(arpackpptest cond2test.cc)
  add_dune_arpackpp_flags(arpackpptest)
  target_link_libraries(arpackpptest PRIVATE dune-istl)
  dune_add_test(TARGET arpackpptest)
endif()

if(ARPACKPP_FOUND AND SuperLU_FOUND)
  dune_add_test(NAME arpackppsuperlutest
                SOURCES cond2test.cc)
endif()
