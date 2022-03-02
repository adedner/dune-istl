# Modules that provides functions for precompiling the solverfactory
#
# .. cmake_function:: dune_precompile_solverfactory
#
#    .. cmake_brief::
#
#       Adds a library with a global variable. It is initialized with a call to
#       initSolversFromFactory<OP> for every operator that is passed to this function
#
#    .. cmake_param:: OPERATOR
#       :multi:
#
#       C++ typenames of Operators for that the solverfactory should be compiled.
#
#    .. cmake_param:: INCLUDES
#       :multi:
#
#       Headers that are included to compile the factory
#

include_guard(GLOBAL)

function(dune_precompile_solverfactory)
  set(MULTI_OPTS OPERATORS INCLUDES)
  cmake_parse_arguments(arg "" "" "${MULTI_OPTS}" ${ARGN})

  string(MD5 HASH "${arg_INCLUDES} ${arg_OPERATORS}")
  set(NAME "solverfactory_${HASH}")

  ## Generate the lib content
  set(LIB_FILE "#include <config.h>\n")
  foreach(inc IN LISTS arg_INCLUDES)
    string(APPEND LIB_FILE "#include <${inc}>\n")
  endforeach()
  string(APPEND LIB_FILE "#include <dune/istl/solverfactory.hh>\n
int ${NAME} = (")
  foreach(op IN LISTS arg_OPERATORS)
    string(APPEND LIB_FILE
      "Dune::initSolverFactories<${op}>(),")
  endforeach()
  string(APPEND LIB_FILE "0)\;\n")
  set(LIB_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${NAME}.cc")
  file(WRITE "${LIB_FILENAME}" ${LIB_FILE})

  ## compile library
  dune_add_library(${NAME} SOURCES "${LIB_FILENAME}")
  add_dune_all_flags(${NAME})
  get_target_property(target_type ${NAME} TYPE)
  if (target_type STREQUAL STATIC_LIBRARY)
    target_link_options(${NAME} PUBLIC "-u${NAME}") # ensure that the library is linked
  else()
    target_link_options(${NAME} PUBLIC "-Wl,--no-as-needed") # ensure that the library is linked
  endif()

  dune_register_package_flags(LIBRARIES "${NAME}")
endfunction(dune_precompile_solverfactory)
