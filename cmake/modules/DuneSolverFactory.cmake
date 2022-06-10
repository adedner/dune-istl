# Modules that provides functions for precompiling the solverfactory
#
# .. cmake_function:: dune_precompile_solverfactory
#
#    .. cmake_brief::
#
#       Adds a library with a global variable. It is initialized with a call to
#       initSolversFromFactory<OP> for every operator that is passed to this function
#
#    .. cmake_param:: REGISTER_GLOBALLY
#       :bool:
#
#       Whether the library should be registered to the dune flags via
#       `dune_register_package_flags`.
#
#    .. cmake_param:: NAME
#       :single:
#
#       Name of the created library target. If omitted the name is generated from
#       a hash of the OPERATORS and INCLUDES arguments.
#
#    .. cmake_param:: OPERATOR
#       :multi:
#
#       List of C++ type of ISTL Operators (e.g. MatrixAdapter<...>) for which
#       the solverfactory should be compiled.
#
#    .. cmake_param:: INCLUDES
#       :multi:
#
#       Headers that are included to compile the factory
#

include_guard(GLOBAL)

function(dune_precompile_solverfactory)
  set(BOOL_OPTS REGISTER_GLOBALLY)
  set(SINGLE_OPTS NAME)
  set(MULTI_OPTS OPERATORS INCLUDES)
  cmake_parse_arguments(arg "${BOOL_OPTS}" "${SINGLE_OPTS}" "${MULTI_OPTS}" ${ARGN})

  if(NOT arg_NAME)
    string(MD5 HASH "${arg_INCLUDES} ${arg_OPERATORS}")
    set(arg_NAME "solverfactory_${HASH}")
    if(NOT arg_REGISTER_GLOBALLY)
      message(WARNING "You dont specified a name for a precompile solverfactory
that is not registered globally. If you want to link explicitly against this library
you should specify a name. Otherwise use REGISTER_GLOBALLY.")
    endif()
  endif()

  ## Generate the lib content
  set(LIB_FILE "#include <config.h>\n")
  foreach(inc IN LISTS arg_INCLUDES)
    string(APPEND LIB_FILE "#include <${inc}>\n")
  endforeach()
  string(APPEND LIB_FILE "#include <dune/istl/solverfactory.hh>\n
int ${arg_NAME} = (")
  foreach(op IN LISTS arg_OPERATORS)
    string(APPEND LIB_FILE
      "Dune::initSolverFactories<${op}>(),")
  endforeach()
  string(APPEND LIB_FILE "0)\;\n")
  set(LIB_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${arg_NAME}.cc")
  file(WRITE "${LIB_FILENAME}" ${LIB_FILE})

  ## compile library
  dune_add_library(${arg_NAME} SOURCES "${LIB_FILENAME}" NO_MODULE_LIBRARY)
  add_dune_all_flags(${arg_NAME})
  get_target_property(target_type ${arg_NAME} TYPE)
  if (target_type STREQUAL STATIC_LIBRARY)
    target_link_options(${arg_NAME} PUBLIC "-u${arg_NAME}") # ensure that the library is linked
  else()
    target_link_options(${arg_NAME} PUBLIC "-Wl,--no-as-needed") # ensure that the library is linked
  endif()

  if(arg_REGISTER_GLOBALLY)
    dune_register_package_flags(LIBRARIES "${arg_NAME}")
  else()
    set_target_properties(${arg_NAME} PROPERTIES EXCLUDE_FROM_ALL 1)
  endif()
endfunction(dune_precompile_solverfactory)
