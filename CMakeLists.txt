cmake_minimum_required(VERSION 3.13)
project("dune-istl" C CXX)

# guess build tree of dune-common
if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})

endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)

# start a dune project with information from dune.module
dune_project(dune-istl INTERFACE)

add_subdirectory(cmake/modules)
add_subdirectory(dune)
add_subdirectory(doc)

# if Python bindings are enabled, include necessary sub directories.
if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory(python)
  dune_python_install_package(PATH python)
endif()

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
