# Defines the functions to use SuperLU
#
# .. cmake_function:: add_dune_superlu_flags
#
#    .. cmake_param:: targets
#       :positional:
#       :single:
#       :required:
#
#       A list of targets to use SuperLU with.
#

# set HAVE_SUPERLU for config.h
set(HAVE_SUPERLU ${SUPERLU_FOUND})

# register all SuperLU related flags
if(SUPERLU_FOUND)
  dune_register_package_flags(
    COMPILE_DEFINITIONS "ENABLE_SUPERLU=1"
    LIBRARIES SuperLU::SuperLU)
endif()

# Provide function to set target properties for linking to SuperLU
function(add_dune_superlu_flags _targets)
  if(SUPERLU_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} SuperLU::SuperLU)
      target_compile_definitions(${_target} PUBLIC ENABLE_SUPERLU=1)
    endforeach()
  endif()
endfunction(add_dune_superlu_flags)

