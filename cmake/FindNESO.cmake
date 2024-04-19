# FindNESO.cmake
#
# Wrapper for finding NESO installations. This imports the provided NESO CMake
# config file and then creates a target for it.
#
# This will define the following variables
#
# NESO_FOUND All other variables provided by the default NESO CMake
#
# and the targets
#
# NESO::NESO
#

find_package(NESO CONFIG)

if(NESO_FOUND AND NOT TARGET NESO::NESO)

  # set(NESO_TP_FILTERED_INCLUDE_DIRS ${NESO_TP_INCLUDE_DIRS}) list(REMOVE_ITEM
  # NESO_TP_FILTERED_INCLUDE_DIRS "include/NESO")

  add_library(NESO::NESO INTERFACE IMPORTED)
  target_link_libraries(NESO::NESO INTERFACE ${NESO_LIBRARIES})
  # target_compile_definitions( NESO::NESO INTERFACE ${NESO_DEFINITIONS}
  # ${NESO_GENERATED_DEFINITIONS}) target_include_directories( NESO::NESO
  # INTERFACE ${NESO_INCLUDE_DIRS} ${NESO_TP_FILTERED_INCLUDE_DIRS}
  # ${SOLVER_INC_DIR})
endif()
