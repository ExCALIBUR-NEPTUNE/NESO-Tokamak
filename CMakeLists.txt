cmake_minimum_required(VERSION 3.21)

project(
  NESOTOKAMAK
  VERSION 0.0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# ##############################################################################
# Find dependencies
# ##############################################################################
find_package(SYCL REQUIRED)
if(SYCL_FOUND)
  message(STATUS "Using SYCL implementation ${SYCL_IMPLEMENTATION}")
  if(hipSYCL_FOUND)
    set(HIPSYCL_TARGETS "omp")
    set(SYCL_FLAG NESO_HIPSYCL)
  elseif(SYCL_IMPLEMENTATION STREQUAL "IntelDPCPP")
    set(SYCL_FLAG NESO_DPCPP)
  endif()
endif()
find_package(NESO REQUIRED)

# Use a custom FindNektar++ script to provide an interface target
find_package(Nektar++ REQUIRED)
# find MPI
find_package(MPI REQUIRED)
# Add the NESO-Particles dependencies
find_package(NESO-Particles REQUIRED)

# Identify source files
set(SRC_FILES
src/EquationSystems/TokamakSystem.cpp
src/BoundaryConditions/TokamakBndCond.cpp
src/BoundaryConditions/ObliqueBC.cpp
)

# ##############################################################################
# Build object library; used by solver exec and tests
# ##############################################################################
set(OBJLIB_NAME ObjLib)
set(SOLVER_LIBS
    ${SOLVER_LIBS} ${OBJLIB_NAME}
    CACHE INTERNAL "")
add_library(${OBJLIB_NAME} OBJECT ${SRC_FILES})
target_compile_options(${OBJLIB_NAME} PRIVATE ${BUILD_TYPE_COMPILE_FLAGS})
target_link_libraries(${OBJLIB_NAME} PRIVATE Nektar++::nektar++ NESO::nesolib)

target_include_directories(${OBJLIB_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR}
                                                 ${NESO_PARTICLES_INCLUDE_PATH})

add_sycl_to_target(TARGET ${OBJLIB_NAME} SOURCES ${SRC_FILES})

# ##############################################################################
# Build solver exec
# ##############################################################################
set(EXEC_NAME tokamak)

add_executable(${EXEC_NAME} src/main.cpp $<TARGET_OBJECTS:${OBJLIB_NAME}>)

# Compile options
target_compile_options(${EXEC_NAME} PRIVATE ${BUILD_TYPE_COMPILE_FLAGS})

target_include_directories(${EXEC_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR}
                                               ${NESO_PARTICLES_INCLUDE_PATH})

# Linker options, target libs
target_link_options(${EXEC_NAME} PRIVATE ${BUILD_TYPE_LINK_FLAGS})
target_link_libraries(
  ${EXEC_NAME} PRIVATE Nektar++::nektar++ ${NESO_PARTICLES_LIBRARIES}
                       NESO::nesolib)

add_sycl_to_target(TARGET ${EXEC_NAME} SOURCES src/main.cpp)

# ##############################################################################
# Install exec
# ##############################################################################
# Set default install path if not provided on the command-line
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
      ${CMAKE_BINARY_DIR}/dist
      CACHE PATH "" FORCE)
endif()

# Install location
install(TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# Exec rpath
set_target_properties(${EXEC_NAME} PROPERTIES INSTALL_RPATH ${INSTALL_RPATH}
                                              ${NEKTAR++_LIBRARY_DIRS})
