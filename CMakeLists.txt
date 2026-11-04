cmake_minimum_required(VERSION 3.16...3.31 FATAL_ERROR)

# HemoCell. Configure out of source, for example:
#   cmake -S . -B build && cmake --build build --target hemocell
# or use a preset (CMake >= 3.21): cmake --preset release
project(hemocell VERSION 2.8 LANGUAGES C CXX)

set(PALABOS_DIR "${PROJECT_SOURCE_DIR}/palabos")
if(NOT EXISTS "${PALABOS_DIR}/src/palabos3D.h")
  message(FATAL_ERROR "Palabos not found in ${PALABOS_DIR}. Run ./setup.sh first.")
endif()

# Release by default (single-configuration generators). Unlike before, an
# explicitly requested build type (e.g. -DCMAKE_BUILD_TYPE=Debug) is honoured.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()
# Write compile_commands.json for tooling unless explicitly disabled (CMake
# creates this cache entry empty when it is not set).
if("${CMAKE_EXPORT_COMPILE_COMMANDS}" STREQUAL "")
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Write compile_commands.json" FORCE)
endif()

# Executables are written next to their sources (cases are run from there).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

# Library and example code is C++11; the tests use C++17 (GoogleTest 1.17).
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(cmake/HemoCellCompilerOptions.cmake)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS C CXX)

set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 REQUIRED COMPONENTS C HL)
# HDF5 as an interface target that is linked *after* MPI, so that MPI's include
# directories come first. An HDF5 installation prefix (e.g. a conda
# environment) can contain the headers of a different MPI, whose mpi.h must not
# shadow the one HemoCell and Palabos are compiled with.
add_library(hemocell_hdf5 INTERFACE)
if(HDF5_INCLUDE_DIRS)
  target_include_directories(hemocell_hdf5 SYSTEM INTERFACE ${HDF5_INCLUDE_DIRS})
endif()
target_link_libraries(hemocell_hdf5 INTERFACE ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})

include(cmake/find_parmetis.cmake)
FindParmetis()

include(CTest) # defines BUILD_TESTING (ON by default)
if(BUILD_TESTING)
  include(cmake/setup_googletest.cmake)
  ConfigureGTest()
endif()

# ---------------------------------------------------------------------------
# Palabos: compiled once (object library); its objects are archived into each
# HemoCell library variant, which therefore stays self-contained.
# ---------------------------------------------------------------------------
file(GLOB_RECURSE PALABOS_SOURCES CONFIGURE_DEPENDS
  "${PALABOS_DIR}/src/*.cpp"
  "${PALABOS_DIR}/externalLibraries/tinyxml/*.cpp")

# The HemoCell patch makes this file call HemoCell's particle data transfer,
# so it is compiled with the HemoCell sources (Palabos objects stay free of
# HemoCell symbols and independent of the variant macros).
set(PALABOS_HEMOCELL_SOURCES "${PALABOS_DIR}/src/parallelism/parallelBlockCommunicator3D.cpp")
list(REMOVE_ITEM PALABOS_SOURCES ${PALABOS_HEMOCELL_SOURCES})

add_library(palabos OBJECT ${PALABOS_SOURCES})
# SYSTEM: warnings from Palabos headers are not reported for HemoCell code.
target_include_directories(palabos SYSTEM PUBLIC
  "${PALABOS_DIR}/src"
  "${PALABOS_DIR}/externalLibraries"
  "${PALABOS_DIR}/externalLibraries/tinyxml")
target_compile_definitions(palabos PUBLIC
  PLB_MPI_PARALLEL
  PLB_USE_POSIX
  $<$<BOOL:${APPLE}>:PLB_MAC_OS_X>
  $<$<BOOL:${HEMOCELL_PALABOS_DEBUG}>:PLB_DEBUG>)
target_link_libraries(palabos PUBLIC MPI::MPI_CXX)
# Vendored code: do not report its warnings.
target_compile_options(palabos PRIVATE -w)

# ---------------------------------------------------------------------------
# HemoCell library variants
# ---------------------------------------------------------------------------
file(GLOB_RECURSE HEMOCELL_SOURCES CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/config/*.cpp"
  "${PROJECT_SOURCE_DIR}/core/*.cpp"
  "${PROJECT_SOURCE_DIR}/helper/*.cpp"
  "${PROJECT_SOURCE_DIR}/io/*.cpp"
  "${PROJECT_SOURCE_DIR}/mechanics/*.cpp"
  "${PROJECT_SOURCE_DIR}/external/tinyxml2/*.cpp")

# Vendored TinyXML-2: do not report its warnings.
set_source_files_properties("${PROJECT_SOURCE_DIR}/external/tinyxml2/tinyxml2.cpp"
  PROPERTIES COMPILE_OPTIONS "-w")

# hemocell_add_library(<name> [<compile definition>...])
# Feature definitions are PUBLIC so that clients see the same class layouts.
function(hemocell_add_library name)
  if(name STREQUAL "hemocell")
    add_library(${name} STATIC ${HEMOCELL_SOURCES} ${PALABOS_HEMOCELL_SOURCES})
  else()
    # Variants are only built when a case asks for them.
    add_library(${name} STATIC EXCLUDE_FROM_ALL ${HEMOCELL_SOURCES} ${PALABOS_HEMOCELL_SOURCES})
  endif()
  target_include_directories(${name} PUBLIC
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/config"
    "${PROJECT_SOURCE_DIR}/core"
    "${PROJECT_SOURCE_DIR}/helper"
    "${PROJECT_SOURCE_DIR}/io"
    "${PROJECT_SOURCE_DIR}/mechanics"
    "${PROJECT_SOURCE_DIR}/external/tinyxml2")
  # COMPILING_HEMOCELL_LIBRARY has always been defined for cases as well; it
  # disables the convenience includes/using-directives at the end of hemocell.h.
  target_compile_definitions(${name} PUBLIC COMPILING_HEMOCELL_LIBRARY ${ARGN})
  # Linking the object library adds the Palabos objects to this archive and
  # passes on its usage requirements (include directories, PLB_* macros).
  target_link_libraries(${name} PUBLIC palabos MPI::MPI_CXX hemocell_hdf5)
  hemocell_enable_warnings(${name})
endfunction()

hemocell_add_library(hemocell)
hemocell_add_library(hemocell_interior_viscosity INTERIOR_VISCOSITY)
hemocell_add_library(hemocell_solidify_mechanics SOLIDIFY_MECHANICS)
if(PARMETIS_FOUND)
  hemocell_add_library(hemocell_parmetis HEMO_PARMETIS)
  target_include_directories(hemocell_parmetis PRIVATE ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR})
  target_link_libraries(hemocell_parmetis PUBLIC ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
endif()

set(_hdf5_parallel no)
if(HDF5_IS_PARALLEL)
  set(_hdf5_parallel yes)
endif()
set(_parmetis no)
if(PARMETIS_FOUND)
  set(_parmetis yes)
endif()
message(STATUS "HemoCell ${PROJECT_VERSION}: ${CMAKE_BUILD_TYPE} build with "
               "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}, "
               "native arch: ${HEMOCELL_NATIVE_ARCH}, sanitizers: '${HEMOCELL_SANITIZERS}', "
               "HDF5 ${HDF5_VERSION} (parallel: ${_hdf5_parallel}), ParMETIS: ${_parmetis}")

# ---------------------------------------------------------------------------
# Cases and tests (not part of the default build; request targets explicitly,
# e.g. `--target pipeflow`, `--target hemocell_examples`, `--target hemocell_test`)
# ---------------------------------------------------------------------------
add_subdirectory(examples EXCLUDE_FROM_ALL)
# `cases/` is not part of every HemoCell distribution.
if(EXISTS "${PROJECT_SOURCE_DIR}/cases/CMakeLists.txt")
  add_subdirectory(cases EXCLUDE_FROM_ALL)
endif()
if(BUILD_TESTING)
  add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()
