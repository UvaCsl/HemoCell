cmake_minimum_required(VERSION 3.16...3.31)
project(packCells CXX)

# packCells is built on its own, not with HemoCell. The binary is written next
# to the sources (tools/packCells/packCells), where the documentation runs it,
# unless -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=<dir> says otherwise: give each of
# several build trees its own, or they overwrite each other's binary (and a
# build tree may then take another's binary as up to date).
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
option(PACKCELLS_NATIVE_ARCH "Optimise for the CPU of the building machine (-march=native)" ON)
option(PACKCELLS_WERROR "Treat compiler warnings as errors" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(packCells packCells.cpp)
# No fused multiply-adds: the compiler fuses them differently with and without
# OpenMP, with -march=native or not, and after small changes to the code, and
# the packing is chaotic, so every build gave different packings for a seed.
# Without them a seed gives the same packing with any of these builds, at no
# cost in speed (measured).
target_compile_options(packCells PRIVATE -Wall -Wextra -ffp-contract=off)
if(PACKCELLS_WERROR)
  target_compile_options(packCells PRIVATE -Werror)
endif()
if(PACKCELLS_NATIVE_ARCH)
  target_compile_options(packCells PRIVATE -march=native)
endif()

# OpenMP parallelises the contact search. Without it (e.g. Apple's clang)
# packCells runs on one thread, with the same results.
find_package(OpenMP COMPONENTS CXX)
if(OpenMP_CXX_FOUND)
  target_link_libraries(packCells PRIVATE OpenMP::OpenMP_CXX)
else()
  message(STATUS "packCells: OpenMP not found, building for one thread")
  target_compile_options(packCells PRIVATE -Wno-unknown-pragmas)
endif()
