cmake_minimum_required(VERSION 2.8)
project(packEllipsoids)

# On MacOSX use Hombrew's g++ instead of the systems bundled clang
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_COMPILER "g++-7")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

# Compile with C++11 and OpenMP support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -Wall -Wextra -march=native -std=c++11 -fopenmp")

# Compile as a standalone application, not a library
add_definitions(-DSTANDALONE)

set(SOURCE_FILES packCells.cpp)
add_executable(packCells ${SOURCE_FILES})