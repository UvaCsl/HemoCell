cmake_minimum_required(VERSION 2.8)
project(packCells)

# On MacOsX use hombrew's g++ instead of clang
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_COMPILER "g++-5")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES ellipsoids.cpp)
add_executable(packCells ${SOURCE_FILES})