cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# hemocell project definition
set(PROJECT_NAME hemocell)
project(${PROJECT_NAME} VERSION 2.6 LANGUAGES CXX C)

# keep track of base directories for source of hemocell and palabos
set(HEMOCELL_BASE_DIR ".")
set(PALABOS_DIR "${HEMOCELL_BASE_DIR}/palabos")

message(STATUS "HEMOCELL_DIR: " ${HEMOCELL_BASE_DIR} )
message(STATUS "PALABOS_DIR : " ${PALABOS_DIR} )

set(DEFAULT_BUILD_TYPE "Release")
set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose build type" FORCE)

# default values of dependency discovery
set(MPI_FOUND NO)
set(HDF5_FOUND NO)
set(PARMETIS_FOUND NO)

# Include CMake functions, mostly related to linking with dependencies
include(cmake/functions.cmake)
include(cmake/setup_mpi.cmake)
include(cmake/setup_hdf5.cmake)
include(cmake/setup_googletest.cmake)
include(cmake/find_parmetis.cmake)
include(cmake/setup_parmetis.cmake)

# define standard installation directories
include(GNUInstallDirs)

include(CTest)
if(BUILD_TESTING)
    # install testing framework if tests are required
    ConfigureGTest()
endif()

# reduced verbosity during compilation
set(CMAKE_VERBOSE_MAKEFILE 0)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

# FIXME: update these flags if needed
# pasumo options
ADD_DEFINITIONS("-DPLB_MPI_PARALLEL")
ADD_DEFINITIONS("-DPLB_USE_POSIX")
ADD_DEFINITIONS("-DPLB_SMP_PARALLEL")
ADD_DEFINITIONS("-DCOMPILING_HEMOCELL_LIBRARY")
IF(APPLE)
    # FIXME: the compilation on Mac OS seems broken at this moment
    ADD_DEFINITIONS("-DPLB_MAC_OS_X")
ENDIF(APPLE)

# all source directories of interest for the `hemocell` libraries
# FIXME: `palabos` should probably be compiled as a library and just linked
set(SRC_DIR_LIST
        ${CMAKE_CURRENT_SOURCE_DIR}/config
        ${CMAKE_CURRENT_SOURCE_DIR}/core
        ${CMAKE_CURRENT_SOURCE_DIR}/helper
        ${CMAKE_CURRENT_SOURCE_DIR}/io
        ${CMAKE_CURRENT_SOURCE_DIR}/mechanics
        ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyxml2
        ${CMAKE_CURRENT_SOURCE_DIR}/${PALABOS_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries/tinyxml
)

# Determine if optional dependencies are present on the system. This sets the
# `{PARMETIS_FOUND}` flag when present, which will be used to toggle optional
# behaviour during configuration and building.
FindParmetis()

# Define multiple variants of the hemocell library that correspond to differetn
# compilation options. These options are toggled by enabling different #define
# macros. The following variants are compiled:
#
# - libhemocell.a: the default hemocell library without any flags set
# - libhemocell_interior_viscosity.a: includes (optional) interior viscosity
# - libhemocell_solidify_mechanics.a: includes (optional) solidification
# - libhemocell_parmetis.a: includes (optional) link to PARMETIS
set(LIBRARY_TARGETS
        "${PROJECT_NAME}"
        "${PROJECT_NAME}_interior_viscosity"
        "${PROJECT_NAME}_solidify_mechanics"
)
if(${PARMETIS_FOUND})
        list(APPEND LIBRARY_TARGETS "${PROJECT_NAME}_parmetis")
endif()

# add all source directories to each library variant
foreach(TARGET ${LIBRARY_TARGETS})
        add_source_to_library(${TARGET} "${SRC_DIR_LIST}")
endforeach()

# assert compile definitions per library targets
target_compile_definitions("${PROJECT_NAME}_interior_viscosity" PUBLIC INTERIOR_VISCOSITY)
target_compile_definitions("${PROJECT_NAME}_solidify_mechanics" PUBLIC SOLIDIFY_MECHANICS)
if(${PARMETIS_FOUND})
        target_compile_definitions("${PROJECT_NAME}_parmetis" PUBLIC HEMO_PARMETIS)
endif()

# include additional directories from palabos/`externalLibraries`
list(APPEND SRC_DIR_LIST
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries
)

foreach(TARGET ${LIBRARY_TARGETS})
        target_include_directories(${TARGET} PUBLIC "${SRC_DIR_LIST}")
endforeach()

# link dependencies
ConfigureMPI("${LIBRARY_TARGETS}")  # updates `{CMAKE -> MPI}_CXX_COMPILER`
ConfigureHDF5("${LIBRARY_TARGETS}")
ConfigureParmetis("${PROJECT_NAME}_parmetis")

if(NOT (${MPI_FOUND} AND ${HDF5_FOUND}))
        message(FATAL_ERROR "\nOne or more required package (MPI, HDF5) not found.")
endif()

message(STATUS "Compiler: \
${CMAKE_CXX_COMPILER}, \
${CMAKE_CXX_COMPILER_ID}, \
${CMAKE_CXX_COMPILER_VERSION}.")

# default flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -Wformat -Wformat-security")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=format-security")

# Specialised flags
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-empty-body")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-qualifiers")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd858")
endif()

# compile examples
add_subdirectory(examples EXCLUDE_FROM_ALL)

# Add `hemocell/cases/` to compiled targets only if the `cases` dir exists to
# avoid CMake errors when the directory does not exist, which applies to the
# hemocell versions released through Github.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cases")
        add_subdirectory(cases EXCLUDE_FROM_ALL)
endif()

# compile tests
if(BUILD_TESTING)
        # to disable compilation provide `-DBUILD_TESTING=OFF` to CMake
        add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()
