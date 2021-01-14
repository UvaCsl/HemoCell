#=======================================
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)

#=======================================
# Build settings for hemocell 
SET(PROJECT_NAME hemocell)

PROJECT(${PROJECT_NAME} VERSION 2.2 LANGUAGES CXX C)

# HemoCell and palabos directories relative to CMakelists.txt
SET(HEMOCELL_BASE_DIR ".")
SET(PALABOS_DIR "${HEMOCELL_BASE_DIR}/palabos")

# TODO: set default build type

# Build settings
SET(ENABLE_MPI 1) # For windows builds set to zero
SET(ENABLE_PARMETIS 0) # Do we enable load_balancing stuff? 
SET(BUILD_TYPE Release)  # Debug or Release

MESSAGE( STATUS "HEMOCELL_DIR:         " ${HEMOCELL_BASE_DIR} )
MESSAGE( STATUS "PALABOS_DIR:          " ${PALABOS_DIR} )

#=======================================
# Set up cmake environment

include(cmake/functions.cmake)

SET(CMAKE_VERBOSE_MAKEFILE 1)
INCLUDE(GNUInstallDirs)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

#===================================
#Pasumo options 
ADD_DEFINITIONS("-DPLB_MPI_PARALLEL")
ADD_DEFINITIONS("-DPLB_USE_POSIX")
ADD_DEFINITIONS("-DPLB_SMP_PARALLEL")
ADD_DEFINITIONS("-DCOMPILING_HEMOCELL_LIBRARY")
IF(APPLE)
    ADD_DEFINITIONS("-DPLB_MAC_OS_X")
ENDIF(APPLE)

# FIXME: clean up dependencies below 

#===================================
# Link HDF5
SET(HDF5_PREFER_PARALLEL 1)
SET(HDF5_USE_STATIC_LIBRARIES 1)
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS C HL)
if(NOT ${HDF5_FOUND})
  message(fatal_error "Hdf5 Libraries not found!")
endif(NOT ${HDF5_FOUND})
INCLUDE_DIRECTORIES(AFTER ${HDF5_INCLUDE_DIRS})
message(STATUS "HDF5 Include dir:     ${HDF5_INCLUDE_DIRS}")

#===================================
# Link parmetis and metis; typically at `hemocell/external/`

IF(ENABLE_PARMETIS)
FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h HINTS
    ${HEMOCELL_BASE_DIR}/external/parmetis-4.0.3/include/)
IF(NOT PARMETIS_INCLUDE_DIR)
  MESSAGE(FATAL_ERROR "Parmetis include path not found!")
ENDIF(NOT PARMETIS_INCLUDE_DIR)
MESSAGE( STATUS "PARMETIS INCLUDE DIR:       " ${PARMETIS_INCLUDE_DIR})
FIND_PATH(METIS_INCLUDE_DIR metis.h HINTS
    ${HEMOCELL_BASE_DIR}/external/parmetis-4.0.3/metis/include/)
IF(NOT METIS_INCLUDE_DIR)
  MESSAGE(FATAL_ERROR "Metis include path not found!")
ENDIF(NOT METIS_INCLUDE_DIR)
MESSAGE( STATUS "METIS INCLUDE DIR:       " ${METIS_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})

FIND_LIBRARY(PARMETIS_LIBRARY libparmetis.a HINTS
    ${HEMOCELL_BASE_DIR}/external/parmetis-4.0.3/build/*/*) 
IF(PARMETIS_LIBRARY STREQUAL "PARMETIS_LIBRARY-NOTFOUND")
  MESSAGE(FATAL_ERROR "Parmetis Library not found, try compiling
  external/parmetis-4.0.3?!")
ENDIF(PARMETIS_LIBRARY STREQUAL "PARMETIS_LIBRARY-NOTFOUND")
MESSAGE( STATUS "PARMETIS LIBRARY DIR:       " ${PARMETIS_LIBRARY})

FIND_LIBRARY(METIS_LIBRARY libmetis.a HINTS
    ${HEMOCELL_BASE_DIR}/external/parmetis-4.0.3/build/*/*) 
IF(METIS_LIBRARY STREQUAL "METIS_LIBRARY-NOTFOUND")
  MESSAGE(FATAL_ERROR "Metis Library not found,try compiling
  external/parmetis-4.0.3? ")
ENDIF(METIS_LIBRARY STREQUAL "METIS_LIBRARY-NOTFOUND")
MESSAGE( STATUS "METIS LIBRARY DIR:       " ${METIS_LIBRARY})

ADD_DEFINITIONS("-DHEMO_PARMETIS")
ENDIF(ENABLE_PARMETIS)

#=======================================
# FIXME: move into separate function call? 
# Create tests
include(CTest)
if(BUILD_TESTING)
        ##################################
        # Download and install GoogleTest
        # see https://github.com/google/googletest/tree/v1.10.x/googletest
        configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
        execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
        if(result)
          message(FATAL_ERROR "CMake step for googletest failed: ${result}")
        endif()
        execute_process(COMMAND ${CMAKE_COMMAND} --build .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
        if(result)
          message(FATAL_ERROR "Build step for googletest failed: ${result}")
        endif()

        # Prevent overriding the parent project's compiler/linker
        # settings on Windows
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

        # Add googletest directly to our build. This defines
        # the gtest and gtest_main targets.
        add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                         ${CMAKE_CURRENT_BINARY_DIR}/googletest-build)

        MESSAGE(STATUS "GoogleTest: done")
endif()

#=======================================
# Detect MPI compiler 
OPTION(ENABLE_MPI "Enable MPI" ${DEFAULT})
INCLUDE(FindMPI)
IF(MPI_CXX_FOUND)
  SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
ELSE(MPI_CXX_FOUND)
  MESSAGE(FATAL_ERROR "MPI-COMPILER NOT found!")
ENDIF(MPI_CXX_FOUND)

#=======================================
# Report compiler and its ID
MESSAGE(STATUS "Compiler: ${MPI_CXX_COMPILER}, ${CMAKE_CXX_COMPILER_ID}, ${CMAKE_CXX_COMPILER_VERSION}")

# Set up correct flags for compiler
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -Wformat -Wformat-security")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=format-security")
IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-empty-body")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-qualifiers")
ELSEIF(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd858")
ENDIF()

#===================================
# all source directories of interest for `hemocell`
# FIXME: `palabos` should probably be compiled as a library and just linked
set(SRC_DIR_LIST 
        ${CMAKE_CURRENT_SOURCE_DIR}/config
        ${CMAKE_CURRENT_SOURCE_DIR}/core
        ${CMAKE_CURRENT_SOURCE_DIR}/helper
        ${CMAKE_CURRENT_SOURCE_DIR}/IO
        ${CMAKE_CURRENT_SOURCE_DIR}/mechanics
        ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyxml2
        ${CMAKE_CURRENT_SOURCE_DIR}/${PALABOS_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries/tinyxml
)

#===================================
# Define `hemocell` as static library
add_source_to_library(${PROJECT_NAME} "${SRC_DIR_LIST}")

#===================================
# Include directories 
LIST(APPEND SRC_DIR_LIST 
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries
)
target_include_directories(${PROJECT_NAME} PUBLIC ${SRC_DIR_LIST})

# link HDF5 dependency 
target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
INCLUDE_DIRECTORIES(AFTER ${HDF5_INCLUDE_DIRS})

# FIXME: create function to include the following dependencies

#IF(ENABLE_PARMETIS)
#  add_custom_target(hemocell_with_parmetis ALL 
#    COMMAND ${CMAKE_AR} crsT lib${PROJECT_NAME}.a
#    $<TARGET_FILE:${PROJECT_NAME}_pre_all_deps> ${PARMETIS_LIBRARY}
#    ${METIS_LIBRARY}
#    DEPENDS hemocell_pre_finalize)
#ELSE(ENABLE_PARMETIS)
#  add_custom_target(hemocell_with_parmetis ALL)
#ENDIF(ENABLE_PARMETIS)
    
#might work someday, not now, must be linked in the case
#add_custom_target(hemocell_with_hdf5 ALL 
#  COMMAND ${CMAKE_AR} crsT lib${PROJECT_NAME}.a
#  $<TARGET_FILE:${PROJECT_NAME}_pre_all_deps> ${HDF5_LIBRARIES})

# finalize library
#add_custom_target(hemocell_finalize ALL 
#  COMMAND ${CMAKE_AR} crsT lib${PROJECT_NAME}.a
#  $<TARGET_FILE:${PROJECT_NAME}_pre_all_deps>
#  DEPENDS hemocell_pre_finalize hemocell_with_parmetis)

#=======================================
# Create examples: 
# This ensures all examples are compiled
add_subdirectory(examples/pipeflow)

#=======================================
# Create tests:
# Builds all tests when testing is enabled, pass `-DBUILD_TESTING=OFF` to 
# disable test generation
if(BUILD_TESTING)
        MESSAGE( STATUS "making tests")
        add_subdirectory(tests)
endif()
