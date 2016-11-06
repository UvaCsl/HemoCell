#=======================================
# Build settings. Set these.

# Project setup
SET(PROJECT_NAME hemocell)

# HemoCell location relative to CMakelists.txt
SET(HEMOCELL_BASE_DIR "./../..")

# Palabos setup
SET(PALABOS_DIR "${HEMOCELL_BASE_DIR}/build/palabos")
SET(PALABOS_BASE_DIR "${HEMOCELL_BASE_DIR}/palabos")

# Build settings
SET(ENABLE_MPI 1)
SET(ENABLE_POSIX 1)
SET(SMP_PARALLEL 1)
SET(BUILD_TYPE Release)  # Debug or Release

# Required additional external libraries
SET(NEEDHDF5 1)  # For highlevel HDF5 routines
SET(NEEDTXML 0)  # Set to 1 to enable external TinyXml libs (Palabos 1.5.0+)

MESSAGE( STATUS "PALABOS_DIR:         " ${PALABOS_DIR} )
MESSAGE( STATUS "HEMOCELL_BASE_DIR:         " ${HEMOCELL_BASE_DIR} )

### Workaround for building with external cmake (e.g., when using an IDE with built in cmake):
# Access system wide libraries in ubuntu 14.04 in case one uses an external cmake ditribution (e.g., CLion)
#SET(CMAKE_MODULE_PATH /usr/share/cmake-2.8/Modules)  # Needed to find installed libraries on the system with the bundled cmake

#=======================================
##### Beginning of build script
##### No need to edit beyond this point

PROJECT(${PROJECT_NAME} CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_VERBOSE_MAKEFILE 1)

INCLUDE(GNUInstallDirs)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

#=======================================
# Set up cmake environment

FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

IF(NEEDTXML)
  FIND_PACKAGE(TinyXML REQUIRED HINTS ${HEMOCELL_BASE_DIR}/cmake )
  INCLUDE_DIRECTORIES(${TINYXML_INCLUDE_DIR})
ENDIF(NEEDTXML)

IF(NOT NEEDTXML)
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/${PALABOS_BASE_DIR}/externalLibraries)
#LIST(APPEND SRC_FILES_SRC ${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries/tinyxml/tinyxml.cpp)
#  LIST(APPEND SRC_FILES_SRC ${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries/tinyxml/tinyxmlerror.cpp)
#  LIST(APPEND SRC_FILES_SRC ${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries/tinyxml/tinyxmlparser.cpp)
ENDIF(NOT NEEDTXML)

IF(NEEDHDF5)
  # This is needed since stupid Lisa has a bug in 3.5.1
  IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
    FIND_PACKAGE(HDF5 COMPONENTS C HL REQUIRED)
  ELSE (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
    FIND_PACKAGE(HDF5 REQUIRED)
  ENDIF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
  
  INCLUDE_DIRECTORIES(AFTER ${HDF5_INCLUDE_DIRS})
ENDIF(NEEDHDF5)

#=======================================

IF (NOT CMAKE_CXX_FLAGS)
  #If flags are not set, add default flags
  SET(CMAKE_BUILD_TYPE ${BUILD_TYPE})
  SET(CMAKE_CXX_FLAGS  " -Wall -Wextra -g -Wformat -Wformat-security -Werror=format-security -O3 -march=native")
ENDIF (NOT CMAKE_CXX_FLAGS)

#=======================================

SET(DEFAULT ON CACHE INTERNAL "Default value for enabled by default options")
OPTION(ENABLE_MPI "Enable MPI" ${DEFAULT})

IF(ENABLE_MPI)
  INCLUDE(FindMPI)
  IF(MPI_C_FOUND AND MPI_CXX_FOUND)
    SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  ELSE(MPI_C_FOUND AND MPI_CXX_FOUND)
    #Check old FindMPI version
    IF(MPI_COMPILER)
      SET(CMAKE_C_COMPILER ${MPI_COMPILER})
      SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
    ELSE(MPI_COMPILER)
      MESSAGE(FATAL_ERROR "MPI-COMPILER NOT found!")
    ENDIF(MPI_COMPILER)
  ENDIF(MPI_C_FOUND AND MPI_CXX_FOUND)
ENDIF(ENABLE_MPI)

#=======================================

OPTION(ENABLE_POSIX "Enable POSIX" ${DEFAULT})

IF(ENABLE_POSIX)
  ADD_DEFINITIONS("-DPLB_USE_POSIX")
ENDIF(ENABLE_POSIX)

#=======================================

OPTION(ENABLE_SMP_PARALLEL "Enable SMP_PARALLEL" ${DEFAULT})

IF(ENABLE_SMP_PARALLEL)
  ADD_DEFINITIONS("-DPLB_SMP_PARALLEL")
ENDIF(ENABLE_SMP_PARALLEL)

#=======================================

IF(APPLE)
    ADD_DEFINITIONS("-DPLB_MAC_OS_X")
ENDIF(APPLE)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/IO)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/core)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/extensions)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/external)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/helper)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/models)
FILE(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/IO/*.cpp"
                     "${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/core/*.cpp"
                     "${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/extensions/*.cpp"
                     "${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/external/*.cpp"
                     "${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/helper/*.cpp"
                     "${CMAKE_SOURCE_DIR}/${HEMOCELL_BASE_DIR}/models/*.cpp")
include(ExternalProject) 
ExternalProject_Add("palabos" PREFIX ${PALABOS_DIR} SOURCE_DIR ${PALABOS_DIR}
    BINARY_DIR ${PALABOS_DIR} INSTALL_COMMAND "")
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/${PALABOS_BASE_DIR}/src)


ADD_library(${PROJECT_NAME} STATIC ${SRC_FILES})
