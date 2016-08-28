#=======================================
# Build settings. Set these.

# Project setup
SET(PROJECT_NAME stretchCell)
SET(PROJECT_SRC "stretchCell.cpp")

# Ficsion location relative to CMakelists.txt
SET(FICSION_BASE_DIR "./../..")

# Palabos setup
SET(PALABOS_DIR "${FICSION_BASE_DIR}/palabos")
SET(PALABOS_MAJOR_VERSION 1)
SET(PALABOS_MINOR_VERSION 5)
SET(PALABOS_PATCH_VERSION 1)

# Build settings
SET(ENABLE_MPI 1)
SET(ENABLE_POSIX 1)
SET(SMP_PARALLEL 1)
SET(BUILD_TYPE Release)  # Debug or Release

# Required additional external libraries
SET(NEEDHDF5 1)  # For highlevel HDF5 routines
SET(NEEDVTK 0)   # For VTK
SET(NEEDTXML 0)  # Set to 1 to enable external TinyXml libs (Palabos 1.5.0+)
SET(NEEDEIGEN 0) # Set to 1 to use external Eigen3 libraries

MESSAGE( STATUS "PALABOS_DIR:         " ${PALABOS_DIR} )
MESSAGE( STATUS "FICSION_BASE_DIR:         " ${FICSION_BASE_DIR} )

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

IF(NEEDTXML)
  FIND_PACKAGE(TinyXML REQUIRED HINTS ${FICSION_BASE_DIR}/cmake )
  INCLUDE_DIRECTORIES(AFTER ${TINYXML_INCLUDE_DIR})
ENDIF(NEEDTXML)

IF(NEEDEIGEN)
  FIND_PACKAGE(Eigen3 REQUIRED)
  INCLUDE_DIRECTORIES(AFTER ${EIGEN3_INCLUDE_DIR})
ENDIF(NEEDEIGEN)

FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(AFTER ${MPI_INCLUDE_PATH})

IF(NEEDHDF5)
  # This is needed since stupid Lisa has a bug in 3.5.1
  IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
    FIND_PACKAGE(HDF5 COMPONENTS CXX HL REQUIRED)
  ELSE (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
    FIND_PACKAGE(HDF5 REQUIRED)
  ENDIF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
  
  INCLUDE_DIRECTORIES(AFTER ${HDF5_INCLUDE_DIRS})
ENDIF(NEEDHDF5)

IF(NEEDVTK)
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${VTK_USE_FILE})
ENDIF(NEEDVTK)

#=======================================

IF (NOT CMAKE_CXX_FLAGS)
  #If flags are not set, add default flags
  SET(CMAKE_BUILD_TYPE ${BUILD_TYPE})
  # GNU C++98 feature set must be forced on new compilers because of Palabos and its dependencies.
  SET(CMAKE_CXX_FLAGS  " -std=gnu++98 -Wall -g -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security")
ENDIF (NOT CMAKE_CXX_FLAGS)

#=======================================

SET(DEFAULT ON CACHE INTERNAL "Default value for enabled by default options")
OPTION(ENABLE_MPI "Enable MPI" ${DEFAULT})

IF(ENABLE_MPI)
  INCLUDE(FindMPI)
  IF(MPI_C_FOUND AND MPI_CXX_FOUND)
    SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    ADD_DEFINITIONS("-DPLB_MPI_PARALLEL")
  ELSE(MPI_C_FOUND AND MPI_CXX_FOUND)
    #Check old FindMPI version
    IF(MPI_COMPILER)
      SET(CMAKE_C_COMPILER ${MPI_COMPILER})
      SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
      ADD_DEFINITIONS("-DPLB_MPI_PARALLEL")
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
    MESSAGE("MacOS X detected: disable OpenMP or set the compiler to g++-6!")
    #SET(CMAKE_C_COMPILER gcc-6)
    #SET(CMAKE_CXX_COMPILER g++-6)
    MESSAGE("Current choice: OpenMP disabled.")
ELSE(APPLE)
  #find_package(OpenMP)
  #if (OPENMP_FOUND)
  #    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   #  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  #endif(OPENMP_FOUND)
ENDIF(APPLE)

INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/src)

INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${FICSION_BASE_DIR})
INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${FICSION_BASE_DIR}/IO)
INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${FICSION_BASE_DIR}/core)
INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${FICSION_BASE_DIR}/extensions)
INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${FICSION_BASE_DIR}/external)
INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${FICSION_BASE_DIR}/helper)
INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${FICSION_BASE_DIR}/models)

IF(NOT NEEDEIGEN)
  #INCLUDE(${PALABOS_DIR}/externalLibraries/Eigen/CMakeLists.txt)
  INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries/Eigen3)
ENDIF(NOT NEEDEIGEN)

FILE(GLOB_RECURSE SRC_FILES_SRC "${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/src/*.cpp")

IF(NOT NEEDTXML)
  INCLUDE_DIRECTORIES(AFTER ${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries)
  LIST(APPEND SRC_FILES_SRC ${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries/tinyxml/tinyxml.cpp)
  LIST(APPEND SRC_FILES_SRC ${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries/tinyxml/tinyxmlerror.cpp)
  LIST(APPEND SRC_FILES_SRC ${CMAKE_SOURCE_DIR}/${PALABOS_DIR}/externalLibraries/tinyxml/tinyxmlparser.cpp)
ENDIF(NOT NEEDTXML)

LIST(APPEND SRC_FILES_SRC ${PROJECT_SRC})

ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_FILES_SRC})

IF(NEEDHDF5)
    # This is needed since stupid Lisa has a bug in 3.5.1
  IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${HDF5_LIBRARIES})
  ELSE (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} hdf5_hl hdf5)
  ENDIF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)

  # Temporary fix for bugous library handling on osx brew
  IF(APPLE)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} /usr/local/lib/libhdf5_hl.a)
  ENDIF(APPLE)
ENDIF(NEEDHDF5)

IF(NEEDVTK)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${VTK_LIBRARIES})
ENDIF(NEEDVTK)

IF(NEEDTXML)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${TINYXML_LIBRARIES})
ENDIF(NEEDTXML)


