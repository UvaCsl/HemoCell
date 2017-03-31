#=======================================
# Build settings. Set these.

# Project setup
SET(PROJECT_NAME pipeflow)
SET(PROJECT_SRC "pipeflow.cpp")

# HemoCell location relative to CMakelists.txt
SET(HEMOCELL_BASE_DIR "./../../")
get_filename_component(HEMOCELL_BASE_DIR "${HEMOCELL_BASE_DIR}" ABSOLUTE)
SET(HEMOCELL_DIR "${HEMOCELL_BASE_DIR}/build/hemocell")
SET(PALABOS_BASE_DIR "${HEMOCELL_BASE_DIR}/palabos")
SET(PALABOS_DIR "${HEMOCELL_BASE_DIR}/build/palabos")


MESSAGE( STATUS "HEMOCELL_BASE_DIR:         " ${HEMOCELL_BASE_DIR} )
MESSAGE( STATUS "HEMOCELL_DIR:         " ${HEMOCELL_DIR} )
MESSAGE( STATUS "PALABOS_BASE_DIR:         " ${PALABOS_BASE_DIR} )
MESSAGE( STATUS "PALABOS_DIR:         " ${PALABOS_DIR} )

### Workaround for building with external cmake (e.g., when using an IDE with built in cmake):
# Access system wide libraries in ubuntu 14.04 in case one uses an external cmake ditribution (e.g., CLion)
#SET(CMAKE_MODULE_PATH /usr/share/cmake-2.8/Modules)  # Needed to find installed libraries on the system with the bundled cmake

SET(ENABLE_MPI 1)  
SET(NEEDTXML 0)  # Set to 1 to enable external TinyXml libs (Palabos 1.5.0+)
SET(NEEDHDF5 1)  # For highlevel HDF5 routines
#SET(HDF5_USE_STATIC_LIBRARIES ON) #does not work anyway :(

#=======================================
##### Beginning of build script
##### No need to edit beyond this point

PROJECT(${PROJECT_NAME} CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_VERBOSE_MAKEFILE 1)

INCLUDE(GNUInstallDirs)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

IF (NOT CMAKE_CXX_FLAGS)
  #If flags are not set, add default flags
  SET(CMAKE_BUILD_TYPE ${BUILD_TYPE})
  SET(CMAKE_CXX_FLAGS  "-Wno-empty-body -Wno-unused-variable -Wno-deprecated-declarations -Wno-unknown-pragmas -Wno-unused-parameter -Wno-unused-result -Wall -Wextra -ggdb -Wformat -Wformat-security -Werror=format-security -O3 -march=native -std=c++11")
ENDIF (NOT CMAKE_CXX_FLAGS)
#SET (CMAKE_CXX_STANDARD 11)

#=======================================
IF(NEEDTXML)
  FIND_PACKAGE(TinyXML REQUIRED HINTS ${HEMOCELL_BASE_DIR}/cmake )
  INCLUDE_DIRECTORIES(${TINYXML_INCLUDE_DIR})
ENDIF(NEEDTXML)

IF(NOT NEEDTXML)
  INCLUDE_DIRECTORIES(${PALABOS_BASE_DIR}/externalLibraries)
  INCLUDE_DIRECTORIES(${PALABOS_BASE_DIR}/src/libraryInterfaces)
ENDIF(NOT NEEDTXML)
IF(NEEDHDF5)
    FIND_PACKAGE(HDF5 COMPONENTS C HL)
  INCLUDE_DIRECTORIES(AFTER ${HDF5_INCLUDE_DIRS})
ENDIF(NEEDHDF5)
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


SET(DEFAULT ON CACHE INTERNAL "Default value for enabled by default options")
OPTION(ENABLE_MPI "Enable MPI" ${DEFAULT})

INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/helper)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/config)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/core)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/models)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/mechanics)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/tools/dynpacking)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/external)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/IO)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/extensions)
INCLUDE_DIRECTORIES(${PALABOS_BASE_DIR}/src)
FILE(GLOB SRC_FILES "${HEMOCELL_BASE_DIR}/IO/*.cpp"
                    "${HEMOCELL_BASE_DIR}/core/*.cpp"
                    "${HEMOCELL_BASE_DIR}/config/*.cpp"
                    "${HEMOCELL_BASE_DIR}/extensions/*.cpp"
                    "${HEMOCELL_BASE_DIR}/external/*.cpp"
                    "${HEMOCELL_BASE_DIR}/helper/*.cpp"
                    "${HEMOCELL_BASE_DIR}/models/*.cpp"
                    "${HEMOCELL_BASE_DIR}/mechanics/*.cpp")

LIST(APPEND SRC_FILES ${PROJECT_SRC})

include(ExternalProject) 
ExternalProject_Add("palabos" PREFIX ${PALABOS_DIR} SOURCE_DIR ${PALABOS_DIR}
    BINARY_DIR ${PALABOS_DIR} INSTALL_COMMAND "")


#SET(CMAKE_EXE_LINKER_FLAGS "-static") !Good luck on finding statically linked mpi libraries, it should be possible though
ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_FILES})
ADD_DEPENDENCIES(${PROJECT_NAME} "palabos")
target_link_libraries(${PROJECT_NAME} ${PALABOS_DIR}/libpalabos.a)
IF(NEEDHDF5)
  # This is needed since stupid Lisa has a bug in 3.5.1
  IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${HDF5_LIBRARIES})
  ELSE (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} hdf5_hl hdf5)
  ENDIF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
ENDIF(NEEDHDF5)


