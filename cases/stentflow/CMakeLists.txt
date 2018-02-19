#=======================================
# Build settings. Set these.

# Project setup
SET(PROJECT_NAME stentflow)
SET(PROJECT_SRC "stentflow.cpp")

# HemoCell location relative to CMakelists.txt
SET(HEMOCELL_BASE_DIR "./../../")
get_filename_component(HEMOCELL_BASE_DIR "${HEMOCELL_BASE_DIR}" ABSOLUTE)
SET(HEMOCELL_DIR "${HEMOCELL_BASE_DIR}/build/hemocell")
SET(PALABOS_BASE_DIR "${HEMOCELL_BASE_DIR}/palabos")


MESSAGE( STATUS "HEMOCELL_BASE_DIR:         " ${HEMOCELL_BASE_DIR} )
MESSAGE( STATUS "HEMOCELL_DIR:         " ${HEMOCELL_DIR} )
MESSAGE( STATUS "PALABOS_BASE_DIR:         " ${PALABOS_BASE_DIR} )

#=======================================
##### Beginning of build script

PROJECT(${PROJECT_NAME} CXX C)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_VERBOSE_MAKEFILE 1)
INCLUDE(GNUInstallDirs)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

#=======================================

ADD_DEFINITIONS("-DPLB_MPI_PARALLEL")
ADD_DEFINITIONS("-DPLB_USE_POSIX")
ADD_DEFINITIONS("-DPLB_SMP_PARALLEL")
IF(APPLE)
  ADD_DEFINITIONS("-DPLB_MAC_OS_X")
ENDIF(APPLE)

#=======================================

OPTION(ENABLE_MPI "Enable MPI" ${DEFAULT})
INCLUDE(FindMPI)
IF(MPI_CXX_FOUND)
  SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
ELSE(MPI_CXX_FOUND)
  MESSAGE(FATAL ERROR "MPI compiler not found!")
ENDIF(MPI_CXX_FOUND)

#=======================================

execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version 
                COMMAND head -n1
                COMMAND "cut" "-d " "-f1"
                OUTPUT_VARIABLE CXX_COMPILER_NAME
                OUTPUT_STRIP_TRAILING_WHITESPACE)
MESSAGE(STATUS "COMPILER: ${MPI_CXX_COMPILER}, TYPE: ${CXX_COMPILER_NAME}")

#Set up correct flags for compiler
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -Wformat -Wformat-security")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=format-security")
IF(${CXX_COMPILER_NAME} STREQUAL "g++")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-empty-body")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-qualifiers")
ENDIF(${CXX_COMPILER_NAME} STREQUAL "g++")
IF(${CXX_COMPILER_NAME} STREQUAL "icpc")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd858")
ENDIF(${CXX_COMPILER_NAME} STREQUAL "icpc")

#=======================================
INCLUDE_DIRECTORIES(${PALABOS_BASE_DIR}/externalLibraries)
INCLUDE_DIRECTORIES(${PALABOS_BASE_DIR}/src/libraryInterfaces)

INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/helper)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/config)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/core)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/models)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/mechanics)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/external)
INCLUDE_DIRECTORIES(${HEMOCELL_BASE_DIR}/IO)
INCLUDE_DIRECTORIES(${PALABOS_BASE_DIR}/src)

LIST(APPEND SRC_FILES ${PROJECT_SRC})

add_custom_target(hemocell_pre COMMAND cd ${HEMOCELL_DIR} && ${CMAKE_COMMAND} .)

include(ExternalProject) 
ExternalProject_Add("hemocell" PREFIX ${HEMOCELL_DIR} SOURCE_DIR ${HEMOCELL_DIR}
    BINARY_DIR ${HEMOCELL_DIR} INSTALL_COMMAND "")
ExternalProject_Add_Step("hemocell" update_custom COMMAND cd ${HEMOCELL_DIR} 
    && ${CMAKE_COMMAND} . && $(MAKE)  ALWAYS 1)

ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_FILES})
add_dependencies("hemocell" hemocell_pre)
ADD_DEPENDENCIES(${PROJECT_NAME} "hemocell")
target_link_libraries(${PROJECT_NAME} ${HEMOCELL_DIR}/libhemocell.a)

SET(HDF5_PREFER_PARALLEL 1)
FIND_PACKAGE(HDF5 COMPONENTS C HL)
if(NOT ${HDF5_FOUND})
   message(fatal_error "Hdf5 Libraries not found!")
endif(NOT ${HDF5_FOUND})
message(STATUS "HDF5 libs:      ${HDF5_LIBRARIES}")
# This is needed because the static hdf5 libraries dont work
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${HDF5_LIBRARIES})
if(HDF5_C_LIBRARY_hdf5)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${HDF5_C_LIBRARY_hdf5})
endif(HDF5_C_LIBRARY_hdf5)
if(HDF5_C_LIBRARY_hdf5_hl)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${HDF5_C_LIBRARY_hdf5_hl})
endif(HDF5_C_LIBRARY_hdf5_hl)

