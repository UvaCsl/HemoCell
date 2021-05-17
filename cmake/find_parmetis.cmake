# This routine attempts to locate the optional dependency to PARMETIS and METIS.
# For both, it searches for its include directories and the libraries to be
# included. We search in the typical paths:
#
# - /usr/local/include
# - /usr/include/
# - /usr/local/lib
# - /usr/lib
#
# Additionally, we hint CMake to search in the `./external` directory that might
# include a local installation of PARMETIS/METIS.
#
# If both PARMETIS and METIS are detected, the variable `PARMETIS_FOUND` is set
# to true. This is done in the parent's scope to allow to toggle the compilation
# behaviour depending on the presence of these libraries. Additionally, we set
# `PARMETIS_DIRS` and `PARAMETIS_LIBS` to the detected library and source paths.
function(FindParmetis)
        find_path(PARMETIS_INCLUDE_DIR parmetis.h
                /usr/local/include
                /usr/include/
                HINTS external/parmetis/build/local
                HINTS ${PARMETIS_DIR} ENV PARMETIS_DIR
                PATH_SUFFIXES include
        )

        find_library(PARMETIS_LIBRARY libparmetis.a
                /usr/local/lib
                /usr/lib
                HINTS external/parmetis/build/local
                HINTS ${PARMETIS_LIBRARY_DIR} ENV PARMETIS_LIBRARY_DIR
                PATH_SUFFIXES lib
        )

        find_path(METIS_INCLUDE_DIR metis.h
                /usr/local/include
                /usr/include/
                HINTS external/parmetis/metis/build/local
                HINTS ${METIS_DIR} ENV METIS_DIR
                PATH_SUFFIXES include
        )

        find_library(METIS_LIBRARY libmetis.a
                /usr/local/lib
                /usr/lib
                HINTS external/parmetis/metis/build/local
                HINTS ${METIS_LIBRARY_DIR} ENV METIS_LIBRARY_DIR
                PATH_SUFFIXES lib
        )

        message(STATUS "parmetis: " ${PARMETIS_INCLUDE_DIR} )
        message(STATUS "parmetis: " ${PARMETIS_LIBRARY} )
        message(STATUS "metis: " ${METIS_INCLUDE_DIR} )
        message(STATUS "metis: " ${METIS_LIBRARY} )

        if(PARMETIS_INCLUDE_DIR AND PARMETIS_LIBRARY)
                IF(METIS_INCLUDE_DIR AND METIS_LIBRARY)
                        set(PARMETIS_FOUND YES PARENT_SCOPE)
                        set(PARMETIS_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR} PARENT_SCOPE)
                        set(PARMETIS_LIBS ${PARMETIS_LIBRARY} ${METIS_LIBRARY} PARENT_SCOPE)
                endif()
        endif()
endfunction(FindParmetis)
