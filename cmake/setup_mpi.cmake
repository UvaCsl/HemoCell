# Finds MPI package. The directories and libraries are then attached to each
# target provided in `TARGETS`. Note: aslo the `CMAKE_CXX_COMPILER` is updated
# to the discovered `MPI_CXX_COMPILER`.
function(ConfigureMPI TARGETS)
        find_package(MPI REQUIRED)

        set(MPI_FOUND YES PARENT_SCOPE)
        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} PARENT_SCOPE)

        foreach(TARGET ${TARGETS})
                target_include_directories(${TARGET} PRIVATE ${MPI_INCLUDE_DIRS})
                target_link_libraries(${TARGET} PRIVATE ${MPI_LIBRARIES})
        endforeach()
endfunction(ConfigureMPI)
