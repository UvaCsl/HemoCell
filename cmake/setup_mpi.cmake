# Finds MPI package. The directories and libraries are then attached to the
# provided `TARGET`. Note: also the `CMAKE_CXX_COMPILER` is updated to the
# discovered `MPI_CXX_COMPILER`.
function(ConfigureMPI TARGET)
        find_package(MPI REQUIRED)

        set(MPI_FOUND YES PARENT_SCOPE)
        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} PARENT_SCOPE)

        target_include_directories(${TARGET} PRIVATE ${MPI_CXX_INCLUDE_DIRS})
        target_link_libraries(${TARGET} PRIVATE MPI::MPI_CXX)
endfunction(ConfigureMPI)
