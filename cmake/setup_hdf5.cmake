# Finds HDF5 package and components `C` and `HL`. The directories and
# libraries are then attached to the provided `TARGET`. The routine sets the
# preference for parallel HDF5 (`HDF5_PREFER_PARALLEL`) and static libraries
# (`HDF5_USE_STATIC_LIBRARIES`) before attempting to find the packages.
function(ConfigureHDF5 TARGET)
        # define preferences
        set(HDF5_PREFER_PARALLEL 1)
        set(HDF5_USE_STATIC_LIBRARIES)

        # obtain the package and its component
        find_package(HDF5 REQUIRED COMPONENTS C HL)
        set(HDF5_FOUND YES PARENT_SCOPE)

        # link to target
        target_include_directories(${TARGET} PRIVATE ${HDF5_INCLUDE_DIRS})
        target_link_libraries(${TARGET} PRIVATE ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
endfunction(ConfigureHDF5)
