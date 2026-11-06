# The stent_strut case for the crossing geometries (config_*.xml): the
# executable is named after this directory, its source is stent_strut.cpp.
get_filename_component(EXEC_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# write the resulting executable in the _current_ directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${EXEC_NAME} "${CMAKE_CURRENT_SOURCE_DIR}/stent_strut.cpp")

# link the executable to `hemocell` and `HDF5` dependencies
target_link_libraries(${EXEC_NAME} ${PROJECT_NAME})
target_link_libraries(${EXEC_NAME} ${HDF5_C_HL_LIBRARIES} ${HDF5_LIBRARIES})
