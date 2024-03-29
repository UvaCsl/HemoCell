# executable will have the same name as its directory
get_filename_component(EXEC_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# write the resulting executable in the _current_ directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# add executable from source files in current directory assuming
# to add more source files, manually register them using `add_executable`
get_filename_component(MAIN ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_executable(${EXEC_NAME} "${CMAKE_CURRENT_SOURCE_DIR}/${MAIN}.cpp")

# link the executable to `hemocell` and `HDF5` dependencies
target_link_libraries(${EXEC_NAME} "${PROJECT_NAME}")
target_link_libraries(${EXEC_NAME} ${HDF5_C_HL_LIBRARIES} ${HDF5_LIBRARIES})
