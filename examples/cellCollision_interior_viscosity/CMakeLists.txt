# executable will have the same name as its directory
get_filename_component(EXEC_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# write the resulting executable in the _current_ directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# add executable from source files in current directory
# this includes all encountered `*.cpp` and `*.h*` files
add_executable_from_source(${EXEC_NAME} ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_options(${EXEC_NAME} PUBLIC -DINTERIOR_VISCOSITY)

# link the executable to `hemocell` and `HDF5` dependencies
target_link_libraries(${EXEC_NAME} ${PROJECT_NAME})
target_link_libraries(${EXEC_NAME} ${HDF5_C_HL_LIBRARIES} ${HDF5_LIBRARIES})
