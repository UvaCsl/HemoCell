# executable will have the same name as its directory
get_filename_component(EXEC_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# write the resulting executable in the _current_ directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# add executable from source files in current directory assuming
# to add more source files, manually register them using `add_executable`
get_filename_component(MAIN ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_executable(${EXEC_NAME} "${CMAKE_CURRENT_SOURCE_DIR}/${MAIN}.cpp")

# Link the executable to the default `hemocell`. Alternative libraries can be
# linked to include specific library features, such as:
# - `${PROJECT_NAME}`: the default library
# - `${PROJECT_NAME}_interior_viscosity`: enables interior viscosity
# - `${PROJECT_NAME}_solidify_mechanics`: enables solidification mechanics
# - `${PROJECT_NAME}_parmetis`: enables linking with PARMETIS/METIS
target_link_libraries(${EXEC_NAME} ${PROJECT_NAME})

# link the executable with `HDF5`
target_link_libraries(${EXEC_NAME} ${HDF5_C_HL_LIBRARIES} ${HDF5_LIBRARIES})
