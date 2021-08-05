# create a test executable for each source file in `tests/`
set (BINARY ${CMAKE_PROJECT_NAME}_test)
file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)
set(SOURCES ${TEST_SOURCES})

# define a test executable per discovered source file
add_executable(${BINARY} ${TEST_SOURCES})

# tests to be runnable with `make test`
add_test(NAME ${BINARY} COMMAND ${BINARY})

# Ensure tests are evaluated from within `tests/` to allow to specify files used
# by tests with relatives paths.
set_tests_properties(${BINARY} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# link with static `hemocell` library, test framework
target_link_libraries(${BINARY} "${PROJECT_NAME}" gtest)
target_link_libraries(${BINARY} ${HDF5_C_HL_LIBRARIES} ${HDF5_LIBRARIES})
