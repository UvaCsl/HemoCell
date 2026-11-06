# Single GoogleTest runner for all tests under tests/ (sources are globbed;
# reconfigure after adding files). tests/main.cpp initialises MPI once.
set(BINARY ${PROJECT_NAME}_test)
file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS LIST_DIRECTORIES false *.h *.cpp)
add_executable(${BINARY} ${TEST_SOURCES})
target_link_libraries(${BINARY} PRIVATE ${PROJECT_NAME} GTest::gtest)
target_include_directories(${BINARY} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
# GoogleTest >= 1.17 requires C++17; the rest of the project is built as C++11.
# The runner stays in its build tree (unlike case executables), so several
# build trees (e.g. presets) do not overwrite each other's test binary.
set_target_properties(${BINARY} PROPERTIES
  CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
hemocell_enable_warnings(${BINARY})

# CTest entries. The tests read their inputs (tests/validation/..., examples/...)
# relative to the repository root, so they run there (also when invoked
# directly: `<build>/hemocell_test` from the repository root). Performance.* expects a
# configuration argument and hardware-dependent timings and is not registered.
#   ctest --test-dir <build> -L unit         # fast checks
#   ctest --test-dir <build> -L validation   # physics validation (minutes)
add_test(NAME unit
  COMMAND ${BINARY} --gtest_filter=-*Validation*:Performance.*)
add_test(NAME validation
  COMMAND ${BINARY} --gtest_filter=*Validation*)
set_tests_properties(unit PROPERTIES LABELS "unit")
set_tests_properties(validation PROPERTIES LABELS "validation")
set(_hemocell_tests unit validation)

# Communication, ownership and periodicity need more than one rank.
if(MPIEXEC_EXECUTABLE)
  set(HEMOCELL_TEST_MPI_RANKS 2 CACHE STRING "Number of MPI ranks for the *_mpi tests")
  # Open MPI derives its prefix from the launcher path; through a /bin -> /usr/bin
  # symlink (merged-/usr systems) that prefix becomes "/", which it rejects.
  get_filename_component(_hemocell_mpiexec "${MPIEXEC_EXECUTABLE}" REALPATH)
  add_test(NAME unit_mpi
    COMMAND ${_hemocell_mpiexec} ${MPIEXEC_NUMPROC_FLAG} ${HEMOCELL_TEST_MPI_RANKS} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:${BINARY}> ${MPIEXEC_POSTFLAGS} --gtest_filter=-*Validation*:Performance.*)
  add_test(NAME validation_mpi
    COMMAND ${_hemocell_mpiexec} ${MPIEXEC_NUMPROC_FLAG} ${HEMOCELL_TEST_MPI_RANKS} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:${BINARY}> ${MPIEXEC_POSTFLAGS} --gtest_filter=*Validation*)
  set_tests_properties(unit_mpi PROPERTIES LABELS "unit;mpi" PROCESSORS ${HEMOCELL_TEST_MPI_RANKS})
  set_tests_properties(validation_mpi PROPERTIES LABELS "validation;mpi" PROCESSORS ${HEMOCELL_TEST_MPI_RANKS})
  list(APPEND _hemocell_tests unit_mpi validation_mpi)
endif()
set_tests_properties(${_hemocell_tests} PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")

# `--target check`: build the test runner and run the unit tests verbosely.
add_custom_target(check
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L unit
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
  DEPENDS ${BINARY})
