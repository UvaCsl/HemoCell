# This configures `GoogleTest` as testing framework. The function is only
# evaluated when `BUILD_TESTING` is true. In case the framework is not
# present, the source code is automatically pulled from their repository.
#
# For details on the testing framework see:
# see https://github.com/google/googletest/tree/v1.10.x/googletest
#
# In addition to `make test` we define a custom command `make check` that
# calls the tests with more verbose output, e.g. to clarify why tests might
# be failing.
function(ConfigureGTest)
    if(NOT BUILD_TESTING)
        return()
    endif()

    set(git_tag release-1.10.0)
    message(STATUS "Setting up `googletests` ${git_tag}...")

    # Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE BOOL "")

    # modern versions `FetchContent` seems most straightforward
    if(CMAKE_VERSION VERSION_GREATER 3.11)
        include(FetchContent)

        FetchContent_Declare(googletest
            GIT_REPOSITORY      https://github.com/google/googletest.git
            GIT_TAG             ${git_tag})

        FetchContent_GetProperties(googletest)

        if(NOT googletest_POPULATED)
            FetchContent_Populate(googletest)
            add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR}
                             EXCLUDE_FROM_ALL)
        endif()
    else()
        configure_file(cmake/googletest_template.cmake googletest-download/CMakeLists.txt)

        # `cmake ..`
        execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
            OUTPUT_QUIET)

        if(result)
            message(FATAL_ERROR "CMake step for googletest failed: ${result}")
        endif()

        # `cmake --build ...`
        execute_process(COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
            OUTPUT_QUIET)

        if(result)
            message(FATAL_ERROR "Build step for googletest failed: ${result}")
        endif()

        # Add googletest directly to our build. This defines
        # the gtest and gtest_main targets.
        add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
            ${CMAKE_CURRENT_BINARY_DIR}/googletest-build EXCLUDE_FROM_ALL)
    endif()

    unset(CMAKE_SUPPRESS_DEVELOPER_WARNINGS)

    # Adds a custom target `make check` that invokes the tests with the
    # additional flag `--output-on-failure` that prints the output to the
    # console for any tests that might fail.
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
        --force-new-ctest-process --output-on-failure)

endfunction(ConfigureGTest)
