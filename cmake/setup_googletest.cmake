# GoogleTest for the `hemocell_test` runner (only when BUILD_TESTING is ON).
#
# An installed GoogleTest (>= 1.17, it needs C++17) is used when found, which
# helps on clusters without internet access (point CMake to it with
# -DGTest_DIR=...). Otherwise the pinned release is downloaded at configure
# time.
function(ConfigureGTest)
  find_package(GTest 1.17 CONFIG QUIET)
  if(GTest_FOUND)
    message(STATUS "Using installed GoogleTest ${GTest_VERSION}")
    return()
  endif()

  message(STATUS "Fetching GoogleTest v1.17.0")
  include(FetchContent)
  set(_exclude_from_all)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.28)
    set(_exclude_from_all EXCLUDE_FROM_ALL) # only built for hemocell_test
  endif()
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        52eb8108c5bdec04579160ae17225d66034bd723 # v1.17.0
    ${_exclude_from_all})
  # Do not install GoogleTest with HemoCell; keep the parent's runtime settings.
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endfunction()
