# Build options and compiler flags of HemoCell.
#
# Optimisation follows CMAKE_BUILD_TYPE (Release: -O3 -DNDEBUG, the default).
# The options below add to it; they apply to every target in the project
# (Palabos, HemoCell, cases and tests), so that templates instantiated in case
# code are compiled with the same settings as the library.

option(HEMOCELL_NATIVE_ARCH
  "Optimise for the CPU of the build machine (-march=native). Turn off for binaries that must run on other CPUs."
  ON)
option(HEMOCELL_DEBUG_INFO "Add debug symbols (-ggdb) to non-Debug builds" ON)
option(HEMOCELL_PALABOS_DEBUG "Enable Palabos run-time checks (PLB_DEBUG); slow" OFF)
option(HEMOCELL_WERROR "Treat warnings in HemoCell code as errors" OFF)
set(HEMOCELL_SANITIZERS "" CACHE STRING
  "Sanitizers to enable (GCC/Clang), e.g. 'address;undefined'. Use with RelWithDebInfo or Debug.")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  if(HEMOCELL_NATIVE_ARCH)
    add_compile_options(-march=native)
  endif()
  if(HEMOCELL_DEBUG_INFO)
    add_compile_options($<$<NOT:$<CONFIG:Debug>>:-ggdb>)
  endif()
  if(HEMOCELL_SANITIZERS)
    string(REPLACE ";" "," _hemocell_sanitizers "${HEMOCELL_SANITIZERS}")
    add_compile_options(-fsanitize=${_hemocell_sanitizers} -fno-omit-frame-pointer)
    add_link_options(-fsanitize=${_hemocell_sanitizers})
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  if(HEMOCELL_NATIVE_ARCH)
    add_compile_options(-xHost)
  endif()
  if(HEMOCELL_DEBUG_INFO)
    add_compile_options($<$<NOT:$<CONFIG:Debug>>:-g>)
  endif()
  if(HEMOCELL_SANITIZERS)
    message(WARNING "HEMOCELL_SANITIZERS is only supported for GCC and Clang")
  endif()
endif()

# Warnings for HemoCell's own code (not for the vendored Palabos sources).
set(HEMOCELL_WARNING_FLAGS)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  list(APPEND HEMOCELL_WARNING_FLAGS
    -Wall -Wextra -Wformat -Wformat-security -Werror=format-security
    -Wno-unused-parameter -Wno-deprecated-declarations -Wno-unknown-pragmas)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    list(APPEND HEMOCELL_WARNING_FLAGS -Wno-empty-body -Wno-unused-result -Wno-ignored-qualifiers)
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  list(APPEND HEMOCELL_WARNING_FLAGS -wd858)
endif()

# hemocell_enable_warnings(<target>): warning flags (and -Werror if
# HEMOCELL_WERROR) for the library and the tests. Case directories use
# add_compile_options(${HEMOCELL_WARNING_FLAGS}) and are not affected by
# HEMOCELL_WERROR.
function(hemocell_enable_warnings target)
  target_compile_options(${target} PRIVATE ${HEMOCELL_WARNING_FLAGS})
  if(HEMOCELL_WERROR AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Werror)
  endif()
endfunction()
