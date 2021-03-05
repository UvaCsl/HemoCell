# Command-line options
#
# These options can set additional arguments/variables for the compilation
# of hemocell and related examples/cases. These options can be provided to
# instruct preprocessing variables. Note, passing external `-D` flags to
# `cmake ..` can collide with the arguments below.

# INTERIOR_VISCOSITY
#
# Enables interior viscosity in hemocell, i.e. enables all code withi
# `#IFDEF INTERIOR_VISCOSITY`. This option comes with a performance hit to
# hemocell and should only be used when interior viscosity is desired.
#
# Values:
# - OFF (default): do not consider interior viscosity
# - ON: do consider interior viscosity
#
# Usage:
# >>> cmake .. `-DINTERIOR_VISCOSITY=ON`
option(INTERIOR_VISCOSITY "Enable interior viscosity in hemocell" OFF)

if(INTERIOR_VISCOSITY)
        add_definitions("-DINTERIOR_VISCOSITY")
endif(INTERIOR_VISCOSITY)

message(STATUS "INTERIOR_VISCOSITY=" ${INTERIOR_VISCOSITY})
