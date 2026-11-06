# Each subdirectory builds one executable named after it (see its CMakeLists.txt).
# `--target hemocell_examples` builds all examples.
add_compile_options(${HEMOCELL_WARNING_FLAGS})

set(HEMOCELL_EXAMPLES
  "capillary"
  "cellCollision_interior_viscosity"
  "cell_shapes"
  "cube"
  "curvedflow_with_preinlet"
  "flowaroundsphere"
  "microcontraction"
  "oneCellShear"
  "parachuting"
  "parallelplanes"
  "pipeflow"
  "pipeflow_with_preinlet"
  "simple"
  "stretchCell"
)
foreach(dir IN LISTS HEMOCELL_EXAMPLES)
  add_subdirectory("${dir}")
endforeach()

add_custom_target(hemocell_examples)
add_dependencies(hemocell_examples ${HEMOCELL_EXAMPLES})
