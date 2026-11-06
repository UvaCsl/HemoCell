# Each subdirectory builds one executable named after it (see its CMakeLists.txt).
# `--target hemocell_cases` builds all research cases.
add_compile_options(${HEMOCELL_WARNING_FLAGS})

set(HEMOCELL_CASES
  "atherosclerosis"
  "AR2"
  "AR2_pulsatile"
  "AR2_stiff"
  "cellCollision"
  "cellCollision_sphere"
  "CEPAC"
  "fc_solidified_stenosis"
  "flowchamber_stenosis"
  "injured_vessel"
  "kolmogorovFlow"
  "leesEdwards"
  "microvessel_bended"
  "performance_testing"
  "pipeflowMalaria"
  "pipeflow_interior_viscosity"
  "pipeflow_wbc"
  "preinlet_shear"
  "solidify_example"
  "stenosis"
  "stentflow"
  "stent_strut"
  "stent_strut/stent_strut_crossings"
  "stl_preinlet"
  "stretchCell_interior_viscosity"
  "stretchMalaria"
  "unbounded"
  "vasoconstriction_pipe"
)
# The target of each directory is named after its last component.
set(HEMOCELL_CASE_TARGETS "")
foreach(dir IN LISTS HEMOCELL_CASES)
  add_subdirectory("${dir}")
  get_filename_component(target "${dir}" NAME)
  list(APPEND HEMOCELL_CASE_TARGETS "${target}")
endforeach()

add_custom_target(hemocell_cases)
add_dependencies(hemocell_cases ${HEMOCELL_CASE_TARGETS})
