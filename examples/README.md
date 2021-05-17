# Examples

This directory contains example simulations using `hemocell`. Each directory
represents a separate case study.

## Compilation

When compiling `hemocell` without providing a specific compilation target, all
examples will be compiled. Each example registers itself with a `CMakeLists.txt`
file, as described in the next section. To compile an individual, specific
example, the directory name can be passed as a compile target to `CMake`. To
compile an example, run the following code from within the `hemocell/build/`
directory and replace `$example` by the desired example's directory name.

```bash
cmake --build . --parallel $(nproc) --target "$example"
```

## Adding examples

To add example case `$new_example`:

- Copy the template directory: `cp -r examples/template examples/$new_example`.
- Add the example directory to `examples/CMakeLists.txt` by appending a line
  containing `add_subdirectory("$new_example")`.
- Include your code that defines the `main` in `$new_example/$new_example.cpp`.
  This file is automatically included for compilation by the template
  `CMakeLists.txt`. If you would like to compile additional `*.cpp` or `*.h`
  files, make sure to register them within `$new_example/CMakeLists.txt.`
- Ensure the right variant of the available `hemocell` libraries is linked
  against. By default the standard library is used, i.e. `libhemocell.a`.
  However, when requiring specific compilation features, such as
  `interior viscosity`, `solidification mechanics`, or `parmetis`, the
  corresponding library needs to be selected. If this applies, change the last
  argument `${PROJECT_NAME}` in following line
  `target_link_libraries(${EXEC_NAME} ${PROJECT_NAME})` to one of the following
  - `${PROJECT_NAME}_interior_viscosity`
  - `${PROJECT_NAME}_solidify_mechanics`
  - `${PROJECT_NAME}_parmetis`
- Build `hemocell` as usual or provide `--target "$new_example$"` to only
  compile the newly added example.
