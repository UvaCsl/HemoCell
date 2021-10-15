# Cases

This directory contains simulations cases using `hemocell`. Each directory
represents a separate case study. Compared to the simulations given under
`examples/` the ones found here (in `cases/`) are considered more experiment,
potentially highlighting unvalidated and/or unfinished parts of the HemoCell
library. Users should proceed with care when using, or extending, examples from
the `cases/` directory and ensure they perform the necessary validation and
verification of their numerical simulations.

## Compilation

Some cases are compiled by default when compiling `hemocell`. To compile a
single, specific case, the directory name can be passed as a compile target
to `CMake`. To compile this case, run the following code from `hemocell/build/`
and replace `"$case"` by the desired case's directory name.

```bash
cmake --build . --parallel $(nproc) --target "$case"
```

## Adding cases

Adding cases is similar to adding examples (see `../examples/README.md`). To add
a simulation case `$new_case`:

- Copy the template directory: `cp -r ../examples/template ./$new_case`.
- Add the case directory to `cases/CMakeLists.txt` by appending a line
  containing `add_subdirectory("$new_case")`.
- Include your code in `$new_case/$new_case.cpp`.
- Build `hemocell` as usual or provide `--target "$new_case$"` to only
  compile the newly added simulation case.
