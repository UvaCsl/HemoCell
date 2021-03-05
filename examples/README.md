# Examples

This directory contains example simulations using `hemocell`. Each directory
represents a separate case study.

## Compilation

The examples are compiled by default when compiling `hemocell`. To compile a
single, specific example, the directory name can be passed as a compile target
to `CMake`. To compile this example, run the following code from
`hemocell/build/` and replace `"$example"` by the desired example's directory
name.

```bash
cmake --build . --parallel $(nproc) --target "$example"
```

## Adding examples

To add example case `$new_example`:

- Copy the template directory: `cp -r examples/template examples/$new_example`.
- Add the example directory to `examples/CMakeLists.txt` by appending a line
  containing `add_subdirectory("$new_example")`.
- Include your code in `$new_example/$new_example.cpp`.
- Build `hemocell` as usual or provide `--target "$new_example$"` to only
  compile the newly added example.
