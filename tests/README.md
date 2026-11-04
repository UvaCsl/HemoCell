# Tests

This directory contains a collection of tests for `HemoCell`. The tests are
implemented using [`GoogleTest`](https://github.com/google/googletest)
([documentation](https://google.github.io/googletest/)).

A brief example how to create simple tests are presented
[here](https://google.github.io/googletest/primer.html#simple-tests). New tests
can be implemented in any `tests/*.cpp` file, which are detected during
compilation.

In contract to _code_ tests, there are _validation_ tests identified using the
test suite `Validation`. This test suite contains a series of tests collected
under `tests/validation`. These integration tests aim to assert previously
validation (and published) behaviour of `HemoCell` and ensure these remain valid
under development.

To build and run tests (from the repository root):

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --target hemocell_test
ctest --test-dir build -L unit --output-on-failure        # fast, 1 and 2 ranks
ctest --test-dir build -L validation --output-on-failure  # minutes, 1 and 2 ranks
```

CTest runs the binary in the repository root, because the tests read their
inputs (`tests/validation/...`, `examples/...`) relative to it. When running
`build/hemocell_test` directly, do so from the repository root as well; use
`--gtest_filter` to select tests, e.g. `--gtest_filter='*Validation*'` (this also
matches the parameterised `Validation/ValidationForceDisplacement` tests, which
`Validation.*` misses). `Performance.*` needs a configuration argument and is
not part of the CTest runs; see `.github/workflows/performance.yml`.
