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

To build and run tests:

```bash
cd hemocell/build
cmake ..
cmake --build . --target hemocell_test
ctest -V
```

To in/exclude the validation tests, you can specify the `GTEST_FILTER`
environment variable. To run validation tests `GTEST_FILTER="Validation.*"` or
to exclude them `GTEST_FILTER="-Validation.*"` (not the leading `-`), e.g.

```bash
GTEST_FILTER="Validation.*" ctest -V
```
