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

The simulations of the tests write their output (HDF5, logs, checkpoints) into
one temporary directory per run (`hemo_test::outputRoot()` in
`test_support.h`), which is removed at the end; set `HEMOCELL_KEEP_TEST_OUTPUT=1`
to keep it. Construct a simulation from `hemo_test::testConfig(path)` (or
`pipeflowConfig()`), never from a configuration that writes into the repository.

`scripts/ci/smoke_drivers.py` runs every built example and case for 20 steps
with output and a checkpoint, in a temporary copy of its directory (`-n` ranks,
default 2; `--slow` adds the large drivers listed in `SLOW`; name drivers such
as `examples/pipeflow` to run only those). CI runs it in the gcc release job.
