# Contributing to HemoCell

Thank you for helping to improve HemoCell. The full guide is part of the
documentation: [Contributing](doc/user_guide/community/contributing.rst) and the
[Developer Guide](doc/user_guide/developer/index.rst) (internals, conventions,
extending HemoCell, known limitations). In short:

1. Branch from `dev` (in your fork) and open your pull request against `dev`.
2. Build and test (see `AGENTS.md` for details):
   ```bash
   ./setup.sh                                   # once: download and patch Palabos
   cmake --preset release                       # or: cmake -S . -B build
   cmake --build --preset release               # builds hemocell_test
   ctest --test-dir build-release -L unit       # seconds, 1 and 2 MPI ranks
   ctest --test-dir build-release -L validation # minutes; for physics/MPI changes
   ```
3. Format the lines you changed with `git clang-format` (style: `.clang-format`).
4. Add tests for new behaviour, and describe behaviour changes in
   `CHANGELOG.md` and the user guide.

Pull requests are checked by CI: builds with GCC and Clang, unit and validation
tests on one and two ranks, unit tests under AddressSanitizer/UBSan, formatting
of the changed lines, and the documentation build.

Contributors are listed in [CONTRIB.md](CONTRIB.md).
