# HemoCell development guide

## Overview and documentation

HemoCell simulates dense suspensions of deformable capsules, especially blood
cells. Palabos supplies the lattice Boltzmann fluid solver (LBM); HemoCell tracks
triangulated cell membranes as Lagrangian surface particles and couples them to
the fluid with the immersed boundary method (IBM). Constitutive models compute
membrane forces; MPI distributes fluid blocks and their particles across ranks.

Start with `doc/README.md`, then `doc/user_guide/QuickStart.rst`, `Case.rst`,
`concepts/`, `xml_files.rst`, and `reference/`. For numerical problems, read
`common_mistakes.rst` and `FAQ.rst`. Check the relevant implementation before
relying on an older example or documentation claim. Some documentation describes
historical releases; the checked-in CMake files and source determine current
build behavior. Keep this guide current when that behavior changes.

## Environment and build

Use an out-of-source CMake build. Required components are a C/C++ toolchain,
MPI (typically OpenMPI), HDF5 with its C and high-level libraries, GNU Patch,
and the HemoCell-patched Palabos source at `palabos/`. CMake prefers parallel
HDF5. ParMETIS/METIS is optional, for load balancing. Python supports documentation,
preprocessing and analysis; it is not the simulation runtime.

Library code currently builds as **C++11**. The GoogleTest runner builds as
**C++17** for GoogleTest v1.17.0. Do not introduce C++14/17 features into library
or example code without deliberately updating the build and supported toolchains.

Dependency setup, only when Palabos needs to be installed or replaced:

```bash
# From the repository root; inspect an existing palabos/ before running.
./setup.sh
```

`setup.sh` deletes an existing `palabos/`, downloads the revision pinned in the
script, and applies `patch/palabos.patch` through `patch/patchPLB.sh`. The guide
mentions Palabos v2.3.0, but the script currently pins commit
`05712164d940a42e06afdd705249912fa0c49f14`. Follow the script for this checkout;
do not substitute an arbitrary upstream release or rerun setup as routine cleanup.
Maintain integration changes in `patch/`, rather than only editing the ignored
dependency tree. Avoid unrelated edits to vendored `external/` code.

```bash
# From the repository root, after dependencies are available:
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --target hemocell --parallel 2
cmake --build build --target pipeflow stretchCell --parallel 2
cmake --build build --target hemocell_test --parallel 2
# or with presets (CMake >= 3.21), building into build-<preset>/:
cmake --preset release            # also: portable, debug, sanitize
cmake --build --preset release    # builds hemocell_test (and the library)
```

Machine-specific settings (e.g. `HDF5_ROOT`, MPI header paths) belong in a
`CMakeUserPresets.json` (ignored by git) that inherits from these presets.
Choose parallelism to fit available memory: Palabos is compiled once (as an
object library archived into every HemoCell library variant), but it is large.
Initial configuration with `BUILD_TESTING=ON` fetches GoogleTest (pinned
v1.17.0) unless an installed GoogleTest >= 1.17 is found. `-DBUILD_TESTING=OFF`
skips that dependency and the tests; it is useful for a library-only build, not
a substitute for testing a code change.

Build details that matter:

- The default build compiles `hemocell` only. Examples, cases, tests and optional
  library variants are excluded from `all`; request their targets explicitly
  (`hemocell_examples` and `hemocell_cases` build all examples/cases).
- Example/case targets match their directory names and put executables into
  those source directories (so the last build tree built wins). `hemocell_test`
  is written to the build directory.
- Library `.cpp` files and test sources are discovered recursively with CMake
  globs (`CONFIGURE_DEPENDS`). Additional source files for an example must be
  registered explicitly in its `add_executable` call, and the directory added
  to the list in `examples/CMakeLists.txt` or `cases/CMakeLists.txt`.
- The build type is honoured (default `Release`: `-O3 -DNDEBUG`). Options:
  `HEMOCELL_NATIVE_ARCH` (default ON, `-march=native`; binaries may not run on
  other CPUs, use the `portable` preset for those), `HEMOCELL_DEBUG_INFO`
  (`-ggdb`, default ON), `HEMOCELL_PALABOS_DEBUG` (`PLB_DEBUG` run-time checks),
  `HEMOCELL_SANITIZERS` (e.g. `address;undefined`) and `HEMOCELL_WERROR`
  (library and tests only). These flags apply to all targets; warnings are
  only enabled for HemoCell code, not for Palabos or vendored tinyxml2.
- MPI and HDF5 are linked through `MPI::MPI_CXX` and an interface target that
  keeps MPI's include directories before HDF5's: an HDF5 prefix such as a conda
  environment may ship another MPI's `mpi.h`, which must not shadow the MPI the
  code is built with.

### Compile-time variants

| Library target | CMake feature definition |
|---|---|
| `hemocell` | Base target |
| `hemocell_interior_viscosity` | `INTERIOR_VISCOSITY` |
| `hemocell_solidify_mechanics` | `SOLIDIFY_MECHANICS` |
| `hemocell_parmetis` | `HEMO_PARMETIS`; available only when ParMETIS is found |

Select the appropriate library in the case's `target_link_libraries`, following
`examples/cellCollision_interior_viscosity/CMakeLists.txt` or
`cases/solidify_example/CMakeLists.txt`. Feature definitions are public so clients
and libraries agree on class layouts. Do not enable a layout-changing macro only
in an executable or link incompatible variants together. Combined features need
an explicitly consistent build; the existing targets do not represent every
combination.

Also inspect `config/constant_defaults.h`: it currently defines
`SOLIDIFY_MECHANICS` even for the base target. Target names alone therefore do not
tell you all enabled macros. Runtime switches in XML/per-cell setup are separate
from compiled support. Changes to these defaults affect every target and can
change particle layout, performance and checkpoint compatibility.

## Source layout and responsibilities

| Location | Responsibility |
|---|---|
| `hemocell.h`, `core/hemoCell.cpp` | Public `hemo::HemoCell` interface, initialization, stepping, outputs, checkpoints and case-level controls |
| `core/hemoCellFields.*` | All cell types, shared multi-block particle field, communication orchestration and checkpoint save/load |
| `core/hemoCellField.*` | One cell type: material configuration, reference mesh, mechanical model, output selection and update interval |
| `core/hemoCellParticle.h` | Surface vertex state, serialized values, force pointers and particle integration |
| `core/hemoCellParticleField.*` | Atomic-block particles, spatial/cell lookup structures, force evaluation, repulsion and local coupling |
| `core/hemoCellParticleDataTransfer.*` | Particle serialization, MPI/local transfers and periodic position/cell-ID offsets |
| `core/hemoCellFunctional.*`, `core/immersedBoundaryMethod.h` | Palabos processing functionals and IBM interpolation kernels |
| `mechanics/` | `CellMechanics`, RBC/PLT/WBC/malaria force laws, reference cell constants, SI–LBM conversion in `Parameters` |
| `config/` | XML access, compile-time constants, global runtime settings and logging |
| `helper/` | Arrays/geometry, mesh generation/metrics, cell/fluid/particle statistics, voxelization, stretching, pre-inlets, load balancing, profiling and optional physics helpers |
| `io/` | Initial positions, HDF5 fields, CSV summaries and output field extraction |
| `examples/` | User-facing simulation drivers and their XML, position and geometry inputs |
| `cases/` | Research cases, including experimental or incompletely validated features; validate before treating them as reference behavior |
| `tests/` | GoogleTest runner, validation cases and performance-related inputs |
| `tools/` | Standalone `packCells`, cell material templates/tester and position conversion tools |
| `scripts/` | Postprocessing, analysis, cluster environments and batch templates; older Python 2 variants also exist |
| `patch/`, `cmake/`, `external/` | Palabos integration patch, build/dependency helpers and vendored dependencies |
| `doc/`, `legacy/` | Sphinx/Doxygen sources and historical implementations respectively |

Keep reusable simulation behavior in the library. Case drivers should assemble
geometry, boundary conditions, materials, forcing and run/output schedules.
Keep postprocessing in `scripts/` or `tools/`; do not move solver behavior into
analysis scripts. `legacy/` is a historical reference, not the default extension
point.

## Case workflow and configuration

There is no single solver CLI: each case is a C++ executable using `HemoCell`.
Build from the repository root, then **run from the case's input directory**:

```bash
cmake --build build --target pipeflow --parallel 2
cd examples/pipeflow
mpirun -n 4 ./pipeflow config.xml
# Resume from a checkpoint produced by this case (substitute its actual path):
mpirun -n 4 ./pipeflow output/checkpoint/checkpoint.xml
```

The checked-in pipeflow config requests `output/`; absent an explicit
`<parameters><outputDirectory>`, the default is `tmp/`. Existing output
directories receive suffixes such as `_0`. Use the actual run directory for
restart/analysis. Paths to `.pos`, cell XML and geometry are generally relative
to the working directory, not to the main XML file.

Typical initialization order (see `examples/pipeflow/pipeflow.cpp`):

1. Construct `HemoCell` with the configuration and command-line arguments.
2. Compute `param::lbm_*_parameters` and print them. Create the domain/lattice,
   boundary conditions and equilibrium, then initialize the lattice.
3. Call `initializeCellfield()`. Register each type with
   `addCellType<Mechanics>(name, shape)`, configure outputs and update intervals,
   and set periodicity through `setSystemPeriodicity()` so both fields agree.
4. Set initial clearance controls before `loadParticles()`. On a restart follow
   the case's `cfg->checkpointed`/`loadCheckPoint()` branch instead of loading
   fresh positions. Enable interactions before stepping.
5. Run `iterate()` with case-specific forcing and output/checkpoint intervals.

`iterate()` spreads particle forces, advances the fluid, interpolates velocity
and synchronizes particles at the configured interval, advances particles,
updates constitutive forces and optional physics, cleans nonlocal particles,
clears the lattice force and increments `iter`. Preserve this ordering when
changing coupling code. **External fluid driving force must be reapplied after
each `iterate()`**; also set it before the first driven step when required.

To add a driver, copy `examples/template/`, supply `<directory-name>.cpp` and all
required inputs, and register it with `add_subdirectory` in
`examples/CMakeLists.txt` or `cases/CMakeLists.txt`. The template contains the
build skeleton, not a complete runnable simulation. Add the generated executable's
source-tree path to the root `.gitignore` when adding a target.

Configuration conventions:

- XML has a `<hemocell>` root. Access values with
  `(*hemocell.cfg)["domain"]["dx"].read<T>()`; tag names are case-sensitive.
  Cell files use `<MaterialModel>` (capitalization matters).
- Register cells with, for example,
  `addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE)`. The name determines
  `RBC.xml` and `RBC.pos`. Use the same identifier in subsequent per-type calls.
- A `.pos` file starts with the cell count, followed by one row per cell:
  `x y z rotation_x rotation_y rotation_z`. Positions are in **micrometres**;
  rotations are in **degrees**, applied in X, Y, Z order by the loader.
- Settings such as `tmax`, `tmeas`, `tcsv`, `tcheckpoint` and warmup are consumed
  by individual drivers; adding an XML tag alone does not activate behavior.
- Missing XML children throw `std::invalid_argument`. `XMLElement::read` uses
  stream extraction and does not provide comprehensive value validation. Check
  new required values and cross-field constraints before entering collective or
  numerical work; optional values need explicit defaults.
- Preserve existing key names and public API spellings, including
  `TimeScaleSeperation` where used. For configuration/API changes, update readers,
  affected example/test inputs and `doc/user_guide/` together. Explain any
  compatibility impact rather than silently changing old input semantics.

## Scientific and parallel correctness

- The solver uses lattice units internally. Use `mechanics/constantConversion.*`
  (`hemo::Parameters`, aliased as `param`) for `dx`, `dt`, `dm`, `df` and related
  conversions. Document each API's units explicitly: XML mixes SI values,
  micrometres, dimensionless coefficients and `kBT`-scaled material constants.
- `dx` is metres per lattice spacing and `dt` is seconds per iteration. A
  negative input `dt` selects `tau = 1` and computes `dt`. Check low lattice
  velocity (`u_lbm < 0.1`) and `tau > 0.5`; these are stability guidelines, not
  proof of physical validity. Never conceal a regression by retuning material
  constants, force caps or resolution.
- `FORCE_LIMIT` is specified in pN per surface point and converted to lattice
  force in `Parameters`. Check the actual use site for new terms rather than
  assuming every model applies the same scaling. Repulsion constants supplied
  to `setRepulsion`/`enableBoundaryParticles` are pN and their cutoffs
  micrometres, both converted there; the repulsion force per vertex is capped
  at `FORCE_LIMIT` separately from the membrane force.
- Update intervals must be positive. The material intervals and the enabled
  repulsion/interior-viscosity intervals checked in `HemoCell::sanityCheck()`
  must be divisible by the particle-velocity update interval. Preserve these
  synchronization constraints and assess numerical effects when changing them.
- Preserve closed mesh topology, vertex/triangle indexing, normal orientation,
  reference lengths/areas/volume and force balance when changing mechanics.
  Add models through `CellMechanics` (`ParticleMechanics`, `statistics`) and
  register them with `addCellType<Model>`. Start from the nearest existing model
  and include a small case or validation for the new force response.
- Particle force-component pointers normally alias `sv.force`; output separates
  them and later reunifies them. Respect that lifecycle and particle copy/repoint
  behavior. Do not assume each force pointer owns a distinct allocation.
- Particle storage changes can invalidate pointers, references and lookup
  caches. Inspect add/remove/copy and communication paths before keeping a
  particle pointer across container mutation or synchronization.
- `serializeValues_t` is transferred as raw bytes. Changing member types,
  ordering, precision or feature macros affects MPI transfer and potentially
  checkpoint compatibility. Review send/receive, local/periodic/pre-inlet paths
  and save/load together; do not serialize owning pointers or dynamic containers
  by copying their object bytes.
- Distinguish absolute coordinates from atomic-block-local coordinates and real
  owned particles from halo copies. Preserve periodic cell-ID mapping and
  center-based ownership; avoid double-counting ghost cells in statistics.
  Particle envelopes must cover stretched cells and coupling support. A single
  rank run cannot establish communication correctness.
- Keep MPI collectives in matching order on all participating ranks. A routine
  that gathers data or writes coordinated output may need to be called by every
  rank even if only rank zero writes. Use `pcout`/`hlog` and existing profiler
  conventions; avoid per-particle logging and extra barriers in hot loops.
  Report errors with `hemo::fatal()` (`config/error.h`; prints with the rank on
  any rank, ends the job with `MPI_Abort`), never `exit()`: `hlog` prints on
  rank 0 only.
- Runs with the same rank count and decomposition are bit-for-bit
  reproducible, checkpoint files included (received particles are unpacked by
  source rank, not arrival order). Keep received data processed in a fixed
  order, and compare outputs exactly when a change should not alter numerics.
  Do not use Palabos' `duplicateOverlaps(modif::hemocell)` for the particle
  field: it unpacks in arrival order.
- Output and statistics must not change the state (currently they do not, bit
  for bit): read complete cells inside a `HemoCellFields::CompleteCellsScope`
  (a copy of the particles, discarded afterwards), never by calling
  `syncEnvelopes()` or `deleteIncompleteCells()` directly; `iterate()` deletes
  broken cells.
- Cell IDs are `plint` and grow with every periodic crossing (never folded
  back); do not store them in an `int`.
- Global/static simulation state means only one `HemoCell` instance may exist
  at a time. Do not assume reentrancy or thread safety. `~HemoCell` resets all
  of it (`resetGlobalState()` in `core/hemoCell.cpp`), so sequential
  simulations in one process start alike; reset any new static state there,
  or keep per-simulation data in an object.

## Outputs and supporting tools

`setOutputs(name, {OUTPUT_*...})` selects cell fields;
`setFluidOutputs({OUTPUT_*...})` selects fluid fields. Constants live in
`config/constant_defaults.h`. Output is normally in SI units;
`outputInSiUnits = false` requests lattice units. Check dataset-level precision
when selecting tolerances: double-precision computation does not imply every
HDF5 dataset is stored as doubles.

Runs write `hdf5/<iteration>/`, XDMF `*.xmf` descriptors, `csv/`, logs and
`checkpoint/` under the selected output directory. `writeOutput()` includes
HDF5, XDMF and CSV; drivers can call
`writeCellInfo_CSV` separately for a different CSV frequency. Checkpoints include
lattice/particle `.dat` files and, written last, `checkpoint.xml`, with previous
copies renamed `.old`; the XML records a build fingerprint (particle layout,
cell types, dx, dt) that a restart checks.
Restarting may use a different number of ranks, but requires compatible state
and a driver that implements restart handling.

HDF5 files are written per atomic block; every cell is written once (by the
block holding its vertex 0) and fluid files hold bulk nodes only. Rank 0 writes
`<name>.<iteration>.xmf` per step and `<name>.xmf` with all steps (XDMF 3,
readable by ParaView's XDMF3 and legacy readers); `scripts/batchPostProcess.sh`
(`scripts/HDF5toXDMF.py`) writes the same for output of older versions. When
changing output constants, dataset names/shapes, units or CSV columns, update
writers, the XDMF records (`io/xdmfOutput.h`), `scripts/HDF5toXDMF.py` and
`doc/user_guide/reference/outputs.rst` together. Do not commit generated
results, binaries, build trees or dependency downloads.

`tools/packCells` has an independent CMake build and writes its binary into that
directory. It produces `.pos` inputs with domain dimensions in micrometres:

```bash
cmake -S tools/packCells -B tools/packCells/build
cmake --build tools/packCells/build --parallel 2
tools/packCells/packCells --help
```

## Testing and verification

Build `hemocell_test` explicitly. All discovered tests share `tests/main.cpp`,
which initializes MPI once. Construct simulations in tests with
`HemoCell::MPIHandle::External`; do not initialize/finalize MPI for each test.
Use meaningful assertions on behavior, not just successful execution.

Run the tests with CTest, or run the binary **from the repository root** (the
tests read their inputs relative to it):

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --target hemocell_test --parallel 2
ctest --test-dir build -L unit --output-on-failure        # unit, 1 and 2 ranks
ctest --test-dir build -L validation --output-on-failure  # validation, 1 and 2 ranks
# or directly:
build/hemocell_test --gtest_list_tests
build/hemocell_test --gtest_filter='-*Validation*:Performance.*'
build/hemocell_test --gtest_filter='*Validation*'
# For changes to communication, ownership, periodicity or decomposition:
mpirun -n 2 build/hemocell_test --gtest_filter='*Validation*'
```

CTest registers `unit`, `validation`, `unit_mpi` and `validation_mpi` (labels
`unit`, `validation`, `mpi`; rank count `HEMOCELL_TEST_MPI_RANKS`, default 2),
all running in the repository root.

The broad `*Validation*` pattern includes parameterized stretch tests named
`Validation/ValidationForceDisplacement...` as well as the `Validation` suite.
The documented `Validation.*` pattern misses those parameterized tests.
Exclude `Performance.*` from ordinary test runs: `Performance.cube` expects a
configuration as its first positional argument and asserts hardware-dependent
timings. Run it deliberately with a matching rank count and configuration, e.g.
`mpirun -n 2 build/hemocell_test tests/performance/configs/config-2.xml --gtest_filter='Performance.*'`
(GitHub workflow `performance.yml`, manual, self-hosted runner).
The unit tests cover selected components (particle grid and repulsion, config
reader, arrays, particle copies, logging, checkpoint restart); passing them is
not sufficient verification for solver changes.

Before submitting larger changes, also build and run the unit tests with the
`sanitize` (AddressSanitizer + UBSan; set `ASAN_OPTIONS=detect_leaks=0`) and
`debug` (Palabos assertions) presets. CI (`.github/workflows/ci.yml`) builds
everything with GCC (`-Werror`) and Clang and runs unit and validation tests
on 1 and 2 ranks, plus the unit tests under the sanitizers.

Match verification to the change:

- Mechanics, scaling or integration: relevant stretch/shear/pipeflow validations
  with justified physical tolerances. Do not widen bounds merely to get a pass.
- Communication/coupling: compare a suitable small case on one and multiple
  ranks, including block/periodic crossings and cell counts. Exercise checkpoint
  restart, including a changed rank count, when modifying state or transfer.
- Optional physics: build the affected variant and its example explicitly;
  `hemocell_test` links the base library and does not cover every variant.
- I/O: inspect representative datasets/CSV and exercise affected postprocessing.
  Performance changes: compare the same case/configuration, ranks and workload
  using `hemo::global.statistics`, and verify numerical behavior as well.

Report exactly what was built/run and any missing dependencies or skipped
validation. Do not claim a simulation or build was verified from inspection alone.

## C++ style and code documentation

Formatting is defined by `.clang-format` (clang-format >= 14; `.editorconfig`
for editors); all HemoCell sources follow it. Format the lines you change with
`git clang-format`; CI checks the changed lines of pull requests.
`scripts/format.sh [--check]` formats or checks whole files, but do not
reformat unrelated files in a functional change. Pure reformatting commits are
listed in `.git-blame-ignore-revs` (enable locally with
`git config blame.ignoreRevsFile .git-blame-ignore-revs`). Vendored and
historical code (`external/`, `legacy/`, generated `palabos/`) is excluded.
Follow the surrounding file's naming. Keep edits focused and avoid broad
modernization or unsolicited dependency/toolchain changes. Internals,
conventions, extension recipes and known limitations are documented in
`doc/user_guide/developer/`.

- Keep library symbols, constants included (`OUTPUT_*`, cell shapes), in
  `namespace hemo`; no new macros except compile-time options. Headers contain
  no `using namespace` at all: qualify `std::` and `plb::` there. Library
  `.cpp` files may use `using namespace std; using namespace plb;`; drivers
  and tests declare their own (`hemo`, `plb`, `std`), and then `Array` is
  ambiguous between `hemo::Array` and `plb::Array`. Use include guards and keep template definitions visible in headers
  or included `.hh` files, following the existing layout.
- Use existing `T`, `plint`/`pluint`, `hemo::Array` and Palabos types at their
  interfaces. Check conversions between signed indices, container sizes and
  narrow particle IDs. `hemo::Array` and `plb::Array` are distinct types.
- Prefer const references for read-only large objects, explicit ownership and
  RAII for newly owned resources where compatible with Palabos. Existing raw
  pointers may be owned by `HemoCell`, a field or Palabos; trace destruction
  before changing ownership. Mark new virtual overrides with `override`.
- Avoid repeated allocation, copies and invariant calculations in per-vertex,
  per-triangle and per-step loops. Preserve numerical ordering/precision unless
  a measured, validated change justifies altering it.
- Document public APIs and nontrivial numerical helpers with Doxygen comments
  (`@brief`, `@param[in/out]`, `@return`, and preconditions as applicable).
  Explain units, coordinate frame, array dimensions/index meaning, valid ranges,
  ownership/lifetime, mutation and whether all MPI ranks must participate.
  Describe equations or cite the relevant model/paper for scientific algorithms.
- Explain the reason for numerical guards, caps and synchronization decisions.
  Document real error behavior; do not promise exceptions or validation the
  implementation does not provide. Update docs alongside behavioral changes.

Documentation build (Doxygen and Graphviz must be installed; Python dependencies
are in `doc/requirements.txt`):

```bash
python3 -m venv doc/.venv
source doc/.venv/bin/activate
python -m pip install -r doc/requirements.txt
make -C doc/user_guide html
```

Sphinx invokes Doxygen via `doc/user_guide/conf.py`. Edit source `.rst`/`.md`
pages and C++ comments, not generated `doc/doxygen/{html,xml}` or
`doc/user_guide/_build/` output. Register new pages in the appropriate toctree.
