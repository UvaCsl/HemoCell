# TODO

This document contains remaining TODOs and/or remarks of the existing
implementation. This mostly refers to implementations that do resolve issues,
but might not have been tested extensively in other configurations, or other
systems.

## Pre-inlet communication

The messages send between the pre-inlet and main domain have been significantly
reduced by performing bounding-box checks between the involved ranks. This makes
sure that only atomic blocks that are overlapping perform communication. This
has only been tested with pre-inlets aligned along the positive `x` direction,
i.e. `Direction::Xpos`. The other configurations have *not* been tested yet.

## Solidification communication

The original solidification of platelets near binding sites was fully
sequential. Although the approach works suitable for single core simulations,
multi-core simulations ran into issues (see issue #86), where the binding sites
were not properly propagated into neighbouring atomic blocks.

A patch is added that inserts an explicit synchronisation moment between
labelling and removal of tagged particles (e.g. cells that were close to a
binding site), such that information near the boundaries of atomic blocks is
shared with their neighbours. Although this resolves the original issues for a
simple case, such as described in `cases/solidify_example`, more extensive tests
are required for large(r) scale problems. In preliminary tests, issues were
still observed if regions were marked as binding site that were already
overlapping with cells, creating broken or torn cells in the simulation domain.

## Hardening before merging `refactor`

Findings of a review of the `refactor` branch (November 2026). Items marked
*(verified)* were reproduced by a run or a compile; the others follow from
reading the code. Suggested order: A1-A3 and A5 first, then A6-A7, then B10-B12
before others write drivers against the API.

Status: A1-A9 and B11-B16 are done (see CHANGELOG). B10 is deferred (special
case, lower priority).

### A. Silent wrong results or lost runs

1. *Done* (`hemo::fatal()`, `hemo::warning()`, `config/error.h`).
   **Fatal errors lose their message and may not end the job.** There are 45
   `exit(1)` calls in library code, most logging through `hlog`, which writes on
   rank 0 only: an error detected on another rank leaves no reason, and
   launchers other than `mpirun` may leave the other ranks hanging. Add one
   `hemo::fatal(msg)` that prints with the rank to stderr and calls `MPI_Abort`
   (or throws in tests, so error paths can be tested).
2. *Done* (non-finite vertex forces are removed with their cell; repulsion
   skips zero distances). **NaN/infinity pass into the fluid.** No `isfinite` check exists; the force
   cap (`core/hemoCellParticleField.cpp`, `spreadParticleForce()`) uses
   `force_mag > f_limit`, which is false for NaN. Link and bending terms divide
   by edge and normal lengths without a guard, and collapsed RBC triangles occur
   at high shear. One bad vertex corrupts the whole lattice silently. Test
   `!(force_mag <= f_limit)` in the existing branch (no cost in the normal path)
   and report and remove the cell.
3. *Done* (`checkpoint.xml` last and atomic; fingerprint checked on load).
   **Checkpoints are not crash-safe or self-describing.** `checkpoint.xml` is
   written before the data files (`HemoCellFields::save()`): a job killed during
   writing leaves new metadata next to partial `.dat` files. Nothing records
   `sizeof(serializeValues_t)`, the compile-time variant, the cell-type names and
   order, or dx/dt, so a restart with another variant or registration order
   misreads the raw particle bytes. Write the XML last (or write to a temporary
   directory and rename it), and validate a fingerprint on load.
4. *Done* (checked once when the cell type is added).
   **Inner-edge indices are not validated.** `InnerEdges` are read without a
   check against the vertex count (`mechanics/commonCellConstants.cpp`); the WBC
   template's pairs belong to the 642-vertex mesh and index out of bounds with
   another `minNumTriangles`. Check once at setup.
5. *Done* (missing `.pos` is fatal; placed cells are reported; cells across a
   periodic boundary are counted, and placed as periodic images on request,
   `setPlaceCellsAcrossPeriodicBoundaries()`). **The cell loader drops cells silently** *(verified)*. A missing `.pos` file
   is only a warning (`microcontraction`, `stentflow`, `stent_strut` ran without
   cells); the number of cells actually placed is never reported (86 of 100 RBCs
   were dropped without notice in `cases/pipeflow_wbc` before its positions were
   regenerated); cells crossing a periodic boundary are never placed. Make a
   missing file an error, log placed versus requested cells per type with the
   reason, and consider placing periodic images.
6. *Done* (output and statistics read a copy of the particles,
   `HemoCellFields::CompleteCellsScope`; `iterate()` deletes broken cells).
   **Output and statistics change the trajectory** *(verified)*. Two identical
   2-rank pipeflow runs diverge (round-off, about 3e-10 relative, then growing)
   when one takes statistics or writes output every 100 steps. `writeOutput()`
   and the `CellInformationFunctionals` helpers call `syncEnvelopes()`, which
   rebuilds the envelopes in a different particle order, and several helpers
   also call `deleteIncompleteCells()`. The CHANGELOG statement "results are now
   independent of the output interval" holds physically, not bit for bit. Put
   envelope particles in a canonical order after a sync, or let output and
   statistics read without re-syncing; statistics must not modify state.
7. *Done* (IDs are `plint` throughout; the HDF5 `Cell Id` dataset is a 64-bit
   integer).
   **Cell IDs grow with every periodic lap.** A crossing adds `N` (x),
   `100 N` (y) or `10000 N` (z) and nothing folds the ID back, while
   `sanityCheck()` and `developer/architecture.rst` assume a bounded shift. With
   200,000 cells the ID overflows `int` after about 640 laps in x; the cell is
   then reset to its base ID and most likely deleted. Use 64-bit IDs throughout
   (`sv.cellId` is already `plint`; the lookup maps and `base_cell_id()` use
   `int`); checkpoints are unaffected.
8. *Done* (partners from all interface blocks of a rank; every block
   assigned to a rank of its lattice).
   **The pre-inlet is wrong with more than one block per rank.**
   `helper/preInlet.cpp` keeps only the last bulk per rank (FIXME there). Fail
   clearly in that case, or support it.
9. *Done* (`k` in pN, 5 pN in the shipped configurations; repulsion capped
   at `FORCE_LIMIT`). **Repulsion is a no-op in all shipped configurations.** `kRep` is used as a
   lattice-unit force (`HemoCell::setRepulsion()`, `core/particleGrid.cpp`), so
   `2e-22` is about 1e-18 pN; repulsion is also added after the `FORCE_LIMIT`
   cap, so a meaningful value would be uncapped. Give `kRep` units (pN), convert
   it, and cap the repulsion force.

### B. A robust base for other developers

10. *Deferred (special case, lower priority).*
    **The ParMETIS variant does not compile** *(verified)*.
    `helper/loadBalancer.cpp` uses block timers (`pf->timer`, `ff->timer`) that
    do not exist. CI cannot see it (no ParMETIS); build it against a stub
    header in CI, or drop the variant until it is finished.
11. *Done* (`resetGlobalState()` in `~HemoCell`; stretch and Lees-Edwards
    data per object). **Global state survives a simulation.** `~HemoCell`
    does not reset `param::*`, the `global` flags, or the static pointers and
    data of `HemoCellParticleField`, `CellInformationFunctionals`,
    `HemoCellStretch` and the Lees-Edwards helper. Reset them in one place and
    test two sequential simulations.
12. *Done* (no `using namespace` in headers, drivers declare their own; the
    constants keep their names in namespace `hemo`, without macros;
    `unique_ptr`; `const std::string &`). `DESCRIPTOR` and `CEPAC_DESCRIPTOR`
    remain macros by design: they are compile-time options of
    `constant_defaults.h` (an alias template passed as a template template
    argument is not guaranteed to be the same template as the one it names,
    CWG 1286), and the start-up check of B14 covers them. **Public API.** `std::auto_ptr` in
    `getFlagMatrixFromSTL()` (and 10 drivers) was removed in C++17 and only
    compiles with libstdc++; the constructor takes `char *`; public headers
    contain `using namespace std/plb/hemo` (`hemocell.h` makes `Array`
    ambiguous); generic macros (`DESCRIPTOR`, `OUTPUT_*`, `WBC_SPHERE`)
    pollute user code. Use `std::unique_ptr`, `const std::string &`, and
    namespaced constants with the old macros as aliases.
13. *Done* (`<sim><tmax>` optional for `sanityCheck()`; errors name file, line
    and element path). **Configuration friction.** `sanityCheck()` requires `<sim><tmax>` even in
    drivers that do not use it; "XML child X does not exist" names neither the
    parent path nor the file.
14. *Done* (reworded after review: keep `config/constant_defaults.h` as the
    one readable place). **Compile-time configuration.** `SOLIDIFY_MECHANICS`
    was always defined in `config/constant_defaults.h`, so
    `hemocell_solidify_mechanics` equalled the base library; a case compiled
    with other settings than the library it links (layout-changing macros,
    `FORCE_LIMIT`) was not detected; the configuration of a run was not
    recorded; the version was defined both in CMake and in the header. Now:
    the features default to off in the header and the variants switch them
    on; the case's configuration is compared with the library's at start-up;
    it is logged and stored in checkpoints; CMake and the docs read the
    version from the header.
15. *Done* (`getFlagMatrixFromSTL(cfg, ...)` in 14 drivers,
    `HemoCell::initializeLatticeFromFlagMatrix()` in the 5 pipe drivers and the
    tests). **Driver boilerplate.** Seventeen drivers copy the STL voxelisation
    and lattice setup; move it into a library helper.
16. *Done* (removed; `NoOp` fixed and tested). **Dead code.**
    `HemoCellFields::HemoSyncEnvelopes` (declares
    `modif::hemocell`; applying it would bring back arrival-order
    communication), `HemoCellFields::envelope_communicator`,
    `tests/test_example.cpp` (`Dummy.foobar`), `mechanics/NoOp.h` (does not
    compile when used: no matching `CellMechanics` constructor, and its
    `ParticleMechanics()` does not override; documented in
    `mechanical_models.rst` as the null model).

### C. Tests, CI and tooling

17. *Done* (`test_mechanics.cpp`: bending is not torque-free, documented;
    conservation across periodic and block boundaries, pre-inlet transfer,
    determinism; test output in a temporary directory).
    **No unit tests of the mechanics.** Per model (RBC, PLT, WBC, malaria):
    zero force at the rest shape, zero net force and torque of the internal
    forces, symmetry. Also missing: conservation of cells crossing periodic and
    block boundaries, the pre-inlet, and a determinism test (run twice, compare
    the checkpoints with `cmp`). Tests leave their output directories in the
    repository root (`tmp`, `tmp_0`, `tmp_1`, ...: one per simulation a test
    constructs); write them to a temporary directory instead.
18. *Done* (`scripts/ci/smoke_drivers.py`, run in the gcc release CI job for
    34 of the 39 built drivers: `cube` and `performance_testing` need
    generated configurations, `microcontraction`, `stent_strut` and
    `unbounded` take minutes (`--slow`). Not built at all, so not covered:
    `cases/flowchamber_stenosis`, `cases/injured_vessel` and
    `cases/stent_strut/stent_strut_crossings`, which are missing from
    `cases/CMakeLists.txt`.)
    **CI never runs a driver.** A 20-step smoke run of every example and case
    (with output and a checkpoint) would have caught A5 and B13; 31 of 36 pass
    today on 2 ranks, the five large ones (`AR2*`, `stent_strut`, `unbounded`)
    need reduced configurations.
19. *Done* (scripts ported, `py_compile` in CI; `packCells`: CMake >= 3.16,
    exact periodic box, periodic neighbour search fixed, documented).
    **Scripts and tools.** Eight of 42 Python scripts are Python 2 only
    (including `scripts/hemocell_h5utils.py` and
    `scripts/ParticleField3D_HDF5toXMF.py`); add `python3 -m py_compile` to CI.
    `tools/packCells` requires CMake 2.8 and packs periodically (cells cross
    the box faces), which is not documented (such files need
    `setPlaceCellsAcrossPeriodicBoundaries(true)` when packed for the domain).
20. *Done* (`nearestNode()`). **Rounding of negative coordinates.**
    `plint(x + 0.5)` truncates towards zero, mapping (-1, -0.5) to node 0, in
    `advanceParticles()` (wall removal), `nearestCell()`, the particle grid
    and the position loader. Low impact; use one `nearestNode()` helper.

### D. Model limits to document

21. The singular area, volume and bending limits (`MaxCell*` in
    `config/constant_defaults.h`) are shared by all cell types, so lowering
    `kArea` does not let a WBC stretch beyond 30 % local area change; reference
    lengths (7.5 nm persistence length, 0.5 um, 1280 faces) are hard-coded in
    `CellMechanics`; the WBC template's inner-link strengths are 1e-3 and 1e-6
    lattice force units at dx = 0.5 um, dt = 0.1 us (6.4 pN, 0.0064 pN).
22. RBC meshes drift at high shear: beyond a wall strain of about 2 at
    >= 2000 1/s, RBC triangles pass the 30 % local area limit and the fraction
    grows (`cases/pipeflow_wbc`; the same on `master` with the stock pipeflow
    example). Neither updating every step nor a larger `kArea` prevents it.
    Document in `developer/known_limitations.rst`.

### E. Found while testing C

23. **Pressure boundaries on a slab.** Nine drivers pass a box three nodes
    thick to `addPressureBoundary0N()`/`0P()`, whose precondition requires a
    plane: `examples/pipeflow_with_preinlet`,
    `examples/curvedflow_with_preinlet`, `cases/AR2`, `cases/AR2_pulsatile`,
    `cases/AR2_stiff`, `cases/stl_preinlet`, `cases/preinlet_shear`,
    `cases/fc_solidified_stenosis` and `cases/injured_vessel`
    (`bb.x0 .. bb.x0 + 2` or `bb.x1 - 2 .. bb.x1`). With `PLB_DEBUG` (the
    `debug` preset) they abort at set-up; otherwise the Zou/He pressure
    condition is applied on three layers. Use the face plane; this changes
    their results. (`tests/integration/test_preinlet_transfer.cpp` uses the
    plane.)
24. **Unregistered drivers.** `cases/flowchamber_stenosis`,
    `cases/injured_vessel` and `cases/stent_strut/stent_strut_crossings` are
    missing from `cases/CMakeLists.txt`, so they are never built, run or
    smoke-tested; build them or move them to `legacy/`.
