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

Status: A1-A9, B11-B16, C17-C20, E23-E24 and F26-F35 are done (see
CHANGELOG). B10 is deferred (special case, lower priority). D21-D22 and E25
are open, and the vRBC ellipsoid (F27).

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
    35 of the 42 drivers: `cube` and `performance_testing` need generated
    configurations; `microcontraction`, `stent_strut`,
    `stent_strut_crossings`, `flowchamber_stenosis` and `unbounded` take
    minutes (`--slow`). The last three cases were registered in E24.)
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

23. *Done* (one plane each; `curvedflow_with_preinlet` also had the wrong
    orientation, `0N` on its `+x` face, and lost up to 19 % of its flow).
    **Pressure boundaries on a slab.** Nine drivers pass a box three nodes
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
24. *Done* (all three registered, compiled after small fixes, smoke-tested;
    `flowchamber_stenosis` and `stent_strut_crossings` with `--slow`).
    **Unregistered drivers.** `cases/flowchamber_stenosis`,
    `cases/injured_vessel` and `cases/stent_strut/stent_strut_crossings` are
    missing from `cases/CMakeLists.txt`, so they are never built, run or
    smoke-tested; build them or move them to `legacy/`.
25. **The voxeliser opens only one plane at `x0`.**
    `getFlagMatrixFromSTL()` opens both x ends by copying the neighbouring
    plane (`CopyFromNeighbor`), but at the low end it loops in the direction
    of the copy: `x0` gets the old `x0 + 1` (the STL's end cap) before
    `x0 + 1` is opened, so for capped geometries `x0` stays mostly solid
    (`cases/AR2*`: 17 of 882 fluid nodes; `cases/injured_vessel`: 0), while
    the high end is opened fully. The outlets of those drivers therefore sit
    at `x0 + 1`. Fixing the loop order would change their geometry, and would
    open the closed low end of `examples/curvedflow_with_preinlet` (a U-bend
    with both openings at `x1`), which needs a way to keep an end closed.

### F. `tools/packCells`

Review of the packing tool (November 2026), revised with the author's
comments, which every change here has to respect:

* It must stay fast at high hematocrit (0.3 and above); today a 40 um cube at
  Ht 0.3 (198 RBCs) packs in under a second.
* Letting the ellipsoids shrink below their nominal size is a deliberate
  relaxation that makes high volume fractions reachable, at some cost in
  accuracy. Keep it, but make it explicit and bounded.
* The enlarged ellipsoids and the 97 um^3 RBC volume describe the cells'
  *effective* size, which includes the IBM interaction range (cells interact
  from about one kernel width apart). For collision-based viscosity this
  matters more than the mesh.

So scale 1 (the nominal ellipsoids) means "out of each other's interaction
range". Between 1 and the scale at which an ellipsoid still contains its mesh
(the type's *floor*: RBC 0.93, WBC 0.95), cells start within each other's
IBM range; below the floor their meshes may intersect.

Measured with an independent check of the output (the Perram-Wertheim
contact function over all periodic images, and the RBC shape itself) and the
meshes HemoCell builds from the templates. Confirmed correct: the angles
packCells writes and the loader's rotation agree to 1.6e-13 over 1e5
orientations, and the RBC and PLT meshes have their thin axis along y, as
packCells assumes. Suggested order: F26-F28 (small, no cost in run time),
then F29 with the check of F35, then F30-F31, then the rest.

26. *Done* (largest diameter for the contact range and the rotation step;
    the vRBC packing is now clear at its final scale, 1.53).
    *(verified)* **The contact range uses only the x diameter.** `Rcut` and
    `Rc_max` read `getSize()[0]`, so cells whose largest diameter is not
    along x miss contacts: vRBCs (3.5 x 6 x 11 um) packed with `--scale 0.25`
    reported no contact at scale 1.53, yet overlap at scale 1 (min 0.95). Use
    the largest diameter; RBC, PLT and WBC are unaffected (x is their
    largest), so it costs them nothing. `Epsilon_rot` is set from the first
    species only.
27. *Done* (PLT 2.7 x 1.17 x 2.7, PLT_m 2.05 x 1.17 x 2.05, PLT_mko 3.68,
    RBC_m 6.23 x 3.26 x 6.23 um: each mesh fits at 0.93; the sources are
    commented in `packCells.cpp`). Open: the vRBC ellipsoid (3.5 x 6 x 11 um)
    does not contain the gametocyte mesh either (`vRBC_uniform.stl` at radius
    5.4 um, 2.94 x 4.83 x 10.8 um, fit 1.09); with the RBC's margin it would
    be 3.85 x 5.3 x 14.5 um, 28 % more volume and much longer, which changes
    how densely they pack. Left as it is until that is decided.
    *(verified)* **Some ellipsoids are smaller than the cells.** Scale at
    which each ellipsoid contains the mesh HemoCell builds from its template:
    RBC 0.93, WBC 0.95, but PLT 1.04 (2.4 x 1.05 um against 2.50 x 1.09),
    PLT_m 1.03, RBC_m 1.00 (no margin) and PLT_mko 2.0 (1.71 um is the
    template's radius, used as the diameter). Give them the RBC's relative
    margin. Platelets are 0.14 % of the ellipsoid volume at the default PLT
    ratio, so this costs no packing density or time.
28. *Done* (inputs checked once before packing, exit status 1 with a
    message, 2 when not finished; `--seed`, reproducible with one thread;
    help rewritten; `-h` stays the hematocrit, as before).
    *(verified)* **Command line.**
    * Domain sizes are read with `atoi` (`30.5 20.7 20.9` becomes
      30 x 20 x 20).
    * Two crashes: `-s` is declared without an argument (`"h:s"`), and
      `--cell` as the last option reads past `argv` (check `optind + 4`).
    * Zero cells (`--rbc 0`, `--rbc abc`, a hematocrit that rounds to 0) loop
      forever on NaN.
    * Repeated names (`--rbc 3 --rbc 5`, or `--cell RBC ...` with `--rbc`)
      overwrite `RBC.pos`.
    * Hitting `--maxiter` writes the files and exits 0, and its `WARNING`
      (`NumSteps == Max_steps - 1`) never fires.

    Validate all inputs once at the start; exit non-zero when the packing is
    not finished. Add `--seed`: today it is `time(NULL)`, so no packing can be
    reproduced and two runs in the same second are identical. Fix the help:
    * `-h` is the hematocrit;
    * it names `--RBC`, `--PLT` and `--PLT_ratio`, which are rejected;
    * it promises `Cells.pov`, which `savePov()` (commented out) no longer
      writes;
    * `--cell <name> <n> <e1, e2, diameter>` takes the three full diameters
      along body x, y, z in um (y is the thin axis);
    * the "Volume" it prints is the ellipsoid's.
29. *Done* (`--start_scale`, default 1; `--min_scale`, default the largest
    mesh fit of the types packed, below which the packing stops with exit
    status 3; `--scale` removed, the grid set from the largest diameter; the
    final scale reported). With the defaults (40 um cube, 24 threads): Ht 0.3
    at 1.00 in about 120 iterations (< 0.1 s, also with PLTs and in a
    100 x 50 x 50 um box), Ht 0.4 at 1.00 (850, 0.4 s), Ht 0.45 at 0.96
    (7,500, 4 s); Ht 0.5 and 0.55 stop at the RBC's 0.93 after 14,000
    iterations (8-9 s; `--min_scale 0.8` packs Ht 0.55 at 0.90 in 13 s).
    200 RBC + 150 PLT end at 1.00 (were 0.81); vRBCs start at 1.14 (5 %
    above their mesh fit, 1.09), so dense vRBC packings now stop instead of
    returning intersecting cells. An independent overlap check (the script
    F35 would ship) confirmed every final scale.
    *(verified)* **Make the relaxation explicit and bounded.** The ellipsoids
    start at `scale * (mean dx*dy*dz)^(1/3)` and only shrink, so `--scale`,
    documented as a grid setting, sets the start and with it the final margin
    and the run time.
    * The automatic scale (1.2 / the largest diameter of the most numerous
      type) starts RBCs at 0.97, which is fast and above their floor. But a
      mixture's small cells pull it down: 200 RBC + 150 PLT end at 0.81, with
      186 RBC pairs whose shapes intersect; vRBCs start at 0.67.
    * A large scale is slow: `--scale 0.25` starts at 1.69 and takes 79k
      steps (41 s), since each step shrinks by a fixed fraction (`Ntau`).

    Size the neighbour grid from the largest contact range, and give the
    start scale directly (default 1.0). Measured with the start set that way
    (emulated with `--scale`; 40 um cube unless noted; 24 threads):
    * RBC, Ht 0.3, start 1.0 / 1.05: 123 / 230 steps, < 1 s, final scale
      1.00 / 1.05 (100 x 50 x 50 um, 773 RBCs: 146 steps);
    * RBC, Ht 0.4, start 1.0 / 1.05: 700 / 9,300 steps, 1 / 5 s, final
      1.00 / 1.00;
    * RBC, Ht 0.45, start 1.0: 7,600 steps, 5 s, final 0.96;
    * RBC, Ht 0.5, start 1.0: 14,400 steps, 10 s, final 0.93 (the floor);
    * 200 RBC + 150 PLT, start 1.0: 232 steps, < 1 s, final 1.00.

    Keep shrinking below 1 as the relaxation, down to the floor of the types
    present. Report the final scale against 1 and the floor; below the floor
    exit with an error (or require an explicit `--allow_overlap`). Drop the
    "most numerous species" heuristic and its warning.
30. *Done* (cell lists complete before a parallel search in which every
    thread writes only the contacts of its own particles; contacts applied
    and `Force_step` summed in a fixed order; `Rcut` local). Results are
    identical on 1, 4 and 24 threads, and byte-identical to the earlier code
    on one thread in all ten test cases when both are built without FMA
    contraction (with `-march=native` the compiler fuses multiply-adds
    differently, which changes long runs in the last bits). A bounding-sphere
    test skips the contact function for pairs that cannot overlap, and the
    search runs in chunks of 4 (chunks of 16 left most threads idle). On 24
    threads: Ht 0.45 in a 40 um cube 85 instead of 578 us per iteration
    (0.7 s instead of 4.4 s), in 100 x 50 x 50 um 352 instead of 2,768 us
    (2.9 s instead of 22.5 s); one thread is about as fast as 24 were.
    **Data races in the OpenMP loop** (from reading; no missed contact in 16
    repeated runs). `calc_forces()` builds the cell list incrementally inside
    the parallel loop over particles, so which pairs a thread sees depends on
    timing. `Rcut`, `Rcut2` and `Dinner` are shared members written by every
    thread, the partner's force is updated without atomics, and
    `calc_motion()` sums `Force_step` without a reduction. Results are not
    reproducible, even with a seed. Keep the speed-up (24 threads: 41 s,
    1 thread: 179 s, same packing):
    * build the cell list before the force loop;
    * accumulate pair forces in per-thread buffers and reduce them (the same
      pair count as now);
    * make `Rcut` local, and reduce `Force_step` and `Dinner`.

    The grid cell size may win time back (larger cells: fewer cells to visit,
    but more candidate pairs). Compare the same case and thread count before
    and after.
31. *Done* (uniformly random start orientations, Shoemake; spheres and
    `--noRotate` keep theirs). At Ht 0.1 no RBC is left unrotated and the
    thin axes are isotropic (mean |cos| 0.50 / 0.47 / 0.52 along x / y / z);
    run time and final scale unchanged (Ht 0.45, four seeds); `--noRotate`
    gives the same files as before.
    *(verified)* **All cells start with the same orientation** and rotate
    only when they collide. At Ht 0.1, 19 of 66 RBCs come out unrotated and
    the RBCs lie mostly in the xz-plane (mean |cos| of the thin axis with y:
    0.83; isotropic is 0.5). Start from uniform random orientations; tried at
    Ht 0.3 and 0.4, the step count is unchanged and the result isotropic
    (0.45-0.53). Keep `--noRotate` (aligned and fixed) as an option.
    `Quaternion::random()` is unused and not uniform.
32. *Done* (`--walls <axes>`, `--cylinder <axis> <R>`, `--offset`; all
    periodic images within reach in thin boxes, a cell's own included; the
    user guide gives box, radius and offset for walls at given nodes). An
    independent check (pairs through the periodic axes only, walls from 4000
    surface points per cell) confirms the final scale of slabs, channels,
    closed boxes and pipes; the cylinder contact matches brute force to
    1e-5. In a 6 um slab the RBCs turn flat (mean |cos| of the thin axis
    with z 0.86). HemoCell placed every cell of `parallelplanes` (`--walls z`)
    and `pipeflow` (`--cylinder`) packed this way; periodic packings of the
    whole domain lost 63 of 773 RBCs and 62 of 89 RBCs at the walls. The safe
    pipe radius, (r - (k + 0.5) sqrt(2)) dx, was measured on the voxelised
    tube (22.92 of 25 nodes for k = 1). With walls, Ht 0.45 in a 100 x 50 x
    50 um channel or pipe packs at 0.95 in 2.2-2.5 s (about 200 us per
    iteration, as without walls). Both drivers skip the cells crossing their
    periodic x (and y) boundary; they could call
    `setPlaceCellsAcrossPeriodicBoundaries(true)` for packings of their
    domain.
    *(verified estimate)* **Walls.** The packing is periodic along all three
    axes, and the loader drops every cell crossing a wall. In the Case.rst
    tutorial (25 um cube, walls in y and z), 29 of the 48 RBCs touch the
    walls, so the requested Ht 0.3 becomes at most 0.11 (a geometric
    estimate, not run in HemoCell). A box thinner than a cell is never checked
    against the cell's own periodic image. Add per-axis walls with a clearance
    (an ellipsoid-plane distance check, O(N) per step) and an origin offset,
    so the box can be the fluid region; optionally a cylinder for pipes
    (`cases/pipeflow_wbc/README.md`).
33. *Done* (the effective volume explained in the help and the user guide;
    per type the number of cells, the volume fraction by the mesh volume and
    for RBCs the hematocrit by the effective volume are printed).
    **Hematocrit and volume fractions.** `--hematocrit` counts RBCs by their
    effective volume, 97 um^3 (the 81 um^3 mesh plus the IBM interaction
    range). This is intended, but it is explained nowhere: the code only says
    "set to match the model". Document it in the help and the user guide,
    and print the cell counts with the effective and mesh volume fractions.
    `--hematocrit` always adds human RBCs (no `RBC_m`; WBCs are not counted),
    and the "nominal volume fraction" printed is the ellipsoids'.
34. *Done* (CMake: build type, `find_package(OpenMP)`, `PACKCELLS_NATIVE_ARCH`,
    `PACKCELLS_WERROR`, no `g++-7`; the binary stays in `tools/packCells/`
    unless `-DCMAKE_RUNTIME_OUTPUT_DIRECTORY` is given; `register` and the
    dead code removed; centres in [0, L), written with 9 digits). Found on
    the way: build trees sharing the output path overwrote each other's
    binary; FMA contraction made every build (OpenMP or not, `-march=native`
    or not) pack differently for a seed, now off (`-ffp-contract=off`, no
    cost in speed; the gcc variants give identical files); a packing whose
    overlaps balanced (a cell pressed into both walls of a slab) was
    reported done, now it also requires that nothing overlaps. Clang is not
    installed here: the CI job builds packCells with clang.
    **Build and code.**
    * The CMake file forces `g++-7` on macOS (after `project()`) and
      hard-codes `-fopenmp -march=native -O3 -g`; use `find_package(OpenMP)`
      and the build type.
    * The binary always lands in `tools/packCells/`.
    * `register` (removed in C++17), and a non-inline `operator*` in
      `geometry.h`.
    * Dead code: `savePov()`, `getOutput()`, `scene.pov`, `compile.sh`.
    * A centre exactly at L stays there (`pbc()` uses `>`), which the loader
      rejects.
35. *Done* (user guide page `packcells.rst` under How-to Guides, the README,
    the QuickStart; `check_packing.py`; the CI job `packcells`, with gcc and
    clang, packs and checks a periodic box, a slab and the pipeflow pipe in
    about 2 s).
    **Documentation and a check.**
    * QuickStart.rst gives the path `./hemocell/packCells`, and its build
      steps don't put the binary in `build/`.
    * Case.rst shows an old help and output, and its command (`... -r`) now
      fails ("option '-r' is ambiguous").
    * `tools/packCells/README` has only credits. Document:
      * the box: the periodic length is nx*dx um from lattice node 0, and
        positions are relative to the pre-inlet's corner;
      * that cells crossing walls are dropped, so the hematocrit must be
        checked after loading;
      * that the file names must match the `addCellType` names (`RBC_HO`,
        `WBC_HO` and `PLT_HO` in some cases);
      * the rotation convention (body thin axis y; degrees about x, then y,
        then z, about the centre);
      * each type's ellipsoid and what it means (effective size, floor);
      * the progress columns and what to check at the end;
      * the run time against the start scale.
    * Ship the overlap check (a short numpy script, Perram-Wertheim) with the
      tool, and run it in CI on a small RBC + PLT packing at Ht 0.3: the final
      scale at least the floor, no pair below it, and a time limit.
