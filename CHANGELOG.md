HemoCell Changelog
==================

In this changelog, we consider the following classes of changes:
* Features: any new features that are added to HemoCell.
* Structure: any changes to the structure of HemoCell that may break existing cases.
* Fixes: (small) changes that do not fall in the other categories.

Unreleased (refactor)
---------------------
* Fixes that change simulation results
  * Repulsion: particle pairs within the same lattice node received twice the
    force of pairs in neighbouring nodes; every pair is now counted once.
  * Membrane viscosity cap: the per-edge cap compared a lattice-unit force to
    `FORCE_LIMIT` in pN and therefore never acted; it now uses
    `param::f_limit / 4` (12.5 pN with the default limit). Only material files
    with non-zero `eta_m` (the WBC templates) are affected.
  * `setInitialMinimumDistanceFromSolid()` stored the distance in an integer, so
    fractional micrometres (e.g. the 0.5 used by pipeflow) were truncated to 0.
    The distance is now honoured, which removes more cells near walls at
    initialisation.
  * Writing output no longer alters the simulation: the fluid `OUTPUT_FORCE`
    reset the lattice force (removing a body force set before `writeOutput()`),
    and output recomputed repulsion and constitutive forces. Results are now
    independent of the output interval.
  * Output and statistics no longer change the trajectory at all.
    `writeOutput()` and the `CellInformationFunctionals` re-synchronised the
    particle envelopes, which reordered the particles (and with them the
    summation order of force spreading) and deleted incomplete cells: runs with
    and without output differed at round-off level, growing over time. They
    now read complete cells from a copy of the particles
    (`HemoCellFields::CompleteCellsScope`) that is discarded afterwards; a
    pipeflow run with output and statistics every 37 steps ends bit for bit in
    the same state as one without (stored residence times aside, which output
    brings up to date). Output is also faster (one synchronisation instead of
    several). Cells that lose a vertex (to a wall or open boundary, or through
    a non-finite force) are now deleted by `iterate()` at the next particle
    update; before, they were deleted only at output or statistics steps (or
    with `<cellsDeletedInfo>`), and in between their remaining vertices kept
    moving and, with a material time-scale separation, spreading their last
    forces. Runs in which no cell breaks are unchanged without output.
* Fixes
  * Particle neighbour grid: replaced the fixed 10-particles-per-node grid
    (unchecked heap overflow) by a compact cell list without a capacity limit.
  * `syncEnvelopes()` uses a private MPI communicator and receives requests per
    known source rank, so messages of consecutive calls cannot be mismatched.
  * Default periodicity limits are `{100,100,100}` (the y and z limits were 0).
  * Checkpoints additionally store the number of cells, the last output
    iteration (residence time) and the CEPAC field; restarts no longer depend on
    unchanged `.pos` files. Older checkpoints still load.
  * Load balancing (ParMETIS build): keep `domain_lattice` valid after the
    lattice is replaced, rebuild the particle communication structure, read
    the checkpoint from the directory it is written to, and refuse
    combinations with a pre-inlet or the CEPAC field. `doRestructure()` warns
    instead of silently doing nothing without ParMETIS.
  * Configuration: invalid, empty or partially numeric values (e.g. `1e6` for
    an integer) raise `std::invalid_argument` instead of silently becoming 0
    (or 1); malformed XML is reported; XML without a `<?xml?>` prolog is read
    correctly; booleans accept `true`/`false`; strings keep inner spaces.
  * `FORCE_LIMIT` can be undefined again without breaking compilation.
  * A warning is printed when a repulsion cutoff exceeds one lattice spacing
    (the neighbour search only covers adjacent nodes).
  * Exit signals (SIGINT, SIGTERM, SIGHUP, SIGQUIT, SIGUSR1/2) are handled
    collectively every `exitSignalCheckInterval` (100) iterations: all ranks
    save a checkpoint and exit, instead of one rank exiting and the others
    hanging. SIGABRT is no longer intercepted.
  * Log messages written before the log file is opened (e.g. all setup output
    of a restarted run) are printed and kept, instead of being dropped.
  * Pre-inlet: velocities are exchanged with one message per rank pair per
    step from a plan built once (previously one blocking message per inlet
    node, tagged with the node index, which can exceed MPI's tag range), the
    inlet nodes are selected by their fluid flag (previously `isBoundary()` was
    evaluated at unrelated nodes), particles are sent as one message per rank
    pair (several blocks per rank were mismatched), received blocks no longer
    shift a shared buffer twice, and the three global barriers per step are
    gone. Results of the examples are bit-identical.
  * `hemo::Array` no longer uses libstdc++ internals (portability) or writes
    past its end for long initialiser lists; `+= Dot3D` adds the offset once.
  * Particles have move semantics (no deep copies on reallocation/removal);
    copy assignment copies `force_total`; `HemoCellParticleField::operator=`
    copies particles instead of deleting itself (Palabos uses it when copying a
    MultiParticleField3D). The Adams-Bashforth integration option compiles.
  * Debug (`PLB_DEBUG`) builds compile again; non-standard VLAs removed.
  * `deleteIncompleteCells(ctype)` only deletes cells of that type;
    `getNumberOfCells_Global()` returns the actual count (was 0).
  * The binding-field and interior-viscosity helpers are owned by
    `HemoCellFields` instead of being process-wide singletons that outlived it.
* Build system and CI
  * Palabos is compiled once (object library) instead of once per HemoCell
    library variant; the patched `parallelBlockCommunicator3D.cpp`, which calls
    HemoCell code, is compiled with the HemoCell sources.
  * The build type is honoured (default Release); new options
    `HEMOCELL_NATIVE_ARCH`, `HEMOCELL_DEBUG_INFO`, `HEMOCELL_PALABOS_DEBUG`,
    `HEMOCELL_SANITIZERS`, `HEMOCELL_WERROR`; `CMakePresets.json` with release,
    portable, debug and sanitize configurations.
  * MPI and HDF5 via imported/interface targets, with MPI's include directories
    before HDF5's (an HDF5 prefix may contain another MPI's `mpi.h`); the
    ineffective switch to the MPI compiler wrapper was removed.
  * Warnings only for HemoCell code; Palabos and tinyxml2 are compiled quietly.
  * CTest: `unit`, `validation`, `unit_mpi`, `validation_mpi` run in the
    repository root (the old entry ran in `tests/` and could not find its
    inputs); `hemocell_test` is written to the build directory.
  * Aggregate targets `hemocell_examples` and `hemocell_cases`.
  * GoogleTest pinned by commit, an installed GoogleTest >= 1.17 is used when
    found; requires CMake >= 3.16 (presets: 3.21).
  * GitHub Actions: GCC (`-Werror`) and Clang builds of everything, unit and
    validation tests on 1 and 2 ranks, unit tests under ASan/UBSan, and the
    user guide. Performance tests moved to a manual self-hosted workflow.
  * Debug builds work: fixed a const-correctness error that only appeared with
    Palabos assertions, and platelet inner-link lengths no longer use Palabos'
    edge-length function (which asserts that the vertices share an edge).
  * UBSan findings fixed: unaligned record access in the gathering functional,
    indexing an empty receive buffer.
* Formatting and documentation
  * `.clang-format` (derived from the prevailing style), `.editorconfig`,
    `scripts/format.sh`; CI checks the formatting of the lines a pull request
    changes.
  * All HemoCell sources reformatted (Palabos, `external/` and `legacy/`
    untouched); the library, variants, examples and cases compile to
    identical code. The commit is listed in `.git-blame-ignore-revs`.
  * New Developer Guide in the user guide: architecture and internals, coding
    conventions (including performance rules), extending HemoCell, known
    limitations. `CONTRIBUTING.md`; contribution workflow targets `dev`.
  * Units, preconditions and collectiveness documented for the public
    `HemoCell` setters; the repulsion page states the actual units.
* Fixes found while documenting
  * Wall repulsion without cell-cell repulsion accumulated in
    `force_repulsion` at every update (it was never reset); it now resets the
    field itself. Cell-cell and wall repulsion must use the same update
    interval (checked at start-up).
* Performance (results bit-identical; cube benchmark on 4 ranks: `iterate`
  -5 %, `spreadParticleForce` -22 %, `advanceParticles` -35 %)
  * IBM kernel (`phi2`): only the 2x2x2 nodes with non-zero weight are visited
    instead of 27 candidates, with the same products and summation order.
  * Force spreading skips particle copies whose kernel lies outside the bulk
    and the collided envelope layer; the force sum is formed once per particle.
  * Velocity interpolation skips particle copies (the following
    `syncEnvelopes()` replaces their velocity anyway).
  * `advanceParticles()` only runs the removal pass when a particle entered a
    boundary node.
  * Cell-information CSV output gathers to rank 0 instead of to all ranks.
  * `.pos` files are read by rank 0 and broadcast, and parsed once per rank
    instead of once per atomic block (with several XML reads per cell). A
    missing file is reported once, not by every rank for every block; a file
    with fewer complete rows than its declared count uses those rows with a
    warning (previously the missing cells were placed at the origin).
* Reproducibility
  * Multi-rank runs are bit-for-bit reproducible: both particle exchanges of
    `syncEnvelopes()` unpack the received particles in ascending source rank
    instead of MPI arrival order. The arrival order decided the order of
    particles in a block, and with it the summation order of force spreading
    and repulsion, so identical runs differed at round-off level (relative
    1e-8 in forces after 100 steps of pipeflow on 4 ranks), growing over time.
    The envelope exchange, previously Palabos'
    `duplicateOverlaps(modif::hemocell)`, is done by HemoCell
    (`exchangeEnvelopes()`) over the same communication structure. No
    measurable cost.
  * Checkpoints of identical runs are byte-identical: the padding bytes of the
    particle records are zeroed (they held uninitialised memory).
* Output (simulation results unchanged; what the files contain changes).
  HDF5 output of a 4-rank pipeflow step: 6.1 MB -> 4.2 MB.
  * Every cell is written once, by the block whose bulk holds its vertex 0.
    Previously each block holding part of a cell wrote all of it, including
    both images of a cell crossing a periodic boundary (55 copies of 30 RBCs
    in that run), so analyses of the HDF5 files counted cells several times.
  * Fluid and CEPAC files hold the bulk nodes of their block only (previously
    one extra node layer around each block); `relativePosition` is the lower
    corner of the first node's cell (first node - dx/2). Each node is one cell
    of a regular grid, so the blocks tile the domain without overlap.
  * Datasets are byte-shuffled and deflated at level 4 (was level 7 without
    shuffle): about 10 % smaller and half the compression time.
  * SI output: `OUTPUT_DENSITY` is in kg/m^3 (it was multiplied by `df/dx^2`,
    giving three times a pressure); `OUTPUT_OMEGA` is dimensionless (was
    multiplied by `df/dx^2`); the CEPAC density (a concentration) is written
    unscaled.
  * `OUTPUT_CELL_DENSITY` used the block-local node index as array index, which
    placed the counts one node off per direction (with `fluidEnvelope` 2); the
    dataset is named `CellDensity_<name>` with the cell type's output name.
  * Cell types registered with a directory in their name (e.g.
    `addCellType("cells/RBC")`) wrote no HDF5 or CSV output (the files were
    created in non-existent subdirectories); output now uses the name without
    the directory (`HemoCellField::getIdentifier()`).
  * HemoCell writes XDMF 3 files at every output: `<name>.<iteration>.xmf` per
    step and `<name>.xmf` with all steps and their time. ParaView's XDMF3 and
    legacy readers load them (tested with 5.11 and 6.0). Previously
    `batchPostProcess.sh` generated them afterwards, and only the legacy
    reader could open them.
  * `scripts/batchPostProcess.sh` runs the new `scripts/HDF5toXDMF.py`
    (Python 3, h5py), which writes the same files for output of older
    versions. Removed: `FluidHDF5.py`, `CellHDF5toXMF.py`,
    `FluidHDF5toXMF.py`, `batchPostProcess_py2.sh` and the `*_py2.py`
    converters, `CellInfoMergeCSV.sh` (the CSV output is already one file per
    cell type and step). `convert_xmf_to_x3d.py` converts the per-step cell
    files only.
  * HDF5 property lists are closed (one leaked per dataset written).
* Hardening (results bit-identical)
  * Fatal errors go through `hemo::fatal()` (`config/error.h`) instead of 42
    `exit(1)` calls: the message is printed with the MPI rank on whichever rank
    detects the error (previously mostly through `hlog`, which prints on rank 0
    only, so errors found on other ranks were silent) and the job ends with
    `MPI_Abort`. `hemo::warning()` reports from any rank; `FatalErrorsThrow`
    makes `fatal()` throw, for tests of error paths.
  * A vertex with a non-finite force (e.g. from a degenerate triangle) is no
    longer spread into the fluid, where it made the whole lattice NaN: it is
    removed with a warning, which deletes its cell. The `FORCE_LIMIT` check
    catches it at no extra cost. Particle and wall repulsion skip coincident
    points instead of dividing by a zero distance.
  * Checkpoints: `checkpoint.xml` is written last, under a temporary name, so
    a save that is interrupted leaves no `checkpoint.xml` next to partial data
    files (the previous checkpoint is complete in the `.old` files). It records
    the particle record layout, the cell types (name and vertex count, in
    `addCellType()` order), `dx` and `dt`; a restart with a different library
    variant, cell types or resolution is a fatal error instead of misreading the
    raw particle data. Older checkpoints load with a warning.
    `renameFileToDotOld()` reports failed renames (its status was shadowed).
  * Initial positions: a missing `<name>.pos` file is a fatal error (it was a
    warning, and the case ran without those cells; `examples/microcontraction`,
    `cases/stentflow` and `cases/stent_strut` (whose files are named `.POS`)
    did so). `loadParticles()` logs how many cells of each type were placed out
    of those in the file, with a warning when cells were dropped.
  * Cell IDs are 64-bit (`plint`) throughout. An ID grows by a multiple of the
    number of cells at every periodic crossing and is never folded back, but
    the lookup maps, the statistics and the cell ID requests of
    `syncEnvelopes()` (sent as `int`) used `int`: with 200,000 cells an ID
    overflowed after about 640 laps in x, the periodic transfer then reset it
    to the base ID, and the cell was most likely deleted. This also lifts the
    limit of about 212,000 cells with the default periodicity limits;
    `sanityCheck()` requires room for at least 10^6 crossings per cell.
  * Structure: `CellMechanics::ParticleMechanics()` and `solidifyMechanics()`
    take maps keyed by `plint`, as do `get_particles_per_cell()`, `get_lpc()`,
    `CellInformationFunctionals::info_per_cell` and
    `calculateCellInformation()`; `getOutputCells()` returns `plint` IDs.
    Custom cell models need the new signatures; a `solidifyMechanics()`
    override with the old signature (not marked `override`) would silently no
    longer be called. `HemoCellGatheringFunctional::gather()` and
    `gatherToRoot()` accept maps with any integer key.
* Tests
  * Unit tests for the particle grid/repulsion, the configuration reader,
    `hemo::Array`, particle copy/move and log buffering; an integration test
    for checkpoint/restart.
  * Integration test for the output: every cell written once, every lattice
    node once, consistent cell density and XDMF files (1 and 2 ranks).
  * Tests for `fatal()` and fatal configuration errors, for a non-finite
    vertex force (fluid stays finite, the cell is removed), and for the
    checkpoint write order and fingerprint checks.
  * Tests that output and statistics leave the trajectory unchanged, that cell
    IDs beyond 2^32 give the same trajectory, that `iterate()` deletes a broken
    cell, and of `HemoCellGatheringFunctional` with 64-bit keys (1 and 2
    ranks).

2.6 (July 15 2022)
-----------------
* Features
  * Added functionality to hemocell.initializeLattice() to define domain decomposition. It can be read from the config xml (xml tags: mABx, mABy, mABz).
  * Added weak-scaling test config generator for the case ``cube``. The config template contains an example for the previous domain decomposition tags.
  * The cell packing tool has been improved significantly with better automatic resolution to avoid cell overlapping.
  * Added mouse blood cell types to the cell packing tool.
  * The cell templates (along with the new mouse cells) are available under ``tools/cell_templates``.
  * New module shell script for Cirrus (EPCC).
* Fixes
  * Various comment fixes in the code.
  * Better error messages.
  * Improved sanity-check function.
  * ``tools/pos_to_vtk`` can now be used without installation.
  * Changes to the cellPacking default values and messages to make it easier to use.

2.5 (Jan 27 2022)
-----------------

* Features
  * Added ``tools/pos_to_vtk`` for simple cell visualisation.
  * Included section `Common mistakes` in documentation.
* Structure
  * Changed default `make` target to HemoCell library only.
  * Replaces the macro `#HEMOCELL_PARTICLE_FIELD` with the class directly.
  * Moved to Palabos library version `v2.3.0`.
  * Reduced Palabos patch size.
* Fixes
  * Removed packing Palabos source, now pulled from their remote by `setup.sh`.
  * Updated preinlet communication to use a local buffer as temporary storage to
    avoid velocity artefacts due to reuse of the same velocity buffer for
    communication.
  * Updated preinlet communication structure to reduce the number of messages
    send by only considering overlapping and neighbouring atomic blocks as
    targets to send messages.
  * Preliminary fix for solidification mechanics on multi-core machines.

2.4 (Sep 23 2021)
-----------------
* Features
  * Added two capillary flow examples.
  * Added, updated, and extended documentation of most example cases.
  * Added helper functions in `helper/geometry.h` for boolean operations on domain functionals.
  * Added script to convert HDF5 to the X3D format to enable rendering with Blender.
  * Added validation tests based on cell stretching and pipe flow.
* Structure
  * Cell mechanics constants are moved to `config/constant_defaults.h`.
  * The Doxygen pages are automatically generated when building the Sphinx documentation.
  * The continuous integration pipeline is split into four stages: build, test, validation, and documentation.
* Fixes
  * Added guard statement to prevent reading uninitialised Lees-Edwards displacement values.

2.3 (Jun 1 2021)
----------------
* Features
  * Rebuilt the build environment. For more details, please see the quickstart guide in the documentation.
  * Added the cube example in which a cubic domain is subjected to a shear flow on the top/bottom planes.
* Structure
  * HemoCell now depends on Palabos 2.2.1 as the older version used previously is no longer available. The patch and setup.sh file has been updated to reflect this.
* Fixes
  * Remove some useless output in the cellStretch example

2.2 (Dec 7 2020)
----------------
* Features
  * Added support for Lees-Edwards boundary conditions (single core only for now, example soon to come with publication).
  * Significantly improved and updated the documentation.
  * Publications using HemoCell added to documentation.
  * Updated the website with more examples in the gallery.
  * Added further output fields: e.g. shear rate, strain rate, shear stress.
  * Several small updates and QoL improvements in the otput processing scripts.
* Fixes
  * Fix a bug in the velocity output when SI units were used

2.1 (Feb 6 2020)
----------------
* Features
  * Added support for binding sites
  * Added support for interior viscosity
  * Added velocity output for hemoCellParticleField
  * Added curved pipeflow with pre-inlet example
  * Added microcontraction and flowaroundsphere examples
  * Added support for pulsatility
* Structure
  * Rename HemoCell::setMinimumDistanceFromSolid to HemoCell::setInitialMinimumDistanceFromSolid
  * Rename all occurences of RBC_HO to RBC, both in filenames and source code
  * Merge csv and hdf5 reader for scripts
  * The preinlet setup now adjusts the preinlet to always be on the lattice, if the lattice is set up before the preinlet (a warning is given otherwise)
* Fixes
  * Do not override interrupt handlers until the first iteration has started
  * Added missing AGPL license header to several files
  * Added profiler statements to preInlet.cpp
  * Log an error message whenever a nonexisting cell type is requested
  * Fix checkpointing so that it works with the pre-inlet
  * Update scripts so that they are compatible with the new RBC output
  * Fix outputStrainRate in IO/FluidHdf5IO.hh
  * Fix bug in the pre-inlet when a block boundary is spanned
  * Fix computation of the Tresca value
  * Fix location of Palabos download page
  * Fix for the environment script on Cartesius
  * Set the driving force at the warmup part for the pipeflow example

2.0b (Feb 12 2019)
------------------
to be added

2.0 (Dec 17 2018)
-----------------
* Added compile time constants PREINLET_MECHANICS and SOLIDIFY_MECHANICS to constant_defaults
* Added runtime preinlets, see stl_preinlet for examples, instead of preinletFromSlice() you can call autoPreinletFromBoundary() if the boundary happens to be the inlet as well
* Interior viscosity edge tracking works again, enable interior viscosity through CELL.xml file, see cases/RBC_template.xml for an example
* Output is saved in a seperate folder for each run, see pipeflow/config.xml for the options
* CSV output is already concatenated in csv folder, now also includes base_cell_id, correct cell_id and velocity

1.4 (Jul 9 2018)
----------------
* Added hemocell.statistics as an advanced (fully recursive) profiling tool
* Added OUTPUT_INNER_LINKS to output possibility

0.5
---
* Inner links are specified in the `<CellType>.xml` file through "CellMechanics"->"InnerEdges"->"Edge"
* Added MESH_FROM_STL for addCellType. Takes a triangulized STL file from xml-file->"MaterialModel"->"StlFile".
* Added "Volume" option for Celltype xml files
* Added Celdensity output for fluid field
* Added preInlet possibility and cases
* Default is now HDF5 with openmpi support, if not present, using a preInlet will generate an error when using multiple processors
* variables checkpointDirectory and outputdirectory can now be specified in the config.xml

0.4
---
* Periodicity limit can now be set (hemocell::setSystemPeriodicityLimit) per axis, default is 100, this means that a cell can wrap around a maximum of 100 times over a given axis
* Loadbalancer now reloads the old (smallest possible) atomic blocks before load balancing, this means timing information is not available as weight, only the number of lsps (or any other measure you come up with)
* The CmakeFile for a case can be generated from cases/CMakeLists_template.txt by the makefile (cases/Makefile). running make in ./cases will create and build all cases.

