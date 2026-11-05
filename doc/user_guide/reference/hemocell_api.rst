The HemoCell API
================

This page is a task-oriented guide to the ``hemo::HemoCell`` class — the object
you drive from a case's ``main()``. It groups the most commonly used methods by
what you do with them. For full signatures, overloads, and the rest of the C++
API, see the auto-generated :ref:`Doxygen reference
<hemocell_doxygen:HemoCell Doxygen Documentation>`.

A minimal case follows the same skeleton in every example (see
:ref:`Case:Creating your own HemoCell Case` for the full walk-through). The
library's headers import no namespace; the examples start with
``using namespace hemo; using namespace plb; using namespace std;``.

.. code-block:: c++

   HemoCell hemocell(argv[1], argc, argv);   // 1. construct
   hemocell.initializeLattice(management);   // 2. set up the fluid
   hemocell.initializeCellfield();           // 3. set up the cells
   hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
   hemocell.setOutputs("RBC", { /* ... */ });
   hemocell.setFluidOutputs({ /* ... */ });
   hemocell.loadParticles();                 // 4. place the cells
   while (hemocell.iter < tmax) {            // 5. run
       hemocell.iterate();
       if (hemocell.iter % tmeas == 0) hemocell.writeOutput();
   }

Construction and lifecycle
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Member
     - Purpose
   * - ``HemoCell(configFile, argc, argv)``
     - Construct the simulation from a ``config.xml``. ``argc``/``argv`` are
       passed through to Palabos for MPI initialisation.
   * - ``HemoCell(configFile, argc, argv, MPIHandle::External)``
     - Same, but for when MPI is initialised outside HemoCell (HemoCell then
       skips ``plbInit``).
   * - ``~HemoCell()``
     - Ends the simulation and resets the process-wide state (``param::``,
       ``global``, the log), so another ``HemoCell`` can be created afterwards.
       Only one may exist at a time.
   * - ``iter``
     - The current iteration counter (read it to drive your main loop).
   * - ``cfg``
     - The parsed configuration; read values with
       ``(*hemocell.cfg)["domain"]["dx"].read<T>()``.

Setting up the fluid
--------------------

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Member
     - Purpose
   * - ``initializeLattice(management)``
     - Create the fluid lattice from a Palabos ``MultiBlockManagement3D``. Call
       after defining any pre-inlets and before ``initializeCellfield()``.
   * - ``initializeLatticeFromFlagMatrix(management, flagMatrix)``
     - Create the fluid lattice of a voxelised geometry (no pre-inlet):
       bounce-back walls where the flag matrix is 0, no periodicity, fluid at
       rest. Then call ``lattice->initialize()``. The geometry itself comes from
       ``getFlagMatrixFromSTL(*cfg, voxelizedDomain, flagMatrix)``, which reads
       ``<domain><geometry>`` and the related settings (see ``examples/pipeflow``).
   * - ``lattice``
     - The underlying Palabos ``MultiBlockLattice3D`` (set periodicity,
       boundaries, driving force, ... directly on it).
   * - ``latticeEquilibrium(rho, vel)``
     - Initialise every fluid node to a density and velocity (lattice units).
   * - ``setSystemPeriodicity(axis, bool)``
     - Make the whole system periodic along an axis.
   * - ``setSystemPeriodicityLimit(axis, limit)``
     - Number of periodic images of a cell distinguishable along an axis
       (default 100); see :ref:`developer/architecture:Cell IDs and
       periodicity`.

Setting up the cells
--------------------

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Member
     - Purpose
   * - ``initializeCellfield()``
     - Create the cell-fields / particle-field container. Call before adding
       cell types.
   * - ``addCellType<Model>(name, constructType)``
     - Register a cell type with a :ref:`mechanical model
       <reference/mechanical_models:Mechanical (cell) models>` and a
       :ref:`construction method <concepts/cell_model:Building the mesh>`. The
       ``name`` must match the ``name.xml`` / ``name.pos`` files.
   * - ``loadParticles()``
     - Load the initial cell positions from the ``*.pos`` files. Call after all
       cell types, outputs, and periodicity are configured.
   * - ``setInitialMinimumDistanceFromSolid(name, d)``
     - Keep cells at least ``d`` (in **micrometres**) away from solids at
       initialisation: vertices within ``floor(d / dx)`` lattice nodes of a
       solid node are rejected, and cells with rejected vertices are removed.
       Must be called **before** ``loadParticles()``.
   * - ``setPlaceCellsAcrossPeriodicBoundaries(bool)``
     - Also place the cells that cross a periodic boundary (as periodic images
       on both sides); only for ``.pos`` files packed periodically for this
       domain. Off by default: such cells are skipped and counted in the log.
       Call after ``setSystemPeriodicity()``, before ``loadParticles()``.

Controlling cost (time-scale separation)
----------------------------------------

See :ref:`concepts/units_and_scaling:Time-scale separation` for the concept.

.. list-table::
   :header-rows: 1
   :widths: 52 48

   * - Member
     - Purpose
   * - ``setMaterialTimeScaleSeparation(name, n)``
     - Re-evaluate this cell type's mechanical model every ``n`` fluid steps.
   * - ``setParticleVelocityUpdateTimeScaleSeparation(n)``
     - Interpolate fluid velocity to particles every ``n`` fluid steps.
   * - ``setRepulsionTimeScaleSeperation(n)``
     - Update cell–cell repulsion every ``n`` fluid steps.
   * - ``setInteriorViscosityTimeScaleSeperation(n, n_grid)``
     - Update interior viscosity every ``n`` steps, with a full-grid pass every
       ``n_grid`` steps.

Optional physics
----------------

.. list-table::
   :header-rows: 1
   :widths: 48 52

   * - Member
     - Purpose
   * - ``setRepulsion(constant, cutoff)``
     - Enable and parameterise short-range cell–cell repulsion (constant in
       pN, cutoff in micrometres; the force is capped at ``FORCE_LIMIT``). See
       :ref:`advanced_cases/repulsion:Repulsive forces`.
   * - ``enableBoundaryParticles(constant, cutoff, timestep)``
     - Enable repulsion between cells and solid boundaries (same units;
       ``timestep`` in iterations, equal to the cell–cell repulsion interval if
       both are enabled).
   * - ``enableSolidifyMechanics(name)``
     - Enable solidification (clotting) for a cell type. Requires the
       ``SOLIDIFY_MECHANICS`` build (see :ref:`constants`).

Output and checkpointing
------------------------

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Member
     - Purpose
   * - ``setOutputs(name, {...})``
     - Choose which per-cell fields to write for a cell type. See
       :ref:`reference/outputs:Output fields`.
   * - ``setFluidOutputs({...})``
     - Choose which fluid fields to write.
   * - ``setCEPACOutputs({...})``
     - Choose outputs for the advection–diffusion (CEPAC) field, when used.
   * - ``outputInSiUnits``
     - ``true`` (default) writes output in SI units; ``false`` writes lattice
       units.
   * - ``writeOutput()``
     - Write the requested HDF5/CSV output for the current iteration. It does
       not change the simulation state (a body force set before it is kept).
   * - ``saveCheckPoint()`` / ``loadCheckPoint()``
     - Write / restore a full checkpoint (see :ref:`QuickStart:Resuming from a
       checkpoint`).

Running and load balancing
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Member
     - Purpose
   * - ``iterate()``
     - Advance the coupled fluid + cells by one fluid time step. If you drive
       the flow with an external force vector, re-apply it after every
       ``iterate()``.
   * - ``checkExitSignals()``
     - Collective. If any rank caught SIGINT, SIGTERM, SIGHUP, SIGQUIT,
       SIGUSR1 or SIGUSR2, all ranks save a checkpoint (unless
       ``checkpointOnExitSignal`` is ``false``) and exit with code 1.
       ``iterate()`` calls it every ``exitSignalCheckInterval`` (default 100)
       iterations, so a scheduler signal such as ``--signal=USR1@300`` leads
       to a final checkpoint.
   * - ``calculateFractionalLoadImbalance()``
     - Return the current fractional load imbalance across processes.
   * - ``doLoadBalance()`` / ``doRestructure()``
     - Rebalance / restructure the domain across processes. See
       :ref:`concepts/parallelism:Load balancing`. Requires a build with
       ParMETIS (``hemocell_parmetis``); without it ``doRestructure()`` only
       prints a warning. Not supported together with a pre-inlet or the
       CEPAC field.
