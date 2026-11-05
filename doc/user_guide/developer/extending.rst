Extending HemoCell
==================

Recipes for the most common extensions. Keep reusable behaviour in the library
(``core/``, ``mechanics/``, ``helper/``, ``io/``); cases only assemble geometry,
boundary conditions, materials and run schedules.

Adding a mechanical model
-------------------------

#. Start from the closest model in ``mechanics/`` (e.g. ``RbcHighOrderModel``)
   and derive from ``hemo::CellMechanics`` (``mechanics/cellMechanics.h``). The
   constructor takes ``(Config & materialConfig, HemoCellField & cellField)``;
   ``HemoCell::addCellType<Model>(name, shape)`` creates it.
#. Implement ``ParticleMechanics(particles_per_cell, lpc, ctype)`` (both maps
   keyed by the ``plint`` cell ID; mark it ``override``): loop over
   the cells in ``lpc`` (complete cells of this block), skip other cell types
   (``cell[0]->sv.celltype != ctype``) and add forces through the component
   pointers (``*particle->force_link += f``), never to ``sv.force`` directly
   (see :ref:`developer/architecture:Force components`). Implement
   ``statistics()`` to log the model parameters.
#. Use ``cellConstants`` (``CommonCellConstants``: triangle and edge lists,
   reference lengths, areas, volume) instead of recomputing reference
   quantities, and read material parameters once in the constructor.
#. Keep the loop allocation-free (see
   :ref:`developer/conventions:Performance rules`). Check force balance:
   internal forces of a free cell must sum to zero (net force and torque).
#. Reconfigure CMake (sources are globbed), and add a small validation (e.g. a
   stretch or shear test, see ``tests/validation``) and a material template in
   ``tools/cell_templates``. Document the model in
   ``doc/user_guide/reference/mechanical_models.rst``.

Adding an output field
----------------------

#. Add an ``OUTPUT_*`` constant (``constexpr int`` in ``namespace hemo``) in
   ``config/constant_defaults.h`` (do not renumber existing ones; case files
   and output files use them).
#. Cell fields: add an ``output*`` member to ``HemoCellParticleField``
   (``io/hemoCellParticleFieldOutputFunctions.cpp``) and register it in
   ``AddOutputMap()``. Fluid fields: add a ``case`` to ``WriteFluidField``
   (``io/FluidHdf5IO.hh``).
#. Output must not change the simulation state; restore anything you modify.
#. Datasets written through these paths are listed in the XDMF files
   automatically (cell fields: one row per vertex of ``getOutputCells()``;
   fluid fields: ``(z, y, x, components)`` over the bulk).
#. State the dataset name, shape, units (SI or lattice units, see
   ``outputInSiUnits``) and precision in
   ``doc/user_guide/reference/outputs.rst``, and update
   ``scripts/HDF5toXDMF.py`` and the analysis scripts in ``scripts/`` if they
   need to know the new dataset.

Adding a configuration setting
------------------------------

* Read it once, at setup, with ``(*cfg)["section"]["name"].read<T>()``. A
  required setting lets the ``std::invalid_argument`` propagate; an optional one
  catches it and uses a documented default.
* Validate ranges and cross-field constraints before the first time step
  (``sanityCheck()`` or the setter), not inside the time loop.
* Keep existing tag names; document the tag (with units) in
  ``doc/user_guide/xml_files.rst`` and update the example/test configurations
  that need it. Note that settings like ``tmax`` or ``tmeas`` are read by the
  case, not by the library.

Adding a case
-------------

Copy ``examples/template`` (or start from a similar example), name the source
after the directory, add the directory to the list in
``examples/CMakeLists.txt`` (or ``cases/CMakeLists.txt``) and its executable
path to ``.gitignore``. Link ``hemocell`` or the variant library the case needs
(``hemocell_interior_viscosity``, ``hemocell_solidify_mechanics``).

Adding tests
------------

All tests live under ``tests/`` and build into one GoogleTest runner
(``hemocell_test``); new files are picked up after reconfiguring.

* Unit tests (``tests/unit``) exercise a component without a full simulation
  and run in seconds; they belong to the CTest label ``unit``.
* Integration and validation tests construct a ``HemoCell`` with
  ``HemoCell::MPIHandle::External`` (MPI is initialised once in
  ``tests/main.cpp``). Name physics validations ``Validation...`` so the
  ``validation`` label picks them up.
* Tests run in the repository root and read their inputs relative to it.
  Construct simulations from ``hemo_test::testConfig(path)``
  (``tests/test_support.h``; ``pipeflowConfig()`` for the pipeflow setup), a
  copy of the configuration whose output goes to ``hemo_test::outputRoot()``,
  one temporary directory per run that is removed afterwards (set
  ``HEMOCELL_KEEP_TEST_OUTPUT`` to keep it). Other files a test writes go there
  as well, with the MPI rank in the name where ranks write separately.
* ``tests/integration/test_mechanics.cpp`` checks the membrane forces of every
  cell model (rest shape, net force and torque, rotation); register a new model
  there.
* Assert behaviour (values, invariants, equality between 1 and N ranks), not
  just that the code runs. Derive tolerances from the physics; do not widen
  them to make a test pass.
