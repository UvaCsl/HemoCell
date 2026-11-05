Known limitations
=================

Behaviour that is known to be incomplete, untested or surprising. Check this
list before relying on a feature in production, and update it when you fix or
discover something.

Features with restrictions
--------------------------

Load balancing
   ``doLoadBalance()``/``doRestructure()`` need the ParMETIS build
   (``hemocell_parmetis``), which currently does not compile: the
   per-block timers it reads were dropped from the Palabos patch. Without
   ParMETIS, ``doRestructure()`` only prints a warning. It is not supported
   together with a pre-inlet or the CEPAC field, and has not been tested since
   its fixes.

Lees-Edwards boundary condition
   Only correct when the blocks on both sides of the shearing (z) boundary are
   on the same rank: the shear offset is applied only to particles copied
   locally, and the population buffers are per rank.

Pre-inlet
   Only the ``Xpos`` direction has been tested. The pre-inlet needs at least two
   ranks.

Solidification
   Binding sites near block boundaries and regions that already overlap cells
   when marked as binding sites have shown problems on several ranks (see
   ``TODO.md``).

Interior viscosity
   Compiled only into ``hemocell_interior_viscosity``; not covered by the test
   suite.

``cases/``
   Research cases, some experimental or incompletely validated. Validate before
   treating them as reference behaviour.

Numerical and modelling notes
-----------------------------

Repulsion
   The force is ``k * r_c / d`` along the line between two vertices of
   different cells (or a vertex and a wall node) for ``d < r_c``; ``k`` is used
   as given, in lattice force units, and ``r_c`` is given in micrometres. The
   neighbour search covers adjacent lattice nodes only, so cutoffs above one
   lattice spacing are only partially effective (a warning is printed). The
   value ``kRep = 2e-22`` in the shipped configurations is negligible in lattice
   units; check the intended magnitude before using repulsion. Particle and
   wall repulsion must use the same update interval.

Force cap
   ``FORCE_LIMIT`` (default 50 pN per vertex) caps the membrane force spread to
   the fluid; a quarter of it caps the membrane-viscosity force per edge. Capping
   keeps runs stable but changes the physics when it is active (a warning is
   printed at start-up).

Cell IDs in HDF5 output
   ``Cell Id`` is written as a 32-bit float, exact up to 16,777,216. IDs grow
   by a multiple of the number of cells at every periodic crossing (see
   :ref:`developer/architecture:Cell IDs and periodicity`), so long periodic
   runs exceed that (200,000 cells after 84 crossings in x) and the written
   IDs are rounded. The CSV files hold exact cell and base cell IDs.

Broken cells
   A cell that loses a vertex (it enters a wall node or leaves the domain, or
   its force is not finite) is deleted as a whole at the next particle update
   and not re-inserted; ``<verbose><cellsDeletedInfo>`` prints each deletion.

Time-scale separation
   Output shows the total force from the last material update but evaluates
   force components for the current shape; with separations > 1 the components
   need not add up to the total.

Initial placement
   ``setInitialMinimumDistanceFromSolid()`` rejects vertices within
   ``floor(d / dx)`` lattice nodes of a wall; distances below one lattice
   spacing have no effect. A cell is dropped as a whole when one of its
   vertices is rejected, lies outside the domain or across a periodic boundary
   (images are not placed); ``loadParticles()`` logs how many cells of each
   type were placed. A missing ``<name>.pos`` file is a fatal error; a file
   containing ``0`` means no cells of that type.

Software constraints
--------------------

* One ``HemoCell`` object per process (process-wide state).
* Library code is C++11; the test runner is C++17.
* ``config/constant_defaults.h`` defines ``SOLIDIFY_MECHANICS`` for every
  build, so ``hemocell_solidify_mechanics`` equals ``hemocell``.
* Output is written as one HDF5 file per atomic block and output step; very
  large runs produce many small files. Most datasets are single precision.
* Restarts need the same case code and cell types; the number of ranks may
  change.
