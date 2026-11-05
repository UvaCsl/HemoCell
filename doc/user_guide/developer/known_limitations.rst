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
   different cells (or a vertex and a wall node) for ``d < r_c``; ``k`` is
   given in pN and ``r_c`` in micrometres. The repulsion force per vertex is
   capped at ``FORCE_LIMIT``, separately from the membrane force, so the total
   force of a vertex can reach twice the limit. The neighbour search covers
   adjacent lattice nodes only, so cutoffs above one lattice spacing are only
   partially effective (a warning is printed). The shipped configurations use
   ``kRep = 5`` pN, a low value chosen to have an effect, not a validated one.
   Particle and wall repulsion must use the same update interval.

Force cap
   ``FORCE_LIMIT`` (default 50 pN per vertex) caps the membrane force spread to
   the fluid; a quarter of it caps the membrane-viscosity force per edge. Capping
   keeps runs stable but changes the physics when it is active (a warning is
   printed at start-up).

Cells across a periodic boundary at load
   By default, a cell of a ``.pos`` file whose surface crosses a periodic
   boundary is not placed; the log counts such cells. With
   ``setPlaceCellsAcrossPeriodicBoundaries(true)`` they are placed as periodic
   images, which is correct only for ``.pos`` files packed periodically for
   this domain: the shipped ``.pos`` files were packed for other box sizes (the
   pipeflow test file for 67 um in x, its domain is 51.5 um), and their cells at
   opposite ends would overlap.

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

* One ``HemoCell`` object at a time per process (process-wide state, reset
  when the object is destroyed).
* Library code is C++11; the test runner is C++17.
* Output is written as one HDF5 file per atomic block and output step; very
  large runs produce many small files. Most datasets are single precision.
* Restarts need the same case code and cell types; the number of ranks may
  change.
