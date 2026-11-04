Architecture and internals
==========================

This page describes how HemoCell works inside: what one time step does, how
particles are stored and owned, how they are communicated, and what a checkpoint
contains. It complements the user-level :ref:`concepts <concepts>` and assumes
familiarity with Palabos multi-blocks and processing functionals. File and
function names refer to the repository layout described in ``AGENTS.md``.

Object structure
----------------

``hemo::HemoCell`` (``hemocell.h``, ``core/hemoCell.cpp``) owns one simulation:

* ``lattice`` -- the Palabos ``MultiBlockLattice3D<T,DESCRIPTOR>`` (D3Q19, BGK
  with Guo forcing). With a pre-inlet there are two lattices
  (``preinlet_lattice``/``domain_lattice``); ``lattice`` is the one this rank
  works on.
* ``cellfields`` -- ``hemo::HemoCellFields`` (``core/hemoCellFields.*``): all
  cell types and the particle field.
* ``cfg`` -- the configuration (``config/config.*``), and ``preInlet``,
  ``loadBalancer``.

``HemoCellFields`` holds

* ``cellFields`` -- one ``HemoCellField`` per cell type (``core/hemoCellField.*``):
  reference mesh, material configuration (``<name>.xml``), constitutive model
  (``mechanics``, a ``CellMechanics``), output selection and material time-scale
  separation.
* ``immersedParticles`` -- a ``MultiParticleField3D<HemoCellParticleField>`` with
  the same block decomposition as the lattice and an envelope of
  ``<particleEnvelope>`` nodes.

``HemoCellParticleField`` (``core/hemoCellParticleField.*``) is one atomic block
of particles. It stores the particles in ``std::vector<HemoCellParticle>
particles``, points to the matching fluid block (``atomicLattice``) and knows
its bulk region (``localDomain``, block-local coordinates).

``HemoCellParticle`` (``core/hemoCellParticle.h``) is one membrane vertex. All
state that travels between ranks is in ``sv`` (``serializeValues_t``): velocity,
position, membrane force, repulsion force, cell ID, vertex ID, residence time,
cell type (and, depending on compile options, ``vPrevious``/``solidify``).
``sv`` is copied as raw bytes for MPI and checkpoints, so it must remain
trivially copyable; changing it changes the checkpoint format.

Positions are absolute lattice coordinates (not block-local); velocities and
forces are in lattice units. Output converts to SI.

Only one ``HemoCell`` object may exist at a time: unit conversion
(``param::``/``hemo::Parameters``), ``hemo::global``, the log (``hlog``) and
``HemoCellParticleField::cellFields`` are process-wide.

One time step
-------------

``HemoCell::iterate()`` performs, in this order (``core/hemoCell.cpp``):

#. **Bookkeeping.** On the first call: install the exit-signal handlers and run
   ``sanityCheck()``; whenever the particle field was (re)created: build the
   large-envelope communication structure (see below). Every
   ``exitSignalCheckInterval`` iterations: collective exit-signal check.
#. **Repulsion** (if enabled and due): particle-particle repulsion resets and
   recomputes ``sv.force_repulsion``; wall repulsion adds to it.
#. **Force spreading.** For every particle (including envelope copies) compute
   the IBM kernel (``kernelMethod``, default ``interpolationCoefficientsPhi2``:
   2x2x2 nodes, wall nodes skipped, weights renormalised), cap ``sv.force`` at
   ``FORCE_LIMIT`` (in place), and add ``(force + force_repulsion) * weight``
   to the external force of the kernel nodes.
#. **Fluid.** ``lattice->collideAndStream()`` (and the CEPAC field, if enabled).
#. **Velocity interpolation** (every ``particleVelocityUpdateTimescale``
   iterations): interpolate ``sv.v`` from the kernel nodes computed in step 3,
   then ``syncEnvelopes()``.
#. **Solidification** (optional, if due).
#. **Advection.** ``position += v`` (Euler; Adams-Bashforth with
   ``HEMOCELL_MATERIAL_INTEGRATION=2``). Vertices whose nearest node is a wall
   node are deleted; their cell becomes incomplete.
#. **Constitutive model.** For each cell type whose material time scale is due,
   reset ``sv.force`` and evaluate ``ParticleMechanics`` for all complete cells
   in the block (envelope copies included).
#. **Interior viscosity** (variant ``hemocell_interior_viscosity``, if due).
#. **Clean-up** (velocity-update steps): drop particles further than 3 nodes
   outside the bulk; with ``<cellsDeletedInfo>`` also delete incomplete cells.
#. **Reset** the lattice external force to zero and increment ``iter``.

Because of the last step a case must re-apply any body force after every
``iterate()``. ``writeOutput()`` does not change the state used by the next
step: it restores the lattice force after writing the fluid force field and
evaluates force components into separate arrays (see `Force components`_).

Particle ownership and caches
-----------------------------

A particle is *local* to the block whose bulk contains its position (with a
half-node tolerance, ``isContainedABS``); copies of it live in the envelopes of
neighbouring blocks. Envelope copies are needed because membrane forces are
spread by every block onto its own nodes, and constitutive forces need complete
cells.

``HemoCellParticleField`` keeps lookup caches over ``particles``:

``particles_per_cell`` (ppc)
   ``cellId -> vector of particle indices`` ordered by vertex ID, ``-1`` for a
   missing vertex. A cell is *complete* in a block if no entry is ``-1``.
``lpc``
   Cells with at least one local (bulk) particle.
``particles_per_type`` (ppt)
   Particle indices per cell type.
particle grid (pg, ``core/particleGrid.*``)
   Particles binned per lattice node (compressed, no capacity limit), used by
   repulsion and solidification.

They hold **indices** into ``particles``. Any insertion, removal, reordering
(``removeParticles`` swaps with the last element) or position change makes them
stale; the ``*_up_to_date`` flags and ``invalidate_*()`` functions track this
and the ``get_*()`` accessors rebuild lazily. Pointers to particles
(``HemoCellParticle *``) are additionally invalidated by vector growth. New code
that mutates ``particles`` must invalidate the affected caches; code that keeps
indices or pointers must not mutate ``particles`` in between.

Force components
----------------

The membrane-force pointers of a particle (``force_area``, ``force_link``,
``force_bending``, ``force_volume``, ``force_visc``, ``force_inner_link``)
normally all point to ``sv.force``, so a model adds every term into the total
force. Mechanical models must therefore write through these pointers (``*p->
force_link += f``), never to ``sv.force`` directly.

For output, ``separateForceVectors()`` stores the total in ``force_total`` and
points each component to its own array; a forced evaluation of the
constitutive model fills them; ``unifyForceVectors()`` frees the arrays and
points everything back to ``sv.force``. ``sv.force`` itself is untouched, so
output does not alter the simulation.

Cell IDs and periodicity
------------------------

Cells get IDs ``0 .. N-1`` in the order of the ``.pos`` files (``N`` is
``HemoCellFields::number_of_cells``, stored in checkpoints). When a cell is
copied across a periodic boundary, its copy gets a shifted ID ``id + k N`` so
that the two images are distinct cells; ``base_cell_id()`` maps back. The
shift per crossing is ``+-1`` (x), ``+-periodicity_limit[0]`` (y) and
``+-periodicity_limit[0] * periodicity_limit[1]`` (z), default ``1``, ``100``,
``10000`` (``setSystemPeriodicityLimit()``). IDs are ``int``-sized in several
places, so ``N`` times the largest shift must fit in an ``int``; with the
defaults this limits a run to about 212,000 cells (checked in
``sanityCheck()``).

Communication
-------------

``HemoCellFields::syncEnvelopes()`` refreshes the envelopes:

#. Every block drops its non-local particles.
#. Palabos' ``duplicateOverlaps`` with ``modif::hemocell`` copies particles in
   the envelope (``HemoCellParticleDataTransfer``: ``sv`` as raw bytes;
   positions and cell IDs shifted across periodic boundaries). Before the first
   time step the envelope is the full ``<particleEnvelope>``; afterwards it is
   3 nodes.
#. After the first time step, cells are completed over the full particle
   envelope: each rank sends the IDs of the cells it holds to its large-envelope
   neighbours (tag 24) and receives all their particles of those cells (tag
   42), on a private communicator, receiving requests per known source rank.
   Copies between blocks on the same rank use ``attribute()``.

HemoCell's processing functionals (``HemoCellFunctional``) declare that they
modify nothing (``modif::nothing``), so Palabos never triggers envelope updates
on its own; HemoCell decides when to synchronise. Palabos is patched
(``patch/palabos.patch``) to add ``modif::hemocell`` and to route particle data
through ``HemoCellParticleDataTransfer``.

Pitfalls that cost time in the past:

* Access through a multi-block is collective in a parallel build:
  ``MultiScalarField3D::get()`` broadcasts the value from its owner, and a cell
  returned by ``MultiBlockLattice3D::get()`` broadcasts on queries such as
  ``computeVelocity()``. Calling them on some ranks only deadlocks. Moreover,
  ``getDynamics()`` of such a cell is a wrapper, so ``isBoundary()`` does not
  describe the node. In per-rank code use the atomic block
  (``getComponent(blockId).get(local coordinates)``).
* Do not match messages with ``MPI_ANY_SOURCE`` on a communicator that other
  code also uses: a rank that finished an exchange early may already send the
  next one. Use a private communicator (``MPI_Comm_dup``) and specific sources.
* Every rank must call collective operations (output, statistics, checkpoints,
  ``CellInformationFunctionals``) in the same order, even if only rank 0
  writes.

The pre-inlet (``helper/preInlet.*``) runs on separate ranks with its own
lattice and particle field; its velocities and cells are passed to the main
domain each step with one message per rank pair, using an exchange plan built
on the first ``applyPreInlet()``.

Checkpoints
-----------

``saveCheckPoint()`` writes to ``<outputDirectory>/<checkpointDirectory>/``:
``checkpoint.xml`` (the configuration under ``<Checkpoint><hemocell>`` plus
``<General>``: iteration, output directory, number of cells, last output
iteration), the lattice and particle field (``lattice.*``,
``particleField.*``, ``PRE_*`` for a pre-inlet), and the CEPAC, binding-site or
interior-viscosity fields when enabled. Previous files are renamed to
``*.old`` first. Case-level state (e.g. a body force) is the case's
responsibility: restart a case by passing ``checkpoint.xml`` as its
configuration and following its ``cfg->checkpointed`` branch.
