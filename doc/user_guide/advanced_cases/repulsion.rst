Repulsive forces
================

Repulsive forces between cells and walls can be enabled through the following
functions:

.. code-block:: c++

   hemocell::setRepulsion(T repulsionConstant, T repulsionCutoff);
   hemocell::enableBoundaryParticles(T boundaryRepulsionConstant, T boundaryRepulsionCutoff);

Details on their implementations can be found in ``hemocell/core/hemocell.cpp``
and ``hemocell/core/hemoCellParticleField.cpp`` under
``HemoCellParticleField::applyRepulsionForce``.

.. note::
   The repulsive forces should be enable before starting the main HemoCell
   iterations.

The first function ``setRepulsion`` enables the repulsion between cells with a
given repulsion constant and cut-off distance for the repulsion to be active.
Similarly, ``enableBoundaryParticles`` applies a repulsion force between the
cells and the boundary walls of the domain. The effect of these forces will
strongly depend on the chosen repulsion constants and cut-off distance
thresholds.

For two vertices of different cells (or a vertex and a wall node) at a distance
``d`` below the cutoff ``r_c``, the force has magnitude ``k * r_c / d`` and acts
along the line between them (equal and opposite on two vertices; at a wall
only the vertex is pushed). ``k`` is thus the force at the cutoff distance.
The repulsion force of a vertex (particle and wall repulsion together) is
capped at ``FORCE_LIMIT`` (50 pN by default), separately from the membrane
force; at ``k = 5 pN`` the cap is reached at a tenth of the cutoff.

.. note::
   The ``repulsionConstant`` and ``boundaryRepulsionConstant`` (``k``) are given
   in **pN** and the cutoffs in **micrometres**; both are converted to lattice
   units (call these functions after ``param::lbm_*_parameters()``). The
   example configurations use ``kRep = 5`` pN. (Before this version ``k`` was
   used as a lattice force; the former value ``2e-22`` had no effect.)

.. note::
   Neighbours are searched in adjacent lattice nodes only, so a cutoff larger
   than one lattice spacing is only partially effective (HemoCell prints a
   warning). Particle and wall repulsion must use the same update interval
   when both are enabled.

Alternatively, HemoCell used to provide more advanced repulsion methods
considering custom repulsion potential forces, see
``legacy/thrombosit/adhesionForces3D.h``. Although this feature is currently not
actively supported in HemoCell, one might pursue such an implementation for more
advanced force potential implementations.

