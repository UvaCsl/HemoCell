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
only the vertex is pushed).

.. note::
   The ``repulsionConstant`` and ``boundaryRepulsionConstant`` (``k``) are used
   as given, in **lattice force units** (one lattice force unit is
   ``rho * dx^4 / dt^2``, e.g. about 6.4e-9 N for ``dx = 0.5 um``,
   ``dt = 1e-7 s`` and plasma density). The cutoffs are given in
   **micrometres** and converted to lattice units. The value ``2e-22`` found in
   several example configurations is negligible in lattice units.

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

