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

.. note::
   The ``repulsionConstant`` and ``boundaryRepulsionConstant`` are to be
   supplied in lattice units and are internally converted to SI units.

Alternatively, HemoCell used to provide more advanced repulsion methods
considering custom repulsion potential forces, see
``legacy/thrombosit/adhesionForces3D.h``. Although this feature is currently not
actively supported in HemoCell, one might pursue such an implementation for more
advanced force potential implementations.

