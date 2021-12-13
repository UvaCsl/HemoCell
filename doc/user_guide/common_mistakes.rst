Common mistakes
===============


Numerical instabilities
-----------------------

Sometimes the simulation outcome might show strange behaviour or unexpected
results. Often these are caused by numerical instabilities occurring somewhere
in the simulation domain. Typical sign of such instabilities are:

- Jiggling cells
- Too large density differences
- Traveling waves in the fluid field
- Strange or extreme cell deformations
- Oscillating (cell) forces or displacements

When observing any, or similar, of those signs, you might consider to
investigate the lattice Boltzmann parameters used during the simulation.
Typically you should consider ``u_lbm < 0.1`` and ``tau > 0.5``.

Cell deletions
--------------

If you observe unexpected cell removals, these might be caused due to
instabilities (see above) as the cell might undergo extreme forces or flow
velocities. When this happens, the cells are removed from the simulation in
attempt to preserve stability. Alternatively, when using too large time-steps,
the cells might leave their domain, without proper communication to the
respective neighbour or periodic in/outlet.

When cells are removed directly after initialisation, this is mostly due to the
fact that the cells extend outside of the simulation domain. For multi-core
simulations, the cells are removed, such that each cell is initialised only
once, i.e. inside the atomic block that contains the cell's center
(see :ref:`FAQ:Frequently Asked Questions`).
