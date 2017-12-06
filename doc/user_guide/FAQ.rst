Frequently Asked Questions
==========================

Q: My cells are not initialized at the right position
-----------------------------------------------------

This is possible due to the lattice spacing usually (unless defined otherwise)
being 0.5 µm while the positions in the ``.pos`` file are in µm. This means
that a cell center of ``10 12 24`` in a ``.pos`` file will have a center of
``20 24 48`` in the simulation.
