Frequently Asked Questions
==========================

Q: My cells are not initialized at the right position
-----------------------------------------------------

This is possible due to the lattice spacing usually (unless defined otherwise)
being 0.5 µm while the positions in the ``.pos`` file are in µm. This means
that a cell center of ``10 12 24`` in a ``.pos`` file will have a center of
``20 24 48`` in the simulation.

Q: Can the cells lie out of the domain?
---------------------------------------

Yes, cells are bound inside the domain by their centerpoints, so parts of the
cell might reach over the domain. These cells are deleted during initialisation,
unless the system is run on a single core.


Q: What is the exact centerpoint of the domain, how big is my domain?
---------------------------------------------------------------------

When you specify a domain of for example 10 cells, the exact
middle will be *BETWEEN* cell 5 and 6. imagine that dx is 0.5µm then the domain
is 9.5µm long, and the middle is 4.75µm. periodic boundaries of course add 1 dx
back to the length.
