Example cases
=============

The following examples illustrate various use cases of HemoCell. The examples
range from relatively simple single cell problem definitions, to more complex
Bulk flow problems. The single cell problems, such as a single
:ref:`shearing <cases/onecellshear:One shearing cell>` or
:ref:`stretching <cases/stretchcell:One stretching cell>` focus more on the
material models and behaviour of individual cells in specific. These can be
excellent starting points when implementing or modifying specific cell
mechanics.

A next step might then be the analysis of single cell's deformation under straight-channel flow conditions, which is illustrated by the :ref:`parachuting
<cases/parachuting:A parachuting cell>` example followed by cell-cell
interactions and collisions between multiple cells
in :ref:`cases/cellCollision_interior_viscosity:Colliding cells with interior viscosity`.
Note that due to the fluid solid coupling scheme, the immersed boundary method, the cell membranes already touch when they are half IBM kernel-width apart.

Afterwards, we present to variations on a pipe flow condition using higher cell
counts and how to use the ``pre-inlet`` boundary conditions. These examples can
be excellent starting positions for typical simulations.

.. toctree::
   :maxdepth: 1

   cases/onecellshear.rst
   cases/stretchcell.rst
   cases/parachuting.rst
   cases/cellCollision_interior_viscosity.rst
   cases/pipeflow
   cases/pipeflow_with_preinlet.rst
