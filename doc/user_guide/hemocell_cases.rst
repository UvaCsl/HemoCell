Cases in HemoCell
=================


Pipeflow
--------

The pipeflow case is a very generic case that can readily be used as a basis for
many other projects. The domain is a tube (loaded from ``tube.stl``) with
periodicity in the x-direction. The domain is (depending on the config files)
initialized with Red Blood Cells and Platelets with an hematocrit of 40%

oneCellShear
------------

The oneCellShear case is used for validation of the material models used in
HemoCell. A single cell is initialized in a domain with periodicity turned on in
all directions. Then the domain is sheared such that the top moves in the
positive x-direction and the bottom in the negative x-direction. The largest
diameter is reported and validated against experimental data.

.. todo::
  Reference

stretchCell
-----------

The stretchCell cases used for validation of the material models used in
HemoCell. A single cell is intitialized in a domain with periodicity in all
directions. The only external force is on the outer points of the opposing sides
of the cell. In this manner the cell is stretched. This is validated against
expermental data.

.. todo::
  Reference
