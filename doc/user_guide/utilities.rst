Utility tools
=============

In addition to helper scripts (:ref:`helper_scripts:Helper scripts`), HemoCell provides
different utility tools under the ``tools/`` directory. Typically these provide
more extensive functionality compared to the more simple helper scripts in
``scripts/``.

.. _pos_to_vtk:

Viewing of cell packings
------------------------

When you are willing to inspect a distribution of cells, for instance as result
of running the ``packCell`` tool, the ``pos_to_vtk`` command-line utility
provides direct visualisation using using the corresponding cell position
``pos``-files only. Default usage illustrates the cells positions, considering
red blood cells, platelets, or a combination thereof, in a small GUI window.
This allows to explore the cell distributions in an easy manner up to 100000s of
cells. Furthermore, the script can convert the ``pos``-files towards output
files readable by Paraview, for when the cell counts are not properly handled
anymore by the simple GUI display.

The utility is easily installed within a virtual Python environment

.. code::

   pip install .

   # or when editing the source code

   pip install --editable .

Usage is straightforward, for instance

.. code::

   # Show full usage and documentation
   pos_to_vtk --help

   # Show a RBC distribution
   pos_to_vtk /path/to/RBC.pos

   # Or pass through a pipe (the dash (-) is required)
   cat /path/to/RBC.pos | pos_to_vtk -

   # Combine RBC distributions
   pos_to_vtk /path/to/RBC.pos /other/path/to/RBC.pos

   # Show RBC and PLT
   pos_to_vtk /path/to/RBC.pos --plt /path/to/PLT.pos

   # Export to XDMF for visualisation in, for example, Paraview
   pos_to_vtk /path/to/RBC.pos --output /path/to/output.xdmf

   # Clip cells to a bounding box
   pos_to_vtk /path/to/RBC.pos --bounding-box 100 50 50 --clip

   # Clip cells to a surface mesh (STL)
   pos_to_vtk /path/to/RBC.pos --stl /path/to/mesh.stl --clip
