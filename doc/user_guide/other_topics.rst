Other topics
============

Loading a cell from an stl file
-------------------------------

Usually a cell is defined to be generated mathematically, this means that the
vertex locations are gathered from describing equations. The place where the
method for generating the cell is defined is in the ``hemocell.addCellType<>()``
function. 

To load a celltype from stl you have to change the second argument to
``MESH_FROM_STL``. This means that HemoCell will look at the stl file located at
the place specified in the ``CELL.xml`` file. The specification of this file
should be in the ``<hemocell><MaterialModel><StlFile>`` tag. 

.. note::
  
  The stl file should only contain triangles. The stl file reader that is used
  by HemoCell is the one already present in  Palabos

  - This reader has no capabilities to triangulate squares and will
    crash on stl files that do not follow this rule.
  - Furthermore the STL file
    needs to be in ASCII format. Nowadays most readers save in Binary format by
    default.

Running a pure fluid flow (without cells)
-----------------------------------------

With HemoCell it is possible to run a pure fluid flow, meaning that no cells are present within the domain.
This is actually the same as just running a fluid simulation with Palabos. The simplest way to remove all the cells is to set the number of cells in the ``CELL.pos`` files to zero (first line). Or to remove the addCellType and all subsequent functions that request that specific celltype (such as setOutputs) from the ``case.cpp`` file. 

However this still leaves a small overhead as the function for the material model are still called (although little work is done inside them as there are no cells). This can be circumvented by replacing 

+----------------------+------+------------------------------------------+
| .. code-block:: c++  |      | .. code-block:: c++                      |
|                      |      |                                          |
|   hemocell.iterate() | With |   hemocell.lattice->collideAndStream();  |
|                      |      |   hemocell.iter++;                       |
+----------------------+------+------------------------------------------+

within the ``case.cpp`` file

Saving only the CSV output and not the HDF5 output
--------------------------------------------------

This is possible with the extra ``<sim><tcsv>`` parameter which is already added in the
:ref:`cases/pipeflow:Pipe flow` case. This parameter then controls a seperate call to writeCellInfo_CSV that writes only CSV output.
To reduce the HDF5 output you can simply increate ``<sim><tmeas>`` in the config
file.

.. note::

  Don't forget to include the right header:

.. code-block:: c++

  #include "writeCellInfoCSV.h"

