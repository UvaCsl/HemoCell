Other topics
============

Loading a cell from an stl file
-------------------------------

Usually a cell is defined to be generated mathematically, this means that the
vertex locations are gathered from describing equations. The place where the
method for generating the cell is defined is in the ``hemocell.addCellType<>()``
function. 

To load a celltype from stl you have to change the second argument to
``MESH_FROM_STL``. This means that hemocell will look at the stl file located at
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
