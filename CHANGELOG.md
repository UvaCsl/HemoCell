Version 0.4
===========
- Periodicity limit can now be set (hemocell::setSystemPeriodicityLimit) per
axis, default is 100, this means that a cell can wrap around a maximum of 100
times over a given axis

- Loadbalancer now reloads the old (smallest possible) atomic blocks before load
balancing, this means timing information is not available as weight, only the
number of lsps (or any other measure you come up with)

- The CmakeFile for a case can be generated from cases/CMakeLists_template.txt
by the makefile (cases/Makefile). running make in ./cases will create and build all cases.

Version 0.5
===========

- Inner links are specified in the <CellType>.xml file through
"CellMechanics"->"InnerEdges"->"Edge"
- Added MESH_FROM_STL for addCellType. Takes a triangulized STL file from xml-file->"MaterialModel"->"StlFile".
