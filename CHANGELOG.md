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
- Added "Volume" option for Celltype xml files
- Added Celdensity output for fluid field
- Added preInlet possibility and cases
- Default is now HDF5 with openmpi support, if not present, using a preInlet
will generate an error when using multiple processors
- variables checkpointDirectory and outputdirectory can now be specified in the
config.xml

Version 1.4
===========

- Added hemocell.statistics as an advanced (fully recursive) profiling tool
- Added OUTPUT_INNER_LINKS to output possibility
- Consistent Changelog (At least I'll try)

Version 1.5
===========

- Added compile time constants PREINLET_MECHANICS and SOLIDIFY_MECHANICS to
  constant_defaults
