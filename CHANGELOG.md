HemoCell Changelog
==================

In this changelog, we consider the following classes of changes:
* Features: any new features that are added to HemoCell.
* Structure: any changes to the structure of HemoCell that may break existing cases.
* Fixes: (small) changes that do not fall in the other categories.

2.6 (Unreleased)
-----------------
* Features
  * Added functionality ti hemocell.initializeLattice() to define domain decomposition from the config xml.
  * Added weak-scaling test configs for the case ``cube``.
* Fixes
  * Various comment fixes in the code.
  * Better error messages.
  * Improved sanity-check function.

2.5 (Jan 27 2022)
-----------------

* Features
  * Added ``tools/pos_to_vtk`` for simple cell visualisation.
  * Included section `Common mistakes` in documentation.
* Structure
  * Changed default `make` target to HemoCell library only.
  * Replaces the macro `#HEMOCELL_PARTICLE_FIELD` with the class directly.
  * Moved to Palabos library version `v2.3.0`.
  * Reduced Palabos patch size.
* Fixes
  * Removed packing Palabos source, now pulled from their remote by `setup.sh`.
  * Updated preinlet communication to use a local buffer as temporary storage to
    avoid velocity artefacts due to reuse of the same velocity buffer for
    communication.
  * Updated preinlet communication structure to reduce the number of messages
    send by only considering overlapping and neighbouring atomic blocks as
    targets to send messages.
  * Preliminary fix for solidification mechanics on multi-core machines.

2.4 (Sep 23 2021)
-----------------
* Features
  * Added two capillary flow examples.
  * Added, updated, and extended documentation of most example cases.
  * Added helper functions in `helper/geometry.h` for boolean operations on domain functionals.
  * Added script to convert HDF5 to the X3D format to enable rendering with Blender.
  * Added validation tests based on cell stretching and pipe flow.
* Structure
  * Cell mechanics constants are moved to `config/constant_defaults.h`.
  * The Doxygen pages are automatically generated when building the Sphinx documentation.
  * The continuous integration pipeline is split into four stages: build, test, validation, and documentation.
* Fixes
  * Added guard statement to prevent reading uninitialised Lees-Edwards displacement values.

2.3 (Jun 1 2021)
----------------
* Features
  * Rebuilt the build environment. For more details, please see the quickstart guide in the documentation.
  * Added the cube example in which a cubic domain is subjected to a shear flow on the top/bottom planes.
* Structure
  * HemoCell now depends on Palabos 2.2.1 as the older version used previously is no longer available. The patch and setup.sh file has been updated to reflect this.
* Fixes
  * Remove some useless output in the cellStretch example

2.2 (Dec 7 2020)
----------------
* Features
  * Added support for Lees-Edwards boundary conditions (single core only for now, example soon to come with publication).
  * Significantly improved and updated the documentation.
  * Publications using HemoCell added to documentation.
  * Updated the website with more examples in the gallery.
  * Added further output fields: e.g. shear rate, strain rate, shear stress.
  * Several small updates and QoL improvements in the otput processing scripts.
* Fixes
  * Fix a bug in the velocity output when SI units were used

2.1 (Feb 6 2020)
----------------
* Features
  * Added support for binding sites
  * Added support for interior viscosity
  * Added velocity output for hemoCellParticleField
  * Added curved pipeflow with pre-inlet example
  * Added microcontraction and flowaroundsphere examples
  * Added support for pulsatility
* Structure
  * Rename HemoCell::setMinimumDistanceFromSolid to HemoCell::setInitialMinimumDistanceFromSolid
  * Rename all occurences of RBC_HO to RBC, both in filenames and source code
  * Merge csv and hdf5 reader for scripts
  * The preinlet setup now adjusts the preinlet to always be on the lattice, if the lattice is set up before the preinlet (a warning is given otherwise)
* Fixes
  * Do not override interrupt handlers until the first iteration has started
  * Added missing AGPL license header to several files
  * Added profiler statements to preInlet.cpp
  * Log an error message whenever a nonexisting cell type is requested
  * Fix checkpointing so that it works with the pre-inlet
  * Update scripts so that they are compatible with the new RBC output
  * Fix outputStrainRate in IO/FluidHdf5IO.hh
  * Fix bug in the pre-inlet when a block boundary is spanned
  * Fix computation of the Tresca value
  * Fix location of Palabos download page
  * Fix for the environment script on Cartesius
  * Set the driving force at the warmup part for the pipeflow example

2.0b (Feb 12 2019)
------------------
to be added

2.0 (Dec 17 2018)
-----------------
* Added compile time constants PREINLET_MECHANICS and SOLIDIFY_MECHANICS to constant_defaults
* Added runtime preinlets, see stl_preinlet for examples, instead of preinletFromSlice() you can call autoPreinletFromBoundary() if the boundary happens to be the inlet as well
* Interior viscosity edge tracking works again, enable interior viscosity through CELL.xml file, see cases/RBC_template.xml for an example
* Output is saved in a seperate folder for each run, see pipeflow/config.xml for the options
* CSV output is already concatenated in csv folder, now also includes base_cell_id, correct cell_id and velocity

1.4 (Jul 9 2018)
----------------
* Added hemocell.statistics as an advanced (fully recursive) profiling tool
* Added OUTPUT_INNER_LINKS to output possibility

0.5
---
* Inner links are specified in the <CellType>.xml file through "CellMechanics"->"InnerEdges"->"Edge"
* Added MESH_FROM_STL for addCellType. Takes a triangulized STL file from xml-file->"MaterialModel"->"StlFile".
* Added "Volume" option for Celltype xml files
* Added Celdensity output for fluid field
* Added preInlet possibility and cases
* Default is now HDF5 with openmpi support, if not present, using a preInlet will generate an error when using multiple processors
* variables checkpointDirectory and outputdirectory can now be specified in the config.xml

0.4
---
* Periodicity limit can now be set (hemocell::setSystemPeriodicityLimit) per axis, default is 100, this means that a cell can wrap around a maximum of 100 times over a given axis
* Loadbalancer now reloads the old (smallest possible) atomic blocks before load balancing, this means timing information is not available as weight, only the number of lsps (or any other measure you come up with)
* The CmakeFile for a case can be generated from cases/CMakeLists_template.txt by the makefile (cases/Makefile). running make in ./cases will create and build all cases.

