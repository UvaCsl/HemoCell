HemoCell Changelog
==================

In this changelog, we consider the following classes of changes:
* Features: any new features that are added to HemoCell.
* Structure: any changes to the structure of HemoCell that may break existing cases.
* Fixes: (small) changes that do not fall in the other categories.

UNRELEASED
----------
* Features
  * Added support for binding sites
  * Added support for interior viscosity
  * Added velocity output for hemoCellParticleField
  * Added curved pipeflow with pre-inlet example
* Structure
  * Rename HemoCell::setMinimumDistanceFromSolid to HemoCell::setInitialMinimumDistanceFromSolid
  * Rename all occurences of RBC_HO to RBC, both in filenames and source code
  * Merge csv and hdf5 reader for scripts
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
to be added

