Scripts within HemoCell
=======================

.. _cases_make:

hemocell/examples/makefile
--------------------------

Multiple helper commands are accessable through ``make``: 

* make **cmakefiles**

  replace ``<case>/CmakeLists.txt`` with ``./CmakeLists_template.txt`` if you
  want to update the build process for all cases.

* make  **executables**

  compile every case

* make **all**
  
  synonym for ``make executables``

* make **clean**

  remove all build directories, CMakeLists.txt and executables


hemocell/patch/patchPLB.sh
--------------------------

Patch Palabos to add some features that HemoCell needs but are by default not
available.

.. _bpp:

hemocell/scripts/batchPostProcess.sh
------------------------------------

This scripts uses the ``*XMF.py`` files to generate all necessary xmf-files so
that the output of a job can be read into ParaView and others. This script
should be run within the ``hemocell/examples/<case>`` or
``hemocell/examples/<case>/tmp`` directory. like so::

  cd hemocell/examples/<case>
  ../../scripts/batchPostProcess.sh

hemocell/scripts/cartesius[_intel]_env.sh
------------------------------------------

These files contain the corresponding build dependencies on the Cartesius system
of SurfSara. They should be sourced instead of executed::

  . ./scripts/cartesius_env.sh

.. _ccsv:

hemocell/scripts/CellInfoMergeCSV.sh
------------------------------------

This script merges the CSV output from multiple processors into a single one in
the current directory. Use it in the ``tmp`` directory like this::

  cd hemocell/examples/<case>/tmp/
  . ./scripts/CellInfoMergeCSV.sh
