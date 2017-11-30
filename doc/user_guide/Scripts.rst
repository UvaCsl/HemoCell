Scripts within HemoCell
=======================

TODO

.. _cases_make:

hemocell/cases/MakeFile
-----------------------

make **cmakefiles**

  replace ``<case>/CmakeLists.txt`` with ``./CmakeLists_template.txt`` if you
  want to update the build process for all cases.

make  **executables**

  compile every case

make **all**
  
  synonym for ``make executables``

make **clean**

  remove all build directories, CMakeLists.txt and executables


hemocell/patch/patchPLB.sh
--------------------------

used to patch Palabos to add some features that HemoCell needs but are by
default not available.

hemocell/scripts/batchPostProcess.sh
------------------------------------

This scripts uses the ``*XMF.py`` files to generate all necessary xmf-files so
that the output of a job can be read into ParaView and others. This script
should be run within the ``hemocell/cases/<case>`` or
``hemocell/cases/<case>/tmp`` directory. like so::

  cd hemocell/cases/<case>
  ../../scripts/batchPostProcess.sh

hemocell/scripts/cartesius[_intel]_env.sh
------------------------------------------

These files contain the corresponding build dependencies on the cartesius system
of SurfSara. They should be sourced instead of executed::

  . ./scripts/cartesius_env.sh
