Helper scripts
==============

A number of helper scripts for working with HemoCell.

.. _patching-palabos:

hemocell/patch/patchPLB.sh
--------------------------

HemoCell is build on top of Palabos and has added a number of additional
features to the source of Palabos. These additional features are contained as a
`patch <https://en.wikipedia.org/wiki/Patch_(Unix)>`_ file in
``hemocell/patch/palabos.patch`` which is applied by the
``hemocell/patch/patchPLB.sh`` script.

The script can be evaluated from the ``hemocell/patch`` directory as

.. code::

   ./patchPLB.sh

.. _bpp:

hemocell/scripts/batchPostProcess.sh
------------------------------------

This scripts uses the ``*XMF.py`` files to generate all necessary xmf-files so
that the output of a job can be read into ParaView and others. This script
should be run within the ``hemocell/examples/<case>`` or
``hemocell/examples/<case>/tmp`` directory. For instance::

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
