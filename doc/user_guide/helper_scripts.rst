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


.. _helper_scripts:xmf_to_x3d:

convert_xmf_to_x3d.py
----------------------------------------------------

This script converts generated ``*.xmf`` particle files (RBC, PLT, etc.) to the
X3D format for rendering purposes, e.g. using :ref:`visualization:Blender`. The
script creates a ``x3d`` directory in the output directory of an example case
and populates this with the ``*.x3d`` files. The script is used as follows

.. code-block:: text

   usage: convert_xmf_to_x3d.py [-h] [--view-size width height]
                                [--smooth SMOOTH]
                                path [path ...]

   positional arguments:
     path                  Directories containing XMF output files to
                           convert to X3D scenes

   optional arguments:
     -h, --help            show this help message and exit
     --view-size width height
                           Set the render view-size in Paraview
     --smooth SMOOTH       Number of mesh smoothing iterations
