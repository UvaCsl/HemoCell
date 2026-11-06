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

HemoCell writes the ``*.xmf`` (XDMF) files that make its HDF5 output readable
by ParaView and others while it runs (see :ref:`reference/outputs:Files`). This
script writes the same files for output of older HemoCell versions, or rewrites
them. It runs ``scripts/HDF5toXDMF.py`` (Python 3 with ``h5py``) on an output
directory, or on every output directory in a case directory (``tmp``,
``tmp_0``, ...). For instance::

  cd hemocell/examples/<case>
  ../../scripts/batchPostProcess.sh

Existing per-step files are kept unless ``--force`` is given; the files with
all steps (``RBC.xmf`` etc.) are always rewritten. Directories can be passed as
arguments.

hemocell/scripts/[_system_]_env.sh
----------------------------------

These files contain the build dependency module setup on various HPC systems.
For instance, ``snellius_env.sh`` contains the module settings to compile on the Snellius system of SURF. These scripts should be sourced instead of executed::

  . ./scripts/snellius_env.sh

.. _helper_scripts:xmf_to_x3d:

convert_xmf_to_x3d.py
----------------------------------------------------

This script converts the per-step ``*.xmf`` particle files (RBC, PLT, etc.) to the
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
