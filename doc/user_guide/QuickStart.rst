HemoCell Getting Started
===================

To get quickly started with HemoCell we recommend using our singularity [#SL]_
which can be downloaded from here: `target_for_download`_

.. todo::
  create target_for_download

Setting up HemoCell with singularity
------------------------------------

Requirements for running HemoCell with singularity:

  =========== =====================
  Dependency  Version
  =========== =====================
  Singularity 2.4 or higher
  =========== =====================

To set up HemoCell with singularity you must first have a machine with
singularity [#SL]_ installed. You need root access to a machine to install
singularity. Pay attention to the
version requirement of singularity. Singularity is still in *beta* and
regularily making breaking changes. 

After you finished setting up singularity you can download the HemoCell
singularity image from here: `target_for_download`_. 

The singularity image supports the following commands which should enable you to
compile, run and parse output from HemoCell:

singularity exec hemocell.img **make_library**
  This will copy the source files of HemoCell to the folder ``./hemocell``.
  Afterwards it will try to compile the HemoCell library within this folder
  while using the dependencies from the singularity image

singularity exec hemocell.img **compile** [*<case>*,...]
  This will try to compile the *case* within the ``./hemocell/cases/<case>`` folder.
  While compiling it will use the dependencies from the singularity image.                                                                
singularity exec hemocell.img **compile_all**
  This will try to compile all the cases within the ``./hemocell/cases/`` folder.
  While compiling it will use the dependencies from the singularity image.                                                                
singularity exec hemocell.img **run** [--optional mpirun arguments] *<case>*
  This runs the HemoCell *case* within ``./hemocell/cases/<case>``. The executable
  is called with mpirun and optional arguments (like -n 4). The output is thus
  stored in ``./hemocell/cases/<case>/tmp/``

singularity exec hemocell.img **post_process** [*<case>*,...]
  This creates the ``*.xmf`` and ``*.csv`` files in the ``./hemocell/cases/<case>/tmp`` 
  directory. These files can be used to analyze the data (e.g. with Paraview [#PF]_)

Setting up HemoCell from source
-------------------------------

Requirements for running HemoCell from source:

  +-------------+---------+
  |Dependency   |Version  |
  +=============+=========+
  |OpenMpi *or* | 1.10.2  |
  |             |         |
  |Intel Mpi    | 17.0.5  |
  +-------------+---------+
  | Gcc         | 5.2.0   |
  +-------------+---------+
  | Cmake       | 3.7.2   |
  +-------------+---------+
  | Hdf5        | 1.8.16  |
  +-------------+---------+
  |Gnu Patch    | 2.7.5   |
  +-------------+---------+
  | h5py        | 2.6.0-1 |
  +-------------+---------+

.. _target_for_download: about:blank 

.. [#PF] `https://paraview.org <https://paraview.org>`_

.. [#SL] `singularity.lbl.gov <http://singularity.lbl.gov/>`_
