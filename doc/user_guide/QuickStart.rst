
HemoCell Getting Started
===================

To get quickly started with HemoCell we recommend using our singularity [#SL]_
which can be downloaded from here: :ref:`downloads`.

HemoCell with singularity
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
singularity image from here: :ref:`downloads`. 

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
singularity exec hemocell.img **run** [--optional mpirun args] *<case>* <config>.xml
  This runs the HemoCell *case* within ``./hemocell/cases/<case>``. The executable
  is called with mpirun and optional arguments (like -n 4). The output is thus
  stored in ``./hemocell/cases/<case>/tmp/``

singularity exec hemocell.img **post_process** [*<case>*,...]
  This creates the ``*.xmf`` and ``*.csv`` files in the ``./hemocell/cases/<case>/tmp`` 
  directory. These files can be used to analyze the data (e.g. with Paraview [#PF]_). This script removes previously generated csv and xmf files

Setting up HemoCell from source
-------------------------------

Requirements for compiling and/or running HemoCell from source:

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
  | (Optional)  | 4.0.3   |
  | ParMetis    |         |
  +-------------+---------+
  | Palabos     | 2.0     |
  +-------------+---------+

On ubuntu 16.04 most of these dependencies can be installed by running::
  
  sudo apt-get install make cmake g++-5 g++ libopenmpi-dev libhdf5-dev patch python-h5py

This leaves the palabos and Parmetis dependencies. Palabos can be downloaded
from `palabos.org`_. After downloading Palabos must be extracted to ``./hemocell/palabos`` which can
be done like so::
  
  tar -xzf palabos-v2.0r0.tgz 
  mv palabos-v2.0r0 ./hemocell/palabos

After this palabos must be patched. This can be done by running
``./patchPLB.sh`` from the ``./hemocell/patch/`` directory, like so::

  cd hemocell/patch && ./patchPLB.sh

The patching should succeed even though there can be an offset in some files.

Parmetis can be downloaded from the `parmetis <http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download>`_ 
site. Due to the license of parmetis we cannot distribute it with HemoCell. 
The parmetis download should be copied to the  ``./hemocell/external/`` directory. 
If you need it
because you want load balancing to be enabled you have to extract it with::

  cd hemocell/external && tar -xzf parmetis-4.0.3.tar.gz 

Compiling HemoCell from source
------------------------------

To compile HemoCell you can simply navigate to a case and invoke ``cmake``
there. We recommend to do it in the following way::

  cd <path/to/case>
  mkdir build
  cd build
  cmake ../
  make

Cmake might sometimes fail while using the -j flag with make. Then simply try again.

Each case depends on the hemocell build located in ``hemocell/build/hemocell``.
Cmake is ran in this directory as a dependency of each case. It is also possible
to first run ``cmake`` in the ``hemocell/build/hemocell`` directory to first
build the library.

Furthermore a ``MakeFile`` is provided in the ``hemocell/cases`` directory. this
makefile can be used to update build files for all cases (see :ref:`cases_make`
for more info)

Generating initial positions for cells
--------------------------------------

At some point you might want to run a slighty different geometry, or run your
simulation with a different concentration of cells. For this we offer the
``packCells`` tool which can be found in the ``./hemocell/packCells`` directory.

This tool has a CMake file and can be build with::
  
  cd ./tools/packCells
  mkdir build && cd build
  cmake ../
  make

The result should be a ``packCells`` binary. This program offers a rich suite of
options to generate initial conditions for cells. Just type ``./packCells --help`` 
to see how it works.

The resulting ``*.pos`` files can be copied to the case where you want to use
them.


Running a HemoCell case
-----------------------

A HemoCell case should be run within the folder containing the ``.xml`` and
``.pos`` files. You can specify the number of desired processors with
``mpirun``. The only argument for the case should be the ``config.xml`` file.
A typical command looks like this::

  cd hemocell/cases/pipeflow
  mpirun -n 4 ./pipeflow config.xml


.. _read_output:

Parsing the output of a HemoCell case
--------------------------------------

a HemoCell case produces multiple types of output. The simplest is the ``csv``
output which consists of all the information about cells in csv files. To merge
the csv files into a single one per timestep you can use the script :any:`ccsv`
in the `tmp` directory. This will generate them for you.

The more detailed ouput on both the fluid field and particle field is stored in
``hdf5`` format. We recommend using the xdmf [#XDMF]_ format to make these
readable for paraview [#PF]_ . To generate ``.xmf`` files run the :any:`bpp`
script.

When you have created the ``.xmf`` files you can load them into paraview, please
select the *Legacy* xdmf file format when loading them in. the HemoCell ``.xmf``
files are not yet Xdmf3 compatible.


.. [#PF] `https://paraview.org <https://paraview.org>`_

.. [#SL] `singularity.lbl.gov <http://singularity.lbl.gov/>`_

.. [#XDMF] `Xdmf.org <www.xdmf.org>`_

.. _palabos.org: http://palabos.org
