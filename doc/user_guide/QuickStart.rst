
HemoCell Getting Started
========================

.. _from_source:

Setting up HemoCell from source
-------------------------------

Requirements for compiling and/or running HemoCell from source:

==========================         ==========================
Dependency                         Version
==========================         ==========================
`OpenMpi`_ or `IntelMPI`_          1.10.2 or 17.0.5
`GCC`_                             5.2.0
`CMake`_                           3.7.2
`HDF5`_                            1.8.16
`GNU Patch`_                       2.7.5
`h5py`_                            2.6.0-1
`Palabos`_                         2.0
`Parmetis`_ (optional)             4.0.3
==========================         ==========================

.. note::

  These are minimal requirements, avoid OpenMPI 2.0.X as in our experience it
  introduces memory leaks.

On Ubuntu 16.04 most of these dependencies can be installed by running::
  
  sudo apt-get install -y \
        make \
        cmake \
        g++-5 \
        g++ \
        libopenmpi-dev \
        libhdf5-dev \
        patch \
        python-h5py

This leaves the `Palabos`_ and `Parmetis`_ dependencies. Palabos can be
downloaded from their `releases
<https://gitlab.com/unigespc/palabos/-/releases>`_. After downloading Palabos
is extracted to ``./hemocell/palabos``::
  
  tar -xzf palabos-v2.0r0.tgz 
  mv palabos-v2.0r0 ./hemocell/palabos

After this Palabos must be patched. This can be done by running
``./patchPLB.sh`` from the ``./hemocell/patch/`` directory, like so::

  cd hemocell/patch && ./patchPLB.sh

The patching should succeed even though there can be an offset in some files.

For convenience, we have added a setup script ``hemocell/setup.sh`` to
automatically download and patch the Palabos library.

`Parmetis`_ can be downloaded from their `downloads
<http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download>`_. Due to the
license of parmetis we cannot distribute it with hemocell. The parmetis
download should be copied to the  ``./hemocell/external/`` directory. If you
need it because you want load balancing to be enabled you have to extract it
with::

  cd hemocell/external && tar -xzf parmetis-4.0.3.tar.gz 

Compiling HemoCell from source
------------------------------

To compile HemoCell you can navigate to a case and invoke ``cmake``
there. We recommend to do it in the following way::

  cd <path/to/case>
  mkdir build
  cd build
  cmake ../
  make

CMake can exploit parallelism when building by providing the ``-j`` flag to
indicate the number of jobs, e.g. ``cmake ../ -j 4``. It might be that CMake
fails when providing this flag, in that case, simply rerunning the command will
typically resolve the errors.

Each case depends on the HemoCell build located in ``hemocell/build/hemocell``.
Cmake is ran in this directory as a dependency of each case. It is also possible
to first run ``cmake`` in the ``hemocell/build/hemocell`` directory to first
build the library.

Furthermore a ``makefile`` is provided in the ``hemocell/examples`` directory. This
makefile can be used to update build files for all cases (see :ref:`cases_make`
for more info)

.. _packcells:

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

  cd hemocell/examples/pipeflow
  mpirun -n 4 ./pipeflow config.xml

Case output folder
------------------

The output of a case is usually written to the ``<case>/tmp`` folder. The
checkpoints are the ``.xml`` and ``.dat`` files. When a new checkpoint is
created they are moved to ``.xml.old and ``.dat.old``. The hdf5 output is stored
per timestep in ``tmp/hdf5`` and the csv output in ``tmp/csv``. See
:any:`read_output` and :any:`bpp` for more info.


.. _read_output:

Parsing the output of a HemoCell case
--------------------------------------

A HemoCell case produces multiple types of output. The simplest is the ``csv``
output which consists of all the information about cells in csv files. To merge
the csv files into a single one per timestep you can use the script :any:`ccsv`
in the ``tmp`` directory. This will generate them for you.

The more detailed ouput on both the fluid field and particle field is stored in
``hdf5`` format. We recommend using the `XDMF`_ format to make these
readable for `Paraview`_ . To generate ``*.xmf`` files run the :any:`bpp`
script.

When you have created the ``*.xmf`` files you can load them into paraview,
please select the *Legacy* XDMF file format when loading them in. The HemoCell
``.xmf`` files are not yet XDMF3 compatible.

Resuming from a checkpoint
--------------------------

To resume from a checkpoint you should run the executable from the directory you
ran it originally from (so the directory with the ``.xml`` and ``.pos`` files
visible. The first argument should be ``tmp{_x}/checkpoint/checkpoint.xml`` instead of
``config.xml``. HemoCell should then automatically resume from the last saved
checkpoint.

.. note::
  
  The number of processors used when reusing from a checkpoint does not need
  to be the same as the number of processors used for the initial run.

.. _Paraview: https://paraview.org
.. _XDMF: http://xdmf.org/index.php/Main_Page
.. _GNU Patch: https://savannah.gnu.org/projects/patch/
.. _IntelMPI: https://software.intel.com/content/www/us/en/develop/tools/mpi-library.html
.. _OpenMPI: https://www.open-mpi.org/
.. _GCC: https://gcc.gnu.org/
.. _CMake: https://cmake.org/
.. _HDF5: https://www.hdfgroup.org/
.. _h5pY: https://www.h5py.org/
.. _Parmetis: http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview
.. _Palabos: https://palabos.unige.ch/
