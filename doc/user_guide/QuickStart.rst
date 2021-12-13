
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

This leaves the `Palabos`_ and `Parmetis`_ dependencies. For Palabos we
currently use the ``v2.2.1`` version with an additional patch. To automatically
download this specific Palabos version and apply our patch, the setup script
``hemocell/setup.sh`` can be evaluated::

  # from ./hemocell/
  ./setup.sh

You can also op to manually download Palabos through their `releases
<https://gitlab.com/unigespc/palabos/-/releases>`_. After downloading Palabos
should be extracted to ``./hemocell/palabos``::

  tar -xzf palabos-v2.2.1.tar.gz
  mv palabos-v2.2.1 ./hemocell/palabos

After this Palabos must be patched, see :ref:`patching-palabos`. This can be
done by running ``./patchPLB.sh`` from the ``./hemocell/patch/`` directory, like
so::

  cd hemocell/patch && ./patchPLB.sh

The patching should succeed even though there might be an offset in some files.

`Parmetis`_ can be downloaded from their `downloads
<http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download>`_. Due to the
license of Parmetis we cannot distribute it with hemocell. The Parmetis
download should be copied to the  ``./hemocell/external/`` directory. If you
need it because you want load balancing to be enabled you have to extract it
with::

  cd hemocell/external && tar -xzf parmetis-4.0.3.tar.gz

.. _compilation:

Compiling HemoCell from source
------------------------------

HemoCell can be compiled from source using ``CMake``. In the root directory
``./hemocell/``, execute the following instructions to configure and compile
the HemoCell library::

  # within ./hemocell/
  mkdir build
  cd build
  cmake ..
  cmake --build .

By default only the standard ``hemocell`` library is compiled. To compile
specific examples, you should specify their corresponding compile targets.
These targets are defined for each example and match the example's directory
name. Thus, to compile the example given in ``hemocell/examples/pipeflow/pipeflow.cpp``
you would compile the target ``pipeflow``. After compilation, the executable
``pipeflow`` is placed under ``hemocell/examples/pipeflow/``. You can indicate
the desired compilation target to ``CMake`` through its ``--target`` flag::

  cmake --build . --target pipeflow

If you intend to compile multiple targets, you can repeat the previous command for
each individual target. Alternatively, if you have ``CMake`` version ``>=3.15``,
you can specify a space-separated list of all targets directly, e.g. to compile both
the ``pipeflow`` and ``parachuting`` examples::

  cmake --build . --target pipeflow parachuting

To speed up the compilation process, ``CMake`` can exploit parallelism by
providing the ``--parallel`` flag. For instance, to use all available cores on
your machine::

  cmake --build . --target pipeflow parachuting --parallel $(nproc)

Afterwards, the ``pipeflow`` and ``parachuting`` executables are placed within
the corresponding example's directories, i.e. ``examples/pipeflow`` and
``examples/parachuting``.


To test if the library is successfully compiled, you can evaluate the defined
tests::

  make test

.. _packcells:


Generating initial positions for cells
--------------------------------------

At some point you might want to run a slightly different geometry, or run your
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
the csv files into a single one per time-step you can use the script :any:`ccsv`
in the ``tmp`` directory. This will generate them for you.

The more detailed output on both the fluid field and particle field is stored in
``hdf5`` format. We recommend using the `XDMF`_ format to make these
readable for `Paraview`_ . To generate ``*.xmf`` files run the :any:`bpp`
script.

When you have created the ``*.xmf`` files you can load them into Paraview,
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
