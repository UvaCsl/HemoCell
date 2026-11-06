
HemoCell Getting Started
========================

.. _from_source:

Setting up HemoCell from source
-------------------------------

The following packages are requirements for compiling and/or running HemoCell from source.

==========================         ==========================
Required package                   Tested versions
==========================         ==========================
`OpenMpi`_                         4.1.6 
`GCC`_                             13.3.0
`CMake`_                           3.28.3
`HDF5`_                            1.10.10
`GNU Patch`_                       2.7.6
`h5py`_                            2.6.0
`Palabos`_                         2.3.0
`Parmetis`_ (optional)             4.0.3
==========================         ==========================

.. note::
  The currently tested versions are listed. Older/newer versions might work as well.
  Avoid OpenMPI 2.0.X as in our experience it introduces memory leaks.

On Ubuntu 24.04 these dependencies can be installed by running::

  sudo apt-get install -y \
        make \
        cmake \
        g++ \
        libopenmpi-dev \
        libhdf5-dev \
        patch \
        python3-h5py

The main external dependency `Palabos`_ has to be added and patched.
We currently support the ``v2.3.0`` version with an additional small patch. This can be done in two ways. 

1. To automatically download this specific Palabos version and apply our patch, run the setup script
``hemocell/setup.sh``::

  # from the hemocell/ root directory execute
  ./setup.sh

At this point HemoCell should be ready for compilation and development, you can move on to the next section below.

2. Alternatively, you can opt to manually download Palabos through their `releases
<https://gitlab.com/unigespc/palabos/-/releases>`_. This is useful if you want to use a different version. After downloading, Palabos
should be extracted to ``./hemocell/palabos``::

  tar -xzf palabos-v2.3.0.tar.gz
  mv palabos-v2.3.0 ./hemocell/palabos

After this Palabos must be patched, see :ref:`patching-palabos`. This can be
done by running ``./patchPLB.sh`` from the ``./hemocell/patch/`` directory, like
so::

  cd hemocell/patch && ./patchPLB.sh

The patching should succeed even though there might be an offset in some files. 
At this point HemoCell should be ready for compilation and development.

(Optional) `Parmetis`_ is used for upcoming load-balancing routines (in development atm.). 
It can be installed via adding ``libparmetis-dev`` on Ubuntu, or it can be downloaded from the
`ParMETIS repository <https://github.com/KarypisLab/ParMETIS>`_. Due to the
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

By default only the standard ``hemocell`` library is compiled, as a ``Release``
build optimised for the CPU of the build machine (``-march=native``). Use
``-DHEMOCELL_NATIVE_ARCH=OFF`` for binaries that must run on other CPUs (e.g.
other nodes of a cluster), and ``-DCMAKE_BUILD_TYPE=Debug`` for a debug build.
This operation can take a while, we suggest a coffee.

With ``CMake >= 3.21`` the predefined presets are a convenient alternative; each
builds into its own ``build-<preset>`` directory::

  cmake --list-presets              # release, portable, debug, sanitize
  cmake --preset release
  cmake --build --preset release

Machine-specific settings (for example ``-DHDF5_ROOT=...``) can be kept in a
``CMakeUserPresets.json`` that inherits from these presets.

To compile specific examples (that link to the library), you should specify their corresponding compile targets.
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
``examples/parachuting``. The targets ``hemocell_examples`` and
``hemocell_cases`` build all examples and all research cases.


To test if the library is successfully compiled, build the test runner and run
the tests (``unit`` takes seconds, ``validation`` a few minutes; both also run
on two MPI ranks)::

  cmake --build . --target hemocell_test
  ctest -L unit --output-on-failure
  ctest -L validation --output-on-failure

Generating initial positions for cells
--------------------------------------

The initial positions of the cells (the ``<type>.pos`` files) for another
geometry or hematocrit are made with the ``packCells`` tool in
``tools/packCells``. It has its own CMake project; from the repository root::

  cmake -S tools/packCells -B tools/packCells/build
  cmake --build tools/packCells/build
  tools/packCells/packCells --help

The binary is written to ``tools/packCells/packCells``. Sizes are in µm, not
lattice units. See :ref:`packcells` for how to pack for a case, its walls
and its hematocrit, and how to check the result.


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

Post-processing the output of a HemoCell case
---------------------------------------------

A HemoCell case produces multiple types of output by default. The simplest is the ``csv``
output which consists of all the information about cells in csv files, one per
cell type and time step.

The more detailed output on both the fluid field and particle field is stored in compressed
``hdf5`` format. HemoCell describes these files in the `XDMF`_ format, which
visualization tools such as `Paraview`_ read: next to the ``hdf5`` directory
the output directory holds ``Fluid.xmf``, ``RBC.xmf``, etc. with all output
steps, and one ``.xmf`` file per step. Open e.g. ``tmp/RBC.xmf`` in Paraview
(any of its XDMF readers) to get the time series. See
:ref:`reference/outputs:Files` for the layout.

For output of older HemoCell versions, generate the ``*.xmf`` files with the
:any:`bpp` script, e.g. from the folder of a case::

  ../../scripts/batchPostProcess.sh

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
.. _XDMF: https://www.xdmf.org/
.. _GNU Patch: https://savannah.gnu.org/projects/patch/
.. _IntelMPI: https://software.intel.com/content/www/us/en/develop/tools/mpi-library.html
.. _OpenMPI: https://www.open-mpi.org/
.. _GCC: https://gcc.gnu.org/
.. _CMake: https://cmake.org/
.. _HDF5: https://www.hdfgroup.org/
.. _h5pY: https://www.h5py.org/
.. _Parmetis: https://github.com/KarypisLab/ParMETIS
.. _Palabos: https://palabos.unige.ch/
