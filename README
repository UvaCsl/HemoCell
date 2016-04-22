ficsion
==========

`ficsion` is a framework for simulating suspensions of deformable cells, focusing on blood. It is based of the combined Immersed boundary-lattice Boltzmann method (IB-LBM) and is built on top of the open source `C++` lattice Boltzmann solver `palabos`.

`ficsion` has been developed as a framework and not as a single monolithic application. The user/developer can choose the geometry and type of flow he/she will use, the type of cells (red blood cells, platelets, etc.), the model which will be applied to them and a variety of actions: from applying forces to a certain number of surface particles to letting platelets stick to the walls or between each other. In that sense, the already implemented files serve as examples where the user can take the modules and built his/her own application, like lego-bricks.

Building Prerequisites
====================

`ficsion` is built on top of `palabos` and it should be possible to be build on any GNU/Linux system with at least the GNU GCC compiler suite and OpenMPI. Make and SCONS or CMake are used as the build automation tools. GCC versions 4.7.2 and 4.8.1 as well as OpenMPI version 1.4.5 and 1.6.5 have been successfully tested. It is very likely that other compilers will also produce favorable results. `ficsion` has been successfully built and ran on x64 and IBM BG/Q systems. 

## `palabos`
The working version of `palabos` is [v1.5r1](http://www.palabos.org/images/palabos_releases/palabos-v1.5r1.zip). Earlier versions of palabos are no longer supported duo to strong dependency on sparse particle classes.

## HDF5
`ficsion` uses the `hdf5` library for the output and post-processing of the results. Debian packeges `h5utils hdf5-tools libhdf5-serial-dev` are known to work for the I/O.

## Post-processing
Almost all the post-processing is performed by python scripts. The necessary libraries are `numpy`, `matplotlib` and `h5py`.

## Debian packages

```bash
sudo apt-get install gcc g++ gdb make scons cmake
sudo apt-get install openmpi-bin openmpi-checkpoint openmpi-common libopenmpi-dev
sudo apt-get install h5utils hdf5-tools libhdf5-serial-dev
sudo apt-get install python python-numpy python-matplotlib python-h5py
```


## Known issues
* `ficsion` compiles just fine on MacOSX v10.11 (El Capitain), yet it break when a case is initiated.
* Ordinary `Makefile`s do not work for IBM BG/Q systems. Instead use the `Makefile.cineca` as a template, found in the parent directory.

Building Instructions
====================

Ensure that the above packages are downloaded and installed on your system. Assuming that ficsion lies in the directory `${FICSION}`,  e.g. 
```shell
export FICSION=${HOME}/ficsion
```
download `palabos` v1.5r1 from its [website](http://www.palabos.org/images/palabos_releases/palabos-v1.5r1.zip) and extract it to the path `${FICSION}/palabos`. Now navigate into one of the `case` subfolders and build the corresponding case:

``` shell
make
```
and it will produce an executable with the name of the case. Should you want to change debugging or profiling options, please alter the `Makefile`. Follow a similar approach to compile the examples.

Alternative building instructions for CMake (preferred)
=======================================================

To build a case, issue the following from its folder:

```shell
mkdir build
cd build
cmake ..
make -j 4
```

The cmake tries to make use of the system libraries instead of building them, thus you might also want to install
development packages for Eigen and TinyXML (Palabos 1.5+).
The cmake files are also suitable to use in IDE-s supporting cmake projects. 
Both ficsion and palabos is handled as an external library relative to the case in the cmake project structure.


Documentation
=============

Program Documentation can be found in `${FICSION}/doc/ficsion_UserGuide.pdf` which provides detailed information on the `ficsion` configuration.
