##########################################################################
## Makefile for the Palabos example program cavity3d.
##
## The present Makefile is a pure configuration file, in which 
## you can select compilation options. Compilation dependencies
## are managed automatically through the Python library SConstruct.
##
## If you don't have Python, or if compilation doesn't work for other
## reasons, consult the Palabos user's guide for instructions on manual
## compilation.
##########################################################################

# USE: multiple arguments are separated by spaces.
#   For example: projectFiles = file1.cpp file2.cpp
#                optimFlags   = -O -finline-functions

# Leading directory of the Palabos source code
# palabosRoot   = ./palabos/	
#palabosRoot   = /home/lmount/devel/palabos/
#palabosRoot   = $(PALABOS_ROOT)
palabosRoot   = ./palabos
ficsionRoot   = .
# Name of source files in current directory to compile and link with Palabos
# projectFiles = losglobulos.cpp
projectFiles = ficsional.cpp
 
# Choose between generic and precompiled mode
precompiled  = false
# Set optimization flags on/off
optimize     = false
# Set debug mode and debug flags on/off
debug        = true
# Set profiling flags on/off
#profile      = false
profile      = true
# Set MPI-parallel mode on/off (parallelism in cluster-like environment)
MPIparallel  = true
# Set SMP-parallel mode on/off (shared-memory parallelism)
SMPparallel  = false
# Decide whether to include calls to the POSIX API. On non-POSIX systems,
#   including Windows, this flag must be false, unless a POSIX environment is
#   emulated (such as with Cygwin).
usePOSIX     = true
# This flag must be defined true if you are using the external library
# cvmlcpp. But first, you need to download cvmlcpp and put it into the
# directory "externalLibraries" of the Palabos root.
useCVMLCPP   = false
# This flag must be defined true if you want to use the parallel I/O 
# library HDF5 for the output of the files. HDF5 library must be 
# installed.
useHDF5 = true

# Path to external libraries (other than Palabos)
libraryPaths = 
# Path to inlude directories (other than Palabos)
includePaths = $(ficsionRoot) $(ficsionRoot)/core/ $(ficsionRoot)/helper/ $(ficsionRoot)/models/ $(ficsionRoot)/cases/ $(ficsionRoot)/external/ $(ficsionRoot)/IO/ $(ficsionRoot)/newMesh/
# Path to inlude directories (other than Palabos)
# ficsionPaths = $(ficsionRoot)/core/ $(ficsionRoot)/helper/ $(ficsionRoot)/models/ $(ficsionRoot)/cases/ $(ficsionRoot)/external/
# Dynamic and static libraries (other than Palabos)
libraries    = 

# Compiler to use without MPI parallelism
serialCXX    = g++
UNAME := $(shell uname)
HOSTNAME := $(shell hostname)
ifeq ($(UNAME),Darwin)
# Compiler to use with MPI parallelism
	parallelCXX  = /usr/local/bin/mpicxx
# General compiler flags (e.g. -Wall to turn on all warnings on g++)
	compileFlags = -DPLB_MAC_OS_X -Wall -fpermissive
else
# Compiler to use with MPI parallelism
    ifeq ($(HOSTNAME), kolmogorov) # WorkMachine needs special treatment!
        parallelCXX  = mpicxx.openmpi
    else
        parallelCXX = mpicxx
    endif
# General compiler flags (e.g. -Wall to turn on all warnings on g++)
    compileFlags = -Wall 
endif

ifeq ($(useHDF5), true)
    compileFlags  += -lhdf5_cpp -lhdf5 -lhdf5_hl_cpp -lhdf5_hl -DFCN_USE_HDF5
    libraries     += -lhdf5 -lhdf5_hl
endif

compileFlags += -Wreorder -Wsign-compare
# General linker flags (don't put library includes into this flag)
linkFlags    = -Wreorder -Wsign-compare
# Compiler flags to use when optimization mode is on
optimFlags   = -O3 -finline-functions
# Compiler flags to use when debug mode is on
debugFlags   = -g
# Compiler flags to use when profile mode is on
profileFlags = -pg


##########################################################################
# All code below this line is just about forwarding the options
# to SConstruct. It is recommended not to modify anything there.
##########################################################################

SCons     = $(palabosRoot)/scons/scons.py -j 4 -f $(ficsionRoot)/SConstruct
SConsArgs = palabosRoot=$(palabosRoot) \
			ficsionRoot=$(ficsionRoot) \
            projectFiles="$(projectFiles)" \
            srcPaths="$(srcPaths)" \
            ficsionPaths="$(ficsionPaths)" \
            precompiled=$(precompiled) \
            optimize=$(optimize) \
            debug=$(debug) \
            profile=$(profile) \
            MPIparallel=$(MPIparallel) \
            SMPparallel=$(SMPparallel) \
            usePOSIX=$(usePOSIX) \
	    useCVMLCPP=$(useCVMLCPP) \
            serialCXX=$(serialCXX) \
            parallelCXX=$(parallelCXX) \
            compileFlags="$(compileFlags)" \
            linkFlags="$(linkFlags)" \
            optimFlags="$(optimFlags)" \
            debugFlags="$(debugFlags)" \
	    profileFlags="$(profileFlags)" \
	    libraryPaths="$(libraryPaths)" \
	    includePaths="$(includePaths)" \
	    libraries="$(libraries)"

all:
	python $(SCons) $(SConsArgs)

compile:
	python $(SCons) $(SConsArgs)

clean:
	python $(SCons) -c $(SConsArgs)
	/bin/rm -vf `find $(palabosRoot) -name '*~'`
