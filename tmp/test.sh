#!/bin/bash 
mpicxx.openmpi -o test.o -c -Wall -lhdf5_cpp -lhdf5 -lhdf5_hl_cpp -lhdf5_hl -DFCN_USE_HDF5 -O3 -finline-functions -g -DPLB_DEBUG -DPLB_MPI_PARALLEL -DPLB_USE_POSIX -Ipalabos/src -Ipalabos/externalLibraries -I. -Icore -Ihelper -Imodels -Icases -Iexternal -IIO $1


