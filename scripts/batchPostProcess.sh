#!/bin/bash
# Write the XDMF files (*.xmf) that describe HemoCell HDF5 output for ParaView.
#
# HemoCell writes them itself while it runs; this is for output of older
# versions or to rewrite them. Run it in a case directory or an output
# directory, or pass directories:
#
#   batchPostProcess.sh [--force] [directory ...]
#
# See HDF5toXDMF.py for the details. Requires Python 3 with h5py.
scriptsDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
exec python3 "${scriptsDir}/HDF5toXDMF.py" "$@"
