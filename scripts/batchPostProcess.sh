#!/bin/bash
# This is a draft script for post-processing. 

# The following does not work on MacOS X due to the lack of readlink command.
#export scriptsDir=$(readlink -f $(dirname $0))

#MacOS X compatible version:
export scriptsDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo ${scriptsDir}
[[ -z $1 ]] || ( echo "Clearing XDMF files..."; rm -rf tmp/*.xmf )
cd tmp; 


# Fluid
${scriptsDir}/FluidHDF5toXMF.py; 

# Cells of different types
${scriptsDir}/CellHDF5toXMF.py PLT RBC_HO RBC_SU SickledRBC WBC WBC_HO TumorCell; 

# vWF
${scriptsDir}/VWFHDF5toXMF.py VWF;

# Boundary particles
${scriptsDir}/ParticleField3D_HDF5toXMF.py BoundaryParticles; 

# trombosit
${scriptsDir}/ParticleField3D_HDF5toXMF.py ActivatedBoundaryParticles;
${scriptsDir}/../trombosit/BondParticleField3D_HDF5toXMF.py BondFieldParticles; 
