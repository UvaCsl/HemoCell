#!/bin/bash
# This is a draft script for post-processing. 

# The following does not work on MacOS X due to the lack of readlink command.
export scriptsDir=$(readlink -f $(dirname $0))
echo ${scriptsDir}
[[ -z $1 ]] || ( echo "Clearing XDMF files..."; rm -rf tmp/*.xmf )
cd tmp; 

# Fluid
${scriptsDir}/HDF5toXMF.py; 
# Cells
${scriptsDir}/CellHDF5toXMF.py PLT RBC SickledRBC TumorCell ; 

# Boundary particles
${scriptsDir}/ParticleField3D_HDF5toXMF.py BoundaryParticles; 

# trombosit
${scriptsDir}/ParticleField3D_HDF5toXMF.py ActivatedBoundaryParticles;
${scriptsDir}/../trombosit/BondParticleField3D_HDF5toXMF.py BondFieldParticles; 
