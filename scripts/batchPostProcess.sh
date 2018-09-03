#!/bin/bash
# This is a draft script for post-processing. 

# The following does not work on MacOS X due to the lack of readlink command.
#export scriptsDir=$(readlink -f $(dirname $0))

#MacOS X compatible version:
export scriptsDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo ${scriptsDir}
[[ -z $1 ]] || ( echo "Clearing XDMF files..."; rm -rf tmp/*.xmf )
if [ -d tmp ]; then
  cd tmp; 
fi

if command -v python; then
python_c=python
elif commond -v python3; then
python_c=python3
else
pyton_c=python2
fi

# Fluid
${python_c} ${scriptsDir}/FluidHDF5.py; 

# Cells of different types
if [ ! -d ./csv ]; then
  echo "No hdf5 folder, exiting"
  exit 0
fi
for dir in ./hdf5/*/
do
  filenames=`ls ${dir}* | cut -d/ -f4 | cut -d. -f1 | sort| uniq`
  for name in $filenames
  do
    if [ $name == "Fluid" ]; then
      continue
    fi
    if [ $name == "Fluid_PRE" ]; then
      ${python_c} ${scriptsDir}/FluidHDF5.py Fluid_PRE; 
      continue
    fi
    if [ $name == "CEPAC" ]; then
      ${python_c} ${scriptsDir}/FluidHDF5.py CEPAC; 
      continue
    fi
    echo ${name}:
    ${python_c} ${scriptsDir}/CellHDF5toXMF.py ${name}; 
  done
  break;
done
# vWF
#${scriptsDir}/VWFHDF5toXMF.py VWF;

# Boundary particles
#${scriptsDir}/ParticleField3D_HDF5toXMF.py BoundaryParticles; 

# trombosit
#${scriptsDir}/ParticleField3D_HDF5toXMF.py ActivatedBoundaryParticles;
#${scriptsDir}/../trombosit/BondParticleField3D_HDF5toXMF.py BondFieldParticles; 
