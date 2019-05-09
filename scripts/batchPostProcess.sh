#!/bin/bash
# This is a draft script for post-processing. 

# The following does not work on MacOS X due to the lack of readlink command.
#export scriptsDir=$(readlink -f $(dirname $0))

#MacOS X compatible version:
export scriptsDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo ${scriptsDir}

if command -v python; then
python_c=python
elif command -v python3; then
python_c=python3
else
python_c=python2
fi

erase_things=$1

function do_work {
  if [ ! -d ./hdf5 ]; then
    return
  fi

  [[ -z $erase_things ]] || ( echo "Clearing XDMF files..."; rm -rf *.xmf )

  # Fluid
  ${python_c} ${scriptsDir}/FluidHDF5.py; 

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
}

do_work;

for dir in */ ; do
  cd $dir;
  do_work;
  cd ../;
done
