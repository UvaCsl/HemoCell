#!/bin/bash
if [ ! -d ./csv ]; then
  echo "No csv folder, exiting"
  exit 0
fi
for dir in ./csv/*/
do
  timestep=`echo ${dir} | cut -d/ -f3`
  filenames=`ls ${dir}* | cut -d/ -f4 | cut -d. -f1 | sort| uniq`

  for name in $filenames
  do 
    echo "Merging ${dir}, Filename ${name}, into ${name}.${timestep}.csv"
    head -1 `ls ${dir}${name}* | head -1` > ${name}.${timestep}.csv
    for file in ${dir}${name}*
    do
      cat $file | sed -n '1!p' >> ${name}.${timestep}.csv
    done
  done 
done
