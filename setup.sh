#!/bin/sh
trap "exit" INT

# supported tag and download target
tag="v2.2.1"
target="https://gitlab.com/unigespc/palabos/-/archive/${tag}/palabos-${tag}.tar.gz"

# clean old palabos
if [ -d "palabos" ]; then
  rm -r ./palabos
fi

# obtain new tag
if [ ! -e "palabos.tar.gz" ]; then
  wget -O palabos.tar.gz ${target} || { echo "Error Downloading palabos, exiting ..."; exit 1;}
fi

# extract source
tar -xzf palabos.tar.gz
mv palabos-* palabos

# apply the patch
(
  cd patch || { echo "Patch directory not present"; exit 1; }
  ./patchPLB.sh
)
