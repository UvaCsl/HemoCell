#!/bin/sh
trap "exit" INT

# This script performs the setup of Palabos: it retrieves a tagged release from
# their Gitlab pages, afterwhich additional hemocell-specific features are
# added by applying the patch at `hemocell/patch`.

# supported tag and download target
tag="v2.3.0"
target="https://gitlab.com/unigespc/palabos/-/archive/${tag}/palabos-${tag}.tar.gz"
archive="palabos.tar.gz"

# clean old palabos
if [ -d "palabos" ]; then
  rm -r ./palabos
fi

# obtain new tag
if [ ! -e "${archive}" ]; then
  wget -qO "${archive}" "${target}" || { echo "Error Downloading palabos, exiting ..."; exit 1;}
fi

# extract source
tar -xzf "${archive}" && rm "${archive}"
mv palabos-* palabos

# apply the patch
(
  cd patch || { echo "Patch directory not present"; exit 1; }
  ./patchPLB.sh
)
