#!/bin/sh
trap "exit" INT

if [ -d "palabos" ]; then
echo "=================== Removing existing palabos source ========================"
  rm -r ./palabos
fi

if [ ! -e "palabos_dev.tgz" ]; then
echo "=================== Downloading palabos =========================="
  wget -O palabos_dev.tgz http://www.palabos.org/images/palabos_releases/palabos-v2.0r0.tgz || { echo "Error Downloading palabos, exiting ..."; exit 1;}
fi

echo "=================== Uncompressing palabos ==========================="
tar -xvzf palabos_dev.tgz
mv palabos-* palabos

echo "=================== Patching palabos ============================="
cd patch
./patchPLB.sh
cd ..

echo "================= Done ===================="
