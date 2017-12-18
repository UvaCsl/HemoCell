#!/bin/sh
trap "exit" INT

if [ -d "palabos" ]; then
echo "=================== Removing existing palabos source ========================"
  rm -r ./palabos
fi

echo "=================== Uncompressing palabos ==========================="
tar -xvzf palabos_dev.tgz
mv palabos-* palabos

echo "=================== Patching palabos ============================="
cd patch
./patchPLB.sh
cd ..

echo "================= Done ===================="
