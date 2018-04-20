#!/bin/bash
set -e

cd ${CI_PROJECT_DIR}/cases/pipeflow
mpirun --allow-run-as-root -n 4 ./pipeflow ../../scripts/ci/config-pipeflow.xml
cd log
if [ -n "`cat logfile | grep "# of cells" | cut -d: -f2 | cut -d" " -f2 | grep -v 32`" ]; then
  echo "Error in number of cells, 32 Expected"
  exit 1
fi
if [ -n "`cat logfile |grep "viscosity" | cut -d: -f 4 |cut -d" " -f2 | xargs -n 1 echo 1.5 \< | bc | grep -v 1`" ]; then
  echo "Error, viscosity goes below 1.5"
  exit 1
fi
if [ -n "`cat logfile |grep "viscosity" | cut -d: -f 4 |cut -d" " -f2 | xargs -n 1 echo 5.0 \> | bc | grep -v 1`" ]; then
  echo "Error, viscosity goes above 5.0"
  exit 1
fi
if [ -n "`cat logfile |grep "Force  -" | cut -d: -f 3 |cut -d" " -f2 | xargs -n 1 echo 6.0 \> | sed -e 's/[eE]+*/\*10\^/' | bc | grep -v 1`" ]; then
  echo "Error, Max force goes above 6.0"
  exit 1
fi
