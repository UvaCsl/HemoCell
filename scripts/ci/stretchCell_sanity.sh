#!/bin/bash
set -e

cd ${CI_PROJECT_DIR}/examples/stretchCell
mpirun --allow-run-as-root -n 4 ./stretchCell ../../scripts/ci/config-stretchCell.xml
cd tmp/log
if [ -n "`cat logfile | grep "diameter" | cut -d: -f2 |cut -d" " -f2 | xargs -n 1 echo 9.6 \> | bc | grep -v 1`" ]; then
  echo "Error stretch is larger then expected (>9.6)"
  exit 1
fi
if [ -n "`cat logfile | grep "Volume:" | cut -d: -f3 |cut -d"(" -f2 | cut -d% -f1 | xargs -n 1 echo 100.1 \> | bc |grep -v 1`" ]; then
  echo "Error, Volume increases above 100.1%"
  exit 1
fi
if [ -n "`cat logfile | grep "Volume:" | cut -d: -f3 |cut -d"(" -f2 | cut -d% -f1 | xargs -n 1 echo 100 \< | bc |grep -v 1`" ]; then
  echo "Error, Volume decreases below 100%"
  exit 1
fi
if [ -n "`cat logfile | grep "Volume:" | cut -d: -f3 |cut -d" " -f2 | xargs -n 1 echo 81.19 \> | bc | grep -v 1`" ]; then
  echo "Error, Volume increases above 81.19µm³"
  exit 1
fi
if [ -n "`cat logfile | grep "Volume:" | cut -d: -f3 |cut -d" " -f2 | xargs -n 1 echo 81.12 \< | bc | grep -v 1`" ]; then
  echo "Error, Volume decreases below 81.12µm³"
  exit 1
fi
if [ -n "`cat logfile | grep "Surface:" | cut -d: -f2 |cut -d" " -f2 | xargs -n 1 echo 133.04 \> | bc | grep -v 1`" ]; then
  echo "Error, Surface above 133.04µm²"
  exit 1
fi
if [ -n "`cat logfile | grep "Surface:" | cut -d: -f2 |cut -d" " -f2 | xargs -n 1 echo 129.34 \< | bc | grep -v 1`" ]; then
  echo "Error, Surface below 129.34µm²"
  exit 1
fi
