#!/bin/sh
# Probe lab/qt.c (quadtree template) variants: qtsweep.sh 'FLAGS1' 'FLAGS2' ...  (archive bytes)
cd "$(dirname "$0")/.."
mkdir -p lab/build
for v in "" "$@"; do
  cc -O2 $v -o lab/build/qtsw$$ lab/qt.c -lm 2>/dev/null || { echo "$v: build failed"; continue; }
  printf '%-36s %s\n' "${v:-(base)}" "$(./lab/build/qtsw$$ < "${IMG:-image.bin}" | wc -c)"
done
rm -f lab/build/qtsw$$
