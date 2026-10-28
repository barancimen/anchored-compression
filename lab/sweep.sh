#!/bin/sh
# Probe lab/anchor.c variants: sweep.sh 'BASE -D FLAGS' 'VARIANT FLAGS' ...  (archive bytes of image.bin)
# IMG=heldout/h2/image.bin sweep.sh ... probes another image.
cd "$(dirname "$0")/.."
mkdir -p lab/build
base=$1; shift
for v in "" "$@"; do
  cc -O2 $base $v -o lab/build/sweep$$ lab/anchor.c -lm || continue
  printf '%-40s %s\n' "${v:-(base)}" "$(./lab/build/sweep$$ "${IMG:-image.bin}" | tail -1)"
done
rm -f lab/build/sweep$$
