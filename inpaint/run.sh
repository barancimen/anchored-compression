#!/bin/sh
# Inpaint an image from its anchors (lattice of step 8) with the quadtree anchor model, then
# make the animation. The model is trained on the 8 generator images in ../lab/train first.
# Usage: ./run.sh [NSAMPLES=100] [SEED=1] [IMAGE=../image.bin]
set -e
cd "$(dirname "$0")"
N=${1:-100}
SEED=${2:-1}
IMG=${3:-../image.bin}
cc -O2 -o inpaint inpaint.c -lm
mkdir -p samples
rm -f samples/sample_*.bin
printf 'sample\tseed\tbpp_drawn\n' > samples/stats.tsv
./inpaint "$IMG" "$N" "$SEED" samples ../lab/train/t0*/image.bin >> samples/stats.tsv
python3 gif.py "$IMG" samples
