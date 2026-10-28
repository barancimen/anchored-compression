# Inpainting from anchors

`image.bin` is reduced to its **anchors**, the lattice of step 8 (x % 8 == 0 and y % 8 == 0:
63 x 63 = 3,969 pixels, 1.6% of the image; `anchors.png`). The quadtree anchor model then fills
in everything else, and `animation.gif` steps through 100 independent fills: one bare 500x500
frame per sample, 25 fps, looping (a 4 s loop).

## How it works (`inpaint.c`)
`inpaint.c` is the quadtree path of `../lab/anchor.c` (the research anchor coder) with the
coder removed; the changed lines are marked `INPAINT`. Settings: the best lab quadtree knobs
(shared counters and variances, residual RLS, directional inputs, radius-6 stencils).
1. **Training.** The model codes the 8 generator images in `../lab/train` completely, learning
   as the compressor does, then the target's anchors. It never sees a hidden pixel.
2. **Filling.** The model goes coarse to fine as the compressor would (for t = 8, 4, 2: cell
   centres, then edge midpoints), but draws every bit from its own probability instead of
   reading it. It is frozen while drawing: counters, mixer, RLS weights and variances keep
   their trained values; only the error buckets follow the drawn pixels.
3. **One sample per `fork()`**, seeded `SEED << 20 | n` (splitmix64).

Port check: with no training images, the model's cost of the true image is 40,299 bytes,
identical to `../lab/anchor.c` with the same knobs.

## Results (seed 1)
| | bits per filled pixel |
|---|---|
| the 100 samples, model's own cost of what it drew | 1.241 - 1.254 (mean 1.247) |
| the true hidden pixels under the same frozen model | 1.297 |
| the true pixels coded with learning on (the compressor, after training) | 1.267 |

- Every sample keeps all 3,969 anchors exactly (checked by `gif.py`).
- The samples keep the layout the anchors pin down: the dark rectangle, the dark ellipse,
  the bright and dark cloud regions.
- Anything thinner than the 8-pixel spacing is redrawn differently every time. Thin lines,
  rings and stripes mostly vanish or become blotches; object outlines turn ragged. The model
  only knows local texture and edges, not shapes.
- The samples are a little more regular than the real image (lower own cost than the true
  pixels' 1.297).

## Files
| path | contents |
|---|---|
| `inpaint.c` | the sampler (usage in its header comment) |
| `run.sh` | builds, draws the samples, makes the pictures: `./run.sh [NSAMPLES=100] [SEED=1] [IMAGE=../image.bin]` |
| `gif.py` | `animation.gif` and `anchors.png`; checks the anchors of every sample |
| `animation.gif` | the samples at 25 fps, looping, no borders or labels |
| `anchors.png` | the input: the anchor pixels on black |
| `samples/stats.tsv` | per sample: index, seed, bits per drawn pixel; then the two `true` lines |

`samples/sample_NNN.bin` and the `inpaint` binary are not tracked; `./run.sh` regenerates
them deterministically (about 3 minutes on one core).
