# Anchor-point compression of generated 4-bit images (autoresearch log)

**Question.** Can a lossless compressor that first stores a set of *anchor* pixels and then
codes the rest of the image *given* those anchors beat a plain raster-scan coder on images
from `gen_image.py`, measured by Hutter Prize rules?

**Short answer: no, not on these images.** Every anchor design tried ties a raster coder at
best in archive size, and the anchor machinery always adds program code, which counts twice.
The best anchor compressor is 313 bytes behind the no-anchor reference (S 41929 vs 41616), and
it gets there by using as few anchors as possible. The best *genuine* anchor algorithm, a
lattice of anchors refined as a quadtree, is 697 bytes behind (S 42313). Both gaps hold on
four held-out images.

## Why
1. **Coding order doesn't change the total information.** By the chain rule, coding anchors
   first and the rest given them has the same ideal total as any other order; an order wins
   only if its online model happens to learn better. Measured on the same pixels, the
   quadtree pays 2341 bytes more than raster for its anchors and coarse levels (they are
   predicted from far-away pixels) and saves 2081 bytes on its finest level (two-sided
   neighbours): net +260. Row-based refinement and anchors-plus-look-ahead land within ~1%
   of raster too.
2. **Selected anchors can't pay for their positions.** Even when anchors are handed to the
   decoder for free, the rest of the image gets cheaper by little more than those pixels'
   own cost. Choosing them by edge strength (Sobel) or by coding cost helps the *other*
   pixels by at most ~1100 bits per 4000 anchors, while sending their positions and values
   needs ~45600 bits. Edge-, corner- or cost-placed anchors miss break-even by 10x or more;
   only positions the decoder derives for free (lattices) are viable.
3. **The cost is diffuse.** The most expensive 1% of pixels carry only 5.4% of the bits; most
   of the cost is sensor noise (0.73 bits/pixel by itself) and fine texture of the random
   field, spread evenly over the image, which no anchor placement can concentrate.

## Results on `image.bin` (Hutter total S, bytes; lower is better)

S = 2 x zipped program + archive + charged command-line characters (50), C = D form.

| candidate | design | archive | zip | S | status |
|---|---|---|---|---|---|
| `000_raster_ref` | reference: earlier project's raster coder, no anchors | 39512 | 1027 | **41616** | ref |
| `001_lattice_lookahead` | A: 8x8 verbatim anchors, model reads ahead into their interpolation | 39566 | 1176 | 41968 | keep |
| `002_quadtree` | B: even/even anchor lattice, then cell centres, then edge midpoints | 40200 | 1197 | 42644 | discard |
| `003_lattice_lookahead_t128` | A with a 4x4 lattice (16 anchors) | 39529 | 1175 | **41929** | keep (best) |
| `004_quadtree_ctx` | B, transposed scan, searched counter contexts | 39939 | 1214 | 42417 | discard |
| `005_quadtree_golf` | B golfed, retuned forgetting / variance constants | 39895 | 1184 | 42313 | discard (best of B) |

Family A: anchors are side information for a raster model. Its look-ahead into the anchors'
interpolation adds nothing measurable (39565 bytes without it vs 39566 with it), and fewer
anchors always scored better (T=32 / 64 / 128: 39721 / 39566 / 39529). Its optimum is the
reference with no anchors at all. Family B is the genuine anchor algorithm.

## Held-out images (`lab/generalization.tsv`; never used for tuning)

| image | reference S | 003 S (diff) | 005 S (diff) |
|---|---|---|---|
| h1 | 41887 | 42200 (+313) | 42611 (+724) |
| h2 | 41174 | 41494 (+320) | 41902 (+728) |
| h3 | 41202 | 41515 (+313) | 41938 (+736) |
| h4 | 40261 | 40573 (+312) | 40935 (+674) |

The gaps on `image.bin` (+313, +697) hold on every held-out image, so the tuning did not
overfit. The quadtree's searched contexts and constants each saved bytes on held-out images
too (`lab/notes.md`).

## What was tried (details and every number in `lab/notes.md`)
- **Number of anchors:** lattice steps T = 1 (plain raster) to 256 for the quadtree, the
  row hierarchy and the look-ahead design; anchors verbatim ("stored in full") vs modelled
  (modelled is ~150 bytes cheaper at T=16).
- **Anchor-conditioned models:** per-level vs shared statistics (sharing counters and
  variances across levels removed most of the hierarchy's handicap), residual prediction from
  the 4 nearest anchors, two-sided directional Gaussians, per-pass mixers, stencil sizes,
  a local search over counter contexts.
- **Selection:** upper-bound tests with free anchors chosen by coding cost (an iterative
  cost-map selection), by edge strength and at random (point 2 above); scan orientation.

## Layout

| path | contents |
|---|---|
| `image.bin`, `image.png` | the benchmark (copied from the earlier project; same image, so S is comparable) |
| `gen_image.py` | the image generator, unchanged |
| `heldout/h1..h4` | four more images from unchanged copies of `gen_image.py` |
| `harness/evaluate.py` | Hutter-rule scorer (sandboxed decompression, charged command lines, `--ref`) |
| `program.md` | the loop protocol and the two design families |
| `results.tsv`, `best.json` | every scored experiment; the current best |
| `candidates/NNN_*` | `src.c` (readable), `p.c` (minified, measured), `cand.json`, `submission/P.zip` |
| `lab/anchor.c`, `lab/sweep.sh` | research coder with every anchor order as a `-D` knob |
| `lab/qt.c`, `lab/qtsweep.sh`, `lab/qtsearch.py` | quadtree template, knob sweeps, counter-context search (`lab/qtsearch.log`) |
| `lab/notes.md` | all probes with numbers, including rejected ones |
| `lab/eval_on.py`, `lab/generalization.tsv` | scoring on held-out images |

## Reproduce

```sh
export HP_SCRATCH=/some/dir/outside/this/repo     # decompression runs with the repo hidden
python3 harness/evaluate.py candidates/005_quadtree_golf          # S on image.bin
python3 lab/eval_on.py heldout/h2/image.bin candidates/005_quadtree_golf
cc -O2 -DT=4 -DRES=1 -DSHC=1 -DSHV=1 -o anchor lab/anchor.c -lm && ./anchor image.bin
```
