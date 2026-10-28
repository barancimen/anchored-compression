# Anchor-point compression of generated 4-bit images (autoresearch log)

**Question.** Can a lossless compressor that first stores a set of *anchor* pixels and then
codes the rest of the image *given* those anchors beat a plain raster-scan coder on images
from `gen_image.py`, measured by Hutter Prize rules?

**Short answer: no, not on these images.** Every anchor design tried ties a raster coder at
best in archive size, and the anchor machinery always adds program code, which counts twice.
The best anchor compressor is 313 bytes behind the no-anchor reference (S 41929 vs 41616), and
it gets there by using as few anchors as possible. The best *genuine* anchor algorithm, a
lattice of anchors refined as a quadtree, is 697 bytes behind (S 42313). Both gaps hold on
four held-out images. A second idea, training the model on images simulated from the same
recipe inside the decompressor, saves 119 archive bytes but costs 1809 bytes of S.

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

## Score table: `image.bin` (Hutter total S, bytes; lower is better)

S = 2 x zipped program + archive + charged command-line characters (34-58), C = D form.
Every row is a verified round trip in this repo's harness (`results.tsv`).
Kinds: **anchor** = this project's anchor candidates; **ref** = the earlier raster project's
algorithms, copied here (`000` is its best, `r00..r16` the rest); **sim** = source simulation.

| S | candidate | kind | archive | zip | what it is |
|---|---|---|---|---|---|
| **41616** | `000_raster_ref` | ref | 39512 | 1027 | raster context mixing, 112-neighbour RLS (earlier project's v15) |
| 41638 | `r16_cm_v14` | ref | 39534 | 1027 | v14: context re-searched under the transposed scan |
| 41652 | `r15_cm_v13` | ref | 39548 | 1027 | v13: transposed + mirrored scan |
| 41821 | `r14_cm_v12` | ref | 39737 | 1017 | v12: micro-golf |
| 41847 | `r13_cm_v11` | ref | 39737 | 1030 | v11: 5th counter context |
| 41919 | `r12_cm_v10` | ref | 39835 | 1017 | v10: golf |
| 41923 | `r11_cm_v9` | ref | 39819 | 1027 | v9: golf |
| **41929** | `003_lattice_lookahead_t128` | anchor | 39529 | 1175 | best anchor compressor: 16 verbatim anchors + look-ahead |
| 41950 | `r10_cm_v8` | ref | 39810 | 1045 | v8: searched counter contexts |
| 41968 | `001_lattice_lookahead` | anchor | 39566 | 1176 | 64 verbatim anchors + look-ahead |
| 42197 | `r09_cm_v7` | ref | 40031 | 1058 | v7: golf step (first scored here) |
| 42250 | `r08_cm_v6` | ref | 40038 | 1081 | v6: RLS radius 6 (84 inputs) |
| **42313** | `005_quadtree_golf` | anchor | 39895 | 1184 | best quadtree: anchor lattice + refinement |
| 42417 | `004_quadtree_ctx` | anchor | 39939 | 1214 | quadtree, transposed, searched contexts |
| 42512 | `r07_cm_v5` | ref | 40204 | 1125 | v5: golfed v4 |
| 42644 | `002_quadtree` | anchor | 40200 | 1197 | first quadtree |
| 42854 | `r06_cm_v4` | ref | 40202 | 1297 | v4: RLS radius 4 (41 inputs) |
| 43425 | `006_simprior` | sim | 39393 | 1987 | `000` trained on 3 simulated images before coding |
| 43803 | `r05_cm_v3` | ref | 40451 | 1647 | v3: 4 contexts + Gaussian/OLS input |
| 44334 | `r03_cm_v1` | ref | 42016 | 1130 | v1: 8 neighbour contexts, no linear predictor |
| 45786 | `r04_cm_ols` | ref | 41628 | 2050 | v2: v1 + OLS predictor (first scored here) |
| 63488 | `r01_py_lzma` | ref | 62870 | 289 | Python standard-library LZMA |
| 125336 | `r00_store` | ref | 125000 | 151 | stored verbatim |

`r04_cm_ols` is copied unchanged except for one removed line that dumped a cost map to a
scratch path (instrumentation). The earlier project's sandbox check (a "compressor" that
reads the original from an absolute path) is not an algorithm and was not copied.

## Held-out images (`lab/generalization.tsv`; never used for tuning)

| image | `000` reference | `003` best anchor | `005` best quadtree | `006` simulation |
|---|---|---|---|---|
| h1 | 41887 | 42200 (+313) | 42611 (+724) | 43748 (+1861) |
| h2 | 41174 | 41494 (+320) | 41902 (+728) | 43017 (+1843) |
| h3 | 41202 | 41515 (+313) | 41938 (+736) | 43109 (+1907) |
| h4 | 40261 | 40573 (+312) | 40935 (+674) | 42113 (+1852) |

The gaps on `image.bin` (+313, +697, +1809) hold on every held-out image.

## Simulating the source (`006_simprior`)
Idea: let the decompressor generate images from the same recipe and train the model on them
before decoding, so the model's learning is paid in CPU time instead of archive bytes.
`lab/sim.c` is a C port of `gen_image.py`: three random fields made by FFT, the light
gradient, the six object kinds, percentile scaling, sensor noise. Its images look like the
real ones and work as well as real ones for training.
- **The ceiling is low.** Training on 1 to 8 *real* generator images (a perfect simulator)
  saves at most ~90 bytes on `image.bin`. The model adapts within a few thousand pixels
  (least-squares memory ~3300 pixels, counters ~80 updates), so earlier images only give it
  a warm start; settings that remember longer get worse, not better.
- **The simulator is expensive.** It adds 960 zipped bytes of program, counted twice.
- **Result:** archive 39393 (-119), S 43425 (+1809 vs the reference); held-out archives
  -21 to -85.

## Inpainting from anchors (`inpaint/`)
The quadtree model used generatively: `image.bin` is reduced to its step-8 anchor lattice
(3,969 pixels, 1.6%) and the model, trained only on the 8 images in `lab/train` and on those
anchors, fills in the rest by drawing every bit from its own probabilities.
`inpaint/animation.gif` loops through 100 samples at 25 fps (bare 500x500 frames). The samples
keep the large-scale layout the anchors fix and redraw everything finer than 8 pixels; thin
lines, rings and stripes mostly vanish. Details: `inpaint/README.md`.

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
- **Source simulation:** priming with real and simulated images, with faster and slower
  learning settings (section above).

## Layout

| path | contents |
|---|---|
| `image.bin`, `image.png` | the benchmark (copied from the earlier project; same image, so S is comparable) |
| `gen_image.py` | the image generator, unchanged |
| `heldout/h1..h4` | four more images from unchanged copies of `gen_image.py` |
| `harness/evaluate.py` | Hutter-rule scorer (sandboxed decompression, charged command lines, `--ref`) |
| `program.md` | the loop protocol, the design families and the reference points |
| `results.tsv`, `best.json` | every scored experiment; the best anchor candidate |
| `candidates/NNN_*` | this project's candidates: `src.c` (readable), `p.c` (minified, measured), `cand.json`, `submission/P.zip` |
| `candidates/r00..r16` | the earlier raster project's algorithms |
| `lab/anchor.c`, `lab/sweep.sh` | research coder with every anchor order as a `-D` knob |
| `lab/qt.c`, `lab/qtsweep.sh`, `lab/qtsearch.py` | quadtree template, knob sweeps, counter-context search (`lab/qtsearch.log`) |
| `lab/prime.c`, `lab/sim.c`, `lab/train/` | priming probe, C port of the generator, 8 generator images for training |
| `lab/notes.md` | all probes with numbers, including rejected ones |
| `lab/eval_on.py`, `lab/generalization.tsv` | scoring on held-out images |
| `inpaint/` | inpainting from the anchors with the quadtree model; `animation.gif` |

## Reproduce

```sh
export HP_SCRATCH=/some/dir/outside/this/repo     # decompression runs with the repo hidden
python3 harness/evaluate.py candidates/005_quadtree_golf          # S on image.bin
python3 harness/evaluate.py candidates/006_simprior --ref         # reference-type rows
python3 lab/eval_on.py heldout/h2/image.bin candidates/005_quadtree_golf
cc -O2 -DT=4 -DRES=1 -DSHC=1 -DSHV=1 -o anchor lab/anchor.c -lm && ./anchor image.bin
cc -O2 -o prime lab/prime.c -lm && ./prime image.bin lab/train/t01/image.bin   # priming
```
