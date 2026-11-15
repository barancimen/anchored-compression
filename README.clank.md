# Anchor-point compression of generated 4-bit images (autoresearch log)

**Question.** Can a lossless compressor that first stores a set of *anchor* pixels and then
codes the rest of the image *given* those anchors beat a plain raster-scan coder on images
from `gen_image.py`, measured by Hutter Prize rules?

**Short answer: no, not on these images.** Every anchor design tried ties a raster coder at
best in archive size, and the anchor machinery always adds program code, which counts twice.
The best anchor compressor, `011_lattice_template` (S 41757), uses only 16 anchors and is 310
bytes behind the raster coder with the same model, `007_template_ref` (S 41447). Before that
model existed the gap was 313 bytes (S 41929 vs 41616). The best quadtree is 697 bytes behind
the old reference (S 42313), and the stronger model does not close the gap. These gaps hold on
four held-out images. Sending the objects' parameters as anchors saves 309 archive bytes but
costs 3523 bytes of S. Training the model on images simulated inside the decompressor saves
119 archive bytes but costs 1809 bytes of S.

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
Kinds: **anchor** = this project's anchor candidates; **model** = this project's changes to the
raster reference's model, without anchors (scored with `--ref`); **ref** = the earlier raster
project's algorithms, copied here (`000` is its best, `r00` and `r03..r16` the rest); **sim** =
source simulation.

| S | candidate | kind | archive | zip | what it is |
|---|---|---|---|---|---|
| **41447** | `007_template_ref` | model | 39051 | 1173 | `000` + non-local means as a second Gaussian input (idea 4 of `ideas.md`) |
| **41616** | `000_raster_ref` | ref | 39512 | 1027 | raster context mixing, 112-neighbour RLS (earlier project's v15) |
| 41616 | `009_fastls_ref` | model | 39040 | 1263 | `007` + a small RLS that forgets fast (idea 2) |
| 41638 | `r16_cm_v14` | ref | 39534 | 1027 | v14: context re-searched under the transposed scan |
| 41652 | `r15_cm_v13` | ref | 39548 | 1027 | v13: transposed + mirrored scan |
| 41703 | `010_selective_ref` | model | 39049 | 1302 | `009` + object pixels kept out of the slow RLS's training (idea 3) |
| **41757** | `011_lattice_template` | anchor | 39065 | 1321 | best anchor compressor: `003` with `007`'s model |
| 41821 | `r14_cm_v12` | ref | 39737 | 1017 | v12: micro-golf |
| 41847 | `r13_cm_v11` | ref | 39737 | 1030 | v11: 5th counter context |
| 41919 | `r12_cm_v10` | ref | 39835 | 1017 | v10: golf |
| 41923 | `r11_cm_v9` | ref | 39819 | 1027 | v9: golf |
| 41929 | `003_lattice_lookahead_t128` | anchor | 39529 | 1175 | 16 verbatim anchors + look-ahead |
| 41950 | `r10_cm_v8` | ref | 39810 | 1045 | v8: searched counter contexts |
| 41968 | `001_lattice_lookahead` | anchor | 39566 | 1176 | 64 verbatim anchors + look-ahead |
| 42044 | `008_kriging_ref` | model | 39026 | 1484 | `007` + a predictor computed from the generator's spectrum (idea 1) |
| 42197 | `r09_cm_v7` | ref | 40031 | 1058 | v7: golf step (first scored here) |
| 42250 | `r08_cm_v6` | ref | 40038 | 1081 | v6: RLS radius 6 (84 inputs) |
| **42313** | `005_quadtree_golf` | anchor | 39895 | 1184 | best quadtree: anchor lattice + refinement |
| 42417 | `004_quadtree_ctx` | anchor | 39939 | 1214 | quadtree, transposed, searched contexts |
| 42460 | `012_quadtree_template` | anchor | 39634 | 1388 | `005` + non-local means |
| 42512 | `r07_cm_v5` | ref | 40204 | 1125 | v5: golfed v4 |
| 42644 | `002_quadtree` | anchor | 40200 | 1197 | first quadtree |
| 42815 | `013_quadtree_edges` | anchor | 39819 | 1473 | `005` + edge-directed refinement (idea 6) |
| 42854 | `r06_cm_v4` | ref | 40202 | 1297 | v4: RLS radius 4 (41 inputs) |
| 43112 | `015_quadtree_kriging` | anchor | 39880 | 1591 | `005` + weights computed from the spectrum for each pass |
| 43425 | `006_simprior` | sim | 39393 | 1987 | `000` trained on 3 simulated images before coding |
| 43803 | `r05_cm_v3` | ref | 40451 | 1647 | v3: 4 contexts + Gaussian/OLS input |
| 44334 | `r03_cm_v1` | ref | 42016 | 1130 | v1: 8 neighbour contexts, no linear predictor |
| 44970 | `014_object_anchors` | anchor | 38742 | 3089 | `007` + object parameters as anchors (idea 5) |
| 45786 | `r04_cm_ols` | ref | 41628 | 2050 | v2: v1 + OLS predictor (first scored here) |
| 125336 | `r00_store` | ref | 125000 | 151 | stored verbatim |

`r04_cm_ols` is copied unchanged except for one removed line that dumped a cost map to a
scratch path (instrumentation). The earlier project's sandbox check (a "compressor" that
reads the original from an absolute path) is not an algorithm and was not copied.

## Held-out images (`lab/generalization.tsv`; never used for tuning)

| image | `000` reference | `003` 16 anchors | `005` best quadtree | `006` simulation | `007` non-local means | `011` anchors + `007`'s model |
|---|---|---|---|---|---|---|
| h1 | 41887 | 42200 (+313) | 42611 (+724) | 43748 (+1861) | 41771 (-116) | 42081 (+310) |
| h2 | 41174 | 41494 (+320) | 41902 (+728) | 43017 (+1843) | 41032 (-142) | 41346 (+314) |
| h3 | 41202 | 41515 (+313) | 41938 (+736) | 43109 (+1907) | 41010 (-192) | 41321 (+311) |
| h4 | 40261 | 40573 (+312) | 40935 (+674) | 42113 (+1852) | 40070 (-191) | 40382 (+312) |

Differences are against `000`, except for `011`, which is compared with `007`, the raster coder
with the same model. The gaps on `image.bin` (+313, +697, +1809, -169 and +310) hold on every
held-out image.

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

## Headroom and next ideas (`ideas.md`)
A later analysis measured how far the reference coder is from the ideal on synthetic images
whose objects are known. With ideal prediction, every coding order costs the same on the
background to within 0.003 bpp, so no choice of anchor pixels can gain there. The coder is
about 3 KB per image above the ideal: about 960 bytes in the background, where one shared
predictor also has to learn the objects, and about 2 KB at outlines, lines, rings and striped
discs. Knowing the additive objects exactly would save about 800 bytes. All six ideas in
`ideas.md` were then built and scored (its Results section). Only template matching paid:
non-local means as a second Gaussian input (`007_template_ref`) saves 461 archive bytes and
169 bytes of S. With that model in both, lattice anchors still trail raster by 310 bytes.

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
| `candidates/r00`, `r03..r16` | the earlier raster project's algorithms |
| `lab/anchor.c`, `lab/sweep.sh` | research coder with every anchor order as a `-D` knob |
| `lab/qt.c`, `lab/qtsweep.sh`, `lab/qtsearch.py` | quadtree template, knob sweeps, counter-context search (`lab/qtsearch.log`) |
| `lab/prime.c`, `lab/sim.c`, `lab/train/` | priming probe, C port of the generator, 8 generator images for training |
| `lab/notes.md` | all probes with numbers, including rejected ones |
| `lab/eval_on.py`, `lab/generalization.tsv` | scoring on held-out images |
| `ideas.md` | headroom analysis, literature and ranked next ideas |
| `lab/gauss_orders.py` | ideal cost of any coding order on the generator's background |
| `lab/gen_labeled.py`, `lab/costmap.c`, `lab/where_bits.py` | labelled synthetic images, per-pixel cost of the reference model with probes, cost by label |
| `lab/model.c` | the reference model with ideas 1-4 of `ideas.md` as `-D` knobs |
| `lab/objects/` | object parameters as anchors: oracle, renderer, fitter |

## Reproduce

```sh
export HP_SCRATCH=/some/dir/outside/this/repo     # decompression runs with the repo hidden
python3 harness/evaluate.py candidates/005_quadtree_golf          # S on image.bin
python3 harness/evaluate.py candidates/006_simprior --ref         # reference-type rows
python3 lab/eval_on.py heldout/h2/image.bin candidates/005_quadtree_golf
cc -O2 -DT=4 -DRES=1 -DSHC=1 -DSHV=1 -o anchor lab/anchor.c -lm && ./anchor image.bin
cc -O2 -o prime lab/prime.c -lm && ./prime image.bin lab/train/t01/image.bin   # priming
cc -O2 -DTM=4 -DTMR=4 -DTMN=12 -DTMH=2 -o model lab/model.c -lm && ./model <image.bin | wc -c  # 007's model
```
