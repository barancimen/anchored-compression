# Anchor-point compression of generated 4-bit images (autoresearch log)

**Question.** Can a lossless compressor that first stores a set of *anchor* pixels and then
codes the rest of the image *given* those anchors beat a plain raster-scan coder on images
from `gen_image.py`, measured by Hutter Prize rules?

**Short answer.** Not on these images. Every anchor design tried ties a raster coder at best
in archive size, and anchors always cost extra program code. Two measurements explain it:
1. **Coding order doesn't change the total information.** By the chain rule, coding the
   anchors first and the rest given them has the same ideal total as any other order. In
   practice, the quadtree pays 2341 bytes more than raster on its anchors and coarse levels
   and saves 2081 bytes on the two-sided finest level (net +260).
2. **Anchors can't pay for transmitted positions.** Handing the decoder anchors *for free*
   (no position or value cost) makes the rest cheaper by little more than those pixels' own
   cost. Choosing them by edge strength or by coding cost helps the other pixels by at most
   ~1100 bits per 4000 anchors, while transmitting them would need ~45600 bits.

## Results on `image.bin` (Hutter total S, bytes; lower is better)

S = 2 x zipped program + archive + charged command-line characters (C = D form).

| candidate | design | archive | program (zip) | S | |
|---|---|---|---|---|---|
| `000_raster_ref` | reference: earlier project's raster coder, no anchors | 39512 | 1027 | **41616** | ref |
| `001_lattice_lookahead` | A: 8x8 verbatim anchor lattice, look-ahead into its interpolation | 39566 | 1176 | 41968 | keep |
| `002_quadtree` | B: even/even anchor lattice + quadtree refinement | 40200 | 1197 | 42644 | discard |

(Full log: `results.tsv`. Probes that never became candidates: `lab/notes.md`.)

## Layout

| path | contents |
|---|---|
| `image.bin`, `image.png` | the benchmark (copied from the earlier project; same image, so S is comparable) |
| `gen_image.py` | the image generator, unchanged |
| `heldout/h1..h4` | four more images from unchanged copies of `gen_image.py`, never used for tuning |
| `harness/evaluate.py` | Hutter-rule scorer (sandboxed decompression, charged command lines, `--ref`) |
| `program.md` | the loop protocol and the two design families |
| `candidates/NNN_*` | `src.c` (readable), `p.c` (minified, measured), `cand.json`, `submission/P.zip` |
| `lab/anchor.c` | research coder with every anchor order as a `-D` knob (`lab/sweep.sh` runs variants) |
| `lab/qt.c`, `lab/qtsearch.py` | quadtree template and its counter-context search |
| `lab/notes.md` | all probes with numbers, including rejected ones |
| `lab/generalization.tsv` | scores on the held-out images |

## Reproduce

```sh
export HP_SCRATCH=/some/dir/outside/this/repo     # decompression runs with the repo hidden
python3 harness/evaluate.py candidates/001_lattice_lookahead       # S on image.bin
python3 lab/eval_on.py heldout/h2/image.bin candidates/001_lattice_lookahead
cc -O2 -DT=4 -DRES=1 -DSHC=1 -DSHV=1 -o anchor lab/anchor.c -lm && ./anchor image.bin
```
