# Lab notebook (quick probes; numbers are archive bytes of image.bin, not official S)

Probe tool: `lab/anchor.c` (research coder, knobs as -D defines), `lab/sweep.sh`.
Reference: the earlier project's raster coder (`candidates/000_raster_ref`) has a 39512-byte
archive (S 41616) with a 112-neighbour RLS and heavily searched contexts.

## 1. Regular lattice anchors + quadtree-style refinement (lab/anchor.c)
Anchors = lattice of step T coded first; then per level t = T..2: cell centres (pass D),
then edge midpoints (pass A). Defaults: RLS stencils 12 (anchors) / 16 (centres) /
26 (axis), 5 counter contexts, per-group RLS, variance buckets and counter hashes.

| probe | result | verdict |
|---|---|---|
| T=16 first version (verify ok) | 40868 | baseline |
| T sweep, same model: T=1 (plain raster) / 2 / 4 / 8 / 32 / 64 / 256 | 40511 / 40529 / 40701 / 40818 / 40885 / 40886 / 40884 | raster wins by ~350 at large T |
| raster T=1 with 7.5-radius stencil (44 causal neighbours) | 40134 | stencil size matters more than order |
| T=16 anchors verbatim (4 bits) instead of modelled | +148 | modelled anchors |
| T=16: RLS on residual from mean of 4 nearest (RES) | -78 | keep |
| per-group mixer weights (MIXG) | +245 | reject: too few samples per group |
| variance bucket + 2 x spread of 4 nearest (KS=2, KQ=1.5) | -46 | keep |
| stencils RD/RX 3.2 -> 4.5 | -32 | keep |
| 6th counter context (a0,a2,a4) / only 4 contexts | -46 / +168 | keep 6 |
| RLS forgetting .9999 / regulariser 30 | +40 / +18 | keep .9997 / 100 |
| combined (RES, KS, NC=6, R 4.5) | 40626 | |
| share across levels: RLS / variances / counters / all three | +38 / -103 / -192 / -201 | share variances + counters |
| T sweep with sharing: T=1 / 2 / 4 / 8 / 16 / 32 / 64 | 40398 / 40380 / 40372 / 40378 / 40383 / 40390 / 40396 | all within 30 bytes; T=4 best |

Reading: by the chain rule every coding order has the same ideal total; an order only
wins if its online model learns better. Separate statistics per level starve the coarse
levels (few pixels); sharing them across levels removes most of the hierarchy's handicap.

## 2. Anchors + raster coding of the rest with look-ahead (RMODE=1)
Anchors = sparse lattice of step T; unknown pixels are filled with the bilinear
interpolation of the anchors; the rest is coded in raster order and its RLS stencil also
reaches ahead (radius RL) into the fill. Noise level of these comparisons: about +-30 bytes
(e.g. an 80-cell cap that silently dropped the look-ahead still moved results by 20-30).

| probe | result | verdict |
|---|---|---|
| raster T=1, causal radius 5.5 | 40074 | baseline for this block |
| look-ahead, contexts from the 4 nearest incl. fill estimates: T=4/8/16/32 | 40406 / 40392 / 40371 / 40346 | worse: fill estimates poison contexts |
| contexts from causal pixels only: T=8 / 16 / 32 / 64 | 40273 / 40119 / 40035 / 40032 | T>=32 ties raster |
| causal radius 7.5 (88 cells) + look-ahead 4.5: raster / T=16 / 32 / 64 | 40000 / 40066 / 39978 / 39969 | within noise of raster |
| T=32 anchors verbatim instead of modelled | +22 | |

## 3. Where the bits are (raster cost map, causal radius 5.5)
Top 0.1% / 1% / 5% / 25% of pixels carry 0.8% / 5.4% / 18% / 53% of the bits: the cost is
diffuse. The cost map shows the cloudy background as uniform noise, object edges as thin
lines, and flat object interiors as cheaper patches. Noise alone in gen_image.py (sigma
.28 before flooring) has an entropy of 0.729 bpp = 22767 bytes, so ~17 KB is field detail.

## 4. Upper bound for selected (irregular) anchors
Anchors given to the decoder for FREE (mask, not transmitted), then raster + look-ahead
(T=64 base). saved = drop in the rest's bits; own = those pixels' own bits in the base run;
look-ahead = saved - own (the benefit to the other pixels); needed = what transmitting
them would cost at least (log2 C(250000, N) for positions + 4 bits per value).

| select | N | saved | own | look-ahead | needed |
|---|---|---|---|---|---|
| highest cost | 250 / 1000 / 4000 | 2513 / 8166 / 26058 | 2695 / 8328 / 24943 | -182 / -162 / +1115 | 3851 / 13406 / 45623 |
| strongest gradient (Sobel) | 250 / 1000 / 4000 | 727 / 2724 / 10313 | 660 / 2469 / 9928 | +67 / +255 / +385 | same |
| random | 250 / 1000 / 4000 | 313 / 1265 / 4878 | 329 / 1336 / 5061 | -16 / -71 / -183 | same |

Knowing a pixel early helps the other pixels almost nothing beyond its own cost (at most
~1100 bits for 4000 anchors), so transmitted positions can never pay: selection by edges or
by cost misses break-even by 10x or more. Only positions the decoder can derive for free
(lattices, or rules on decoded data) are viable.

## 5. Row hierarchy (ROWS=1) and where the quadtree loses
| probe | result | verdict |
|---|---|---|
| whole anchor rows every T rows, rows refined by halving: T=8 / 4 / 2 | 40528 / 40403 / 40276 | worse than raster (40000) |
| T=2 rows with radius 7.5 stencils | 40236 | reject |

Quadtree (T=4, best lab knobs, 40261) vs raster (radius 7.5, 40000), same pixels:

| pixels | n | quadtree bpp | raster bpp | quadtree - raster |
|---|---|---|---|---|
| anchors (x, y % 4 == 0) | 15625 | 1.952 | 1.286 | +1299 B |
| t=4 centres | 15625 | 1.576 | 1.297 | +544 B |
| t=4 edge midpoints | 31250 | 1.413 | 1.286 | +498 B |
| t=2 centres | 62500 | 1.267 | 1.275 | -63 B |
| t=2 edge midpoints | 125000 | 1.149 | 1.278 | -2018 B |
| all | 250000 | 1.288 | 1.280 | +260 B |

The two-sided finest level is much cheaper than raster, but the anchors and coarse levels
pay almost exactly that back: information moves between levels, it does not disappear.

## 6. Candidates (family A: look-ahead on the reference model; family B: quadtree)
| probe | archive | verdict |
|---|---|---|
| 001 build: T=64, look-ahead L=135 | 39566 | scored S 41968 |
| 001 with L=112 (anchors, no look-ahead) / L=120 / L=150 | 39565 / 39564 / 39569 | look-ahead adds nothing |
| 001 with T=32 / T=128 | 39721 / 39529 | fewer anchors is better |
| 002 quadtree, T=2, radius-6 stencils | 40200 | scored S 42644 |
| 002 under the 8 scan orientations | 40064 (transposed) .. 40217 | transposed: -136 |
| quadtree context search (lab/qtsearch.py, 3 rounds, log in lab/qtsearch.log) | 40064 -> 39939 | keep: (v9,v1,v0) (v0,v1,v10) (v2,v3) (m,q) (v7,v11,v4,v5) |
| same contexts on held-out h1 / h2 / h3 | -128 / -159 / -154 | transfers |
| 004 = 002 + transposed + searched contexts | 39939 | scored S 42417 (best of family B) |
| 003 = 001 with T=128 (16 anchors) | 39529 | scored S 41929 (new best) |
| quadtree knobs (lab/qtsweep.sh): radius 6->7 / lr / counter limit / regulariser | within -13..+10 | noise |
| RLS forgetting .9997 -> .9995 / .9993 / .9999 | -27 / -34 / +101 | .9995 |
| forgetting .9995 + bucket scale 3->2 + variance rate 1/500 -> 1/300 | 39939 -> 39895 | keep (no code cost) |
| same on held-out h1 / h2 / h3 / h4 | -44 / -43 / -31 / -30 | transfers |
| 005 = 004 golfed + those constants: p.c 1846 -> 1737 chars, zip 1214 -> 1184 | 39895 | scored S 42313 (best of family B) |

## 7. Simulating the source (learning moved out of the archive)
Idea: the decoder generates images from the same recipe and trains the model on them before
coding, so less learning is paid for inside the archive. Probe: `lab/prime.c` (the raster
reference model, primed on images, then coding the target).

| probe | image.bin archive | verdict |
|---|---|---|
| cold (= 000_raster_ref) | 39512 | |
| primed on 1 / 2 / 4 real generator images (perfect simulator, lab/train) | 39423 / 39422 / 39422 | ceiling ~ -90 |
| counter limit 255: cold / 4 / 8 images | 39531 / 39444 / 39456 | slower learning doesn't use more data |
| counter limit 255 + RLS forgetting .9999: cold / 4 / 8 | 39569 / 39482 / 39491 | |
| RLS forgetting .9999: cold / 4 / 8 | 39550 / 39460 / 39468 | |
| C port of gen_image.py (lab/sim.c), primed on sim image 0 / 1 / 0-2 | 39427 / 39512 / 39421 | sim images work as well as real ones |
| held-out h1: cold / sim-primed / real-primed | 39783 / 39710 / 39667 | |
| held-out h2: cold / sim-primed / real-primed | 39070 / 38983 / 38978 | |

The model adapts fast (RLS memory ~3300 pixels, counters ~80 updates), so earlier images
only give it a warm start: at most ~90-120 bytes. The simulator (3 random fields by FFT,
6 object kinds, percentile scaling) is 960 zipped bytes of extra program, counted twice.
Scored as `006_simprior`: archive 39393 (-119), program 1987 zipped (+960), S 43425 (+1809).

## 8. Headroom: ideal costs, where the bits go, oracle probes (2026-11-15)
Tools: `lab/gauss_orders.py` (ideal cost of any coding order on the generator's Gaussian
background), `lab/gen_labeled.py` (gen_image.py with a seed and per-pixel labels),
`lab/costmap.c` (the reference model with a per-pixel cost dump and -D probes),
`lab/where_bits.py`. Synthetic seeds 3 and 2 (scales 2.15 and 2.58). Write-up: `ideas.md`.

| probe | result | verdict |
|---|---|---|
| ideal background cost, scale 2.5: raster r=2 / 3 / 7 / 15, Szego limit | 1.2279 / 1.2107 / 1.2035 / 1.2024 / 1.2021 bpp | radius 7 is at the limit |
| quadtree from step 2 / 4 / 8, r=6; step 2 with r=3 / 2 | 1.2027 / 1.2029 / 1.2035; 1.2052 / 1.2104 | orders tie; hierarchy helps only with tiny stencils |
| lattice 16 / 4 + raster look-ahead, r=7 | 1.2048 / 1.2027 | ties |
| inside opaque shapes: own weights / background weights | 0.784 / 0.839 bpp | one shared predictor costs 0.055 bpp inside |
| batch least squares, radius 7, clean background (seed 3), first run: 15000-pixel sample in-sample / other 15000 | 1.127 / 1.158 vs 1.122 bpp | the sample was lucky; next rows |
| same fitted on all 95664 clean pixels / weights computed from the spectrum, no fitting | 1.141 / 1.139 vs 1.122 bpp | the ideal is approximate for one image; see the 8-image row |
| seed 2, fitted on all clean pixels / background-only image (seed 3, scale 2.69), fitted / spectral weights / coder | 1.238 vs 1.221; 1.264 / 1.272 / 1.268 vs 1.247 bpp | the gap is ~380 B of background on each of these two seeds; learning online costs ~0.005 bpp |
| fitted radius-7 predictor minus theory, clean background, 8 synthetic images s1 / s2 / s3 / s4 / s5 / s6 / s7 / s8, and b3 (seed 3 without objects) | +0.002 / +0.016 / +0.020 / +0.019 / +0.031 / -0.012 / +0.022 / +0.020 / +0.017 bpp; spectral weights 0.0003 to 0.0016 worse | scatters by image, mean of s1-s8 +0.015: the ideal is good to about 0.03 bpp for one image; negative values rule out rounding alone |
| cost by label (mean of 2 images): excess over ideal, bytes | background 961, inside 257, outlines 731, lines 424, rings 210, stripes 400 | 2983 B = 7.9% over the Gaussian ideal; on these two images ~380 B of the background's share is out of reach of a fixed linear predictor |
| background excess by distance to objects | ~+0.05 bpp even 16+ px away; border 8 px +0.07..0.1 | spread out, not edge pollution |
| Gaussian input alone vs mixed, clean background (seed 3) | 1.244 vs 1.170 bpp | counters + mixer repair the predictor |
| RLS squared error on clean background vs theory: real image / background-only image | +12.6% / +4.6% (+3.4% late in the scan); the best fitted predictor is itself +2.8% / +2.3% | objects drag the shared weights |
| -DVQ=.167 (variance minus rounding) | archive 36156 -> 36155 | Gaussian alone better, mixer already compensates |
| -DGSN=4 (4 Gauss-Seidel sweeps per pixel) | 36156 -> 36140 | not the solver |
| -DROB=4 / -DROB=1 (robust RLS updates) | 36156 -> 36184 / 36219 | reject |
| -DTMPL oracle: additive objects known (seed 3 / 2) | -773 / -818 B | stripes -69%, rings -69%, lines -28%, outlines -23% of excess |
| -DCLS oracle: one RLS weight set per true region (seed 3 / 2) | -319 / -229 B; background excess 961 -> 843 B | perfect segmentation is worth ~0.7% |

## 9. Ideas 1-6 of ideas.md (2026-11-15)
Probe tool: `lab/model.c` (the reference model with -D knobs; without knobs it writes the
reference's archive, 39512). Archive bytes of image.bin, local -O2 builds.

| probe | archive | verdict |
|---|---|---|
| reference (lab/model.c, no knobs) | 39512 | same as 000 |
| idea 1: kriging input (spectral weights, 112 neighbours, scale 2.5) | 39499 | -13 |
| same with scale 2.2 / 2.8 | 39499 / 39499 | the scale hardly matters |
| + mixer weight set chosen by the locally best Gaussian input (MIXC) | 39454 | |
| + its sd fixed at .5 instead of buckets (KFIX) | 39428 | best idea-1 variant, -84 |
| + its counter context (prediction, bucket) / variance minus rounding | 39458 / 39461 | no help |
| RLS regularised towards the kriging weights, strength 3000 / 30000 / 300000 | 39554 / 39748 / 39906 | reject: adaptivity matters more |
| mixer weight sets per bit position (with and without kriging) | 39565 .. 39570 | reject |
| mixer learning rate .0005 / .002 (kriging + MIXC) | 39494 / 39454 | .001 stays |
| diagnostic, seed 3, clean background: kriging Gaussian alone with ideal sd (offline) / in the coder as the only input / full mix | 1.144 / 1.197 / 1.168 bpp | one global mixer weight cannot calibrate it for both background and edges |
| idea 2: one small RLS, 12 nearest @ .99 / 24 @ .998 / two predictors | 39495 / 39488 / 39523 | -24 at best |
| idea 3: slow RLS learns with weight .3 / .1 where the small RLS did better | 39496 / 39518 | reject |
| idea 4: template matching, best match of the 6 nearest within 6 columns, 6 rows | 39380 | -132 |
| template of 10 / 12 nearest, 4x4 window | 39351 / 39351 | |
| mean of the 2 / 4 / 8 best matches (12 nearest, 4x4) | 39262 / 39203 / 39327 | 4 best |
| match on shape (mean removed) / counter context of the match | 39361 / 39349 | no help |
| non-local means, weights exp(-d/h), 12 nearest, 4x4 window, h = 1.4 / 2 / 4 / 8 | 39099 / 39051 / 39078 / 39242 | h = 2: -461, kept as 007 |
| non-local means, window 6x6 / 8x8 | 39046 / 39063 | 4x4 is enough |
| on top of non-local means: kriging (KFIX) / + MIXC / small RLS / + selective training | 39025 / 39029 / 39040 / 39046 | -26 / -22 / -11 / -5 |
| on top of non-local means: shared variance buckets / no RLS input | 39121 / 39308 | own buckets; the RLS still pays |
| cost by label with the 007 model, seeds 3 and 2 | 35488 / 39162 (were 36156 / 39517) | excess per image: background 895 (961), inside 176 (257), outlines 564 (731), lines 349 (424), rings 103 (210), stripes 387 (400) |
| quadtree 005 + non-local means over same-parity pixels within 4 / 8 | 39699 / 39634 (005: 39895) | candidate 012 uses 8 |
| quadtree 005 + edge-directed interpolation (NEDI) input | 39819 | -76: candidate 013 |
| quadtree 005 + kriging weights per pass | 39880 | -15: the T = 2 quadtree has no learning handicap to remove; candidate 015 |
| lattice 003 + non-local means | 39065 (003: 39529) | candidate 011; 14 bytes behind 007 with the same model |

## 10. Object parameters as anchors (idea 5; lab/objects/)
Tools: `lab/objects/oracle.c` (007 with the true layer of lab/gen_layers.py given for free),
`render.c` (object list to layer, gen_image.py's soft masks, opaque objects cover what is
below), `fit.c` (greedy encoder-side fitter), `conv.py`. Archive bytes, 007 model.

| probe | seed 2 | seed 3 | verdict |
|---|---|---|---|
| 007 model, no layer | 39163 | 35491 | |
| true additive layer subtracted for the RLS / for non-local means / for both | 38531 / 38758 / 38434 | 34909 / 34736 / 34553 | both inputs should see it |
| all true objects (opaque ones as additive with fitted amplitude), both inputs | 38099 | 33579 | -1064 / -1912 bound |
| same without covering by opaque objects | +51 | +355 | keep covering |
| parameters quantised to 1/8 px / 1/4 px | +60 / +135 | +135 / +350 | 1/8 px |
| bound minus parameters (about 70 bits per object) | about -790 | about -1540 | |
| full lab fitter (15.5 min at -O2) | | 34405 before parameters | 57% of the oracle, about -650 net |

Candidate 014 (cut-down search so the plain-cc build finishes in about 11 min): image.bin
38742 (-309 vs 007), seed 3 34948 (-543), held-out h1 38875 (-500), 40-50 objects found (16 of
the 43 true objects on seed 3 matched in kind within 3 px; the rest are fragments that still
help). Its p.c deflates to 3035 bytes against 1082 for 007 with zlib at level 9 (2985 and 1069 with
zopfli, as scored): the decoder side alone (renderer,
object coding, layer in the model) is about 540 deflated bytes, i.e. about 1100 bytes of S,
as much as a perfect fitter saves on average before the fitter's own code.

## 11. The generator replayed from its seed (016_seed_replay; not ranked) (2026-11-15)
Question: what would S be if the 32 bytes that gen_image.py draws from os.urandom were known?
The archive is then those bytes (or the generator state numpy derives from them) and the
program is gen_image.py without numpy. No compressor can find the bytes, so this is not a
candidate: image.bin cannot be coded, and the result is never ranked.
Tools: `lab/gen_seeded.py` (the unchanged gen_image.py run with a known seed in place of
os.urandom), `lab/seed_check.py` (the unchanged harness on such an image, with the compress
command replaced by writing the archive from the seed; log in `lab/seed_replay.tsv`). Test
seeds are sha256("seed replay test n"), n = 1..6. S = 2 x zip + 32 + 50.

| probe | zip | S | verdict |
|---|---|---|---|
| first plain port (numpy seed hashing, ziggurat, direct two-pass sum), 6 seeds identical | 1661 | 3404 | works at the first try |
| same, shortened | 1472 | 3026 | |
| archive = PCG64 state and increment after seeding, no seed hashing in the program | 1309 | 2700 | keep |
| plain `double` in place of the D macro, and pass 1 summing straight into B | 1299 | 2680 | keep |
| plain `for` and the number 250000 in place of the F and N macros | 1295 | 2672 | keep: deflate codes the repeats as cheaply |
| tail of the normal sampler uses two free globals / the three terrain fields are added inside g | -5 / -10, each alone | | keep |
| candidate 016: harness check on 6 seeds, all identical, decode 26-45 s | 1285 | 2652 | scored (status oracle) |
| recursive mixed-radix FFT in place of the direct sum (`lab/seed_fft`), 2 seeds identical, decode 5-6 s | 1437 | 2956 | reject: +304 for time that is not scored |
| archive = the raw os.urandom bytes, hashed as numpy's SeedSequence does (`lab/seed_urandom`), 2 seeds in the harness, 6 locally | 1455 | 2992 | variant: +340 |
| the two percentiles stored in the archive (16 bytes) instead of sorting; size only, not built | 1216 | 2530 | not taken: the archive stays the state alone |
| counter-context string, y fold by a conditional, weights in pass 2, hex multiplier | +1 .. +5 | | reject |

Exactness. The uniform and integer draws and every accept decision of the normal sampler have
to match numpy exactly. The normal values and everything after them only have to agree before
the final floor to 16 levels.
- Ziggurat tables rebuilt in the program by the usual recursion from R = 3.654152885361009 and
  the strip area V. With V = 0.00492867323399 (the 12 digits usually quoted) the fast-accept
  thresholds are off by up to 6e5 in 2^52 and a decision changes with probability 1e-12 per
  draw. With V = 0.004928673233974652 (the value implied by numpy's own table) they differ by
  at most 18 in 2^52 (8e-17 per draw), and 37 of the 256 strip widths differ by up to 76 units
  in the last place. So 14% of the normal draws differ from numpy's by up to 1e-14 of their
  value, and the slow-path test against exp(-x^2/2) then flips with probability about 4e-16
  per draw. In all that is about 5e-16 per draw, or 3e-9 per image of roughly 5 million normal
  draws. Over 2 million draws compared with numpy, the uniforms were identical and no decision
  changed.
- Direct sum against np.fft.irfft2, seeds 1 and 5: the values before rounding differ from
  numpy's by at most 5.0e-14 and 1.2e-13 levels (the same at -O0 and -O2). The nearest pixel
  to a rounding boundary is 1.7e-6 and 1.3e-6 away, so the margin is 7 orders of magnitude and
  the chance of one wrong pixel per image is about 1e-7.
- One bug on the way: `W[i] = cos(t = ...) + 1i * sin(t)` read t before it was assigned
  (unsequenced), which changed 41% of the pixels by one level.
- Not handled, each with probability about 2^-32 per image: a second draw in numpy's bounded
  integer (Lemire rejection), and in `lab/seed_urandom` a seed whose top 32 bits are zero.
