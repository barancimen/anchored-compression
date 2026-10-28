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
