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
