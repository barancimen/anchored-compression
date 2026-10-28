# Autoresearch: anchor-point compressor for image.bin

Goal: minimise the Hutter-Prize total size **S** for `image.bin` (500x500, 4-bit pixels,
125000 bytes packed, made once by `gen_image.py` from unrecorded `os.urandom` entropy) with a
compressor built around **anchor points**:
1. **Select** a set of anchor pixels. Their positions are either implied by a rule the
   decoder repeats from data it already has (free), or chosen by the encoder and transmitted
   (paid for in the archive: a sparse map, a tree, ...). The number of anchors may be fixed or
   chosen by the compressor.
2. **Store** the anchor values: verbatim ("in full", 4 bits each) or with a model. The loop
   decides which pays.
3. **Code the rest** of the image with a model conditioned on the anchors, i.e. using the
   patterns that are most informative given them.

The method must generalise to other images from the same procedure: `heldout/h1..h4` are
four more images from unchanged copies of `gen_image.py`. They are never used for tuning;
milestones are checked on them with `lab/eval_on.py`.

## Fixed (never edited during the loop)
- `image.bin`, `gen_image.py`, `heldout/`.
- `harness/evaluate.py` + `harness/minify.py`: the scorer. It implements the Hutter rules:
  `S = L(C)+L(DA)` (self-extracting) | `L(P)+L(P)+L(A)` (one program, C = D) |
  `L(C)+2L(D)+L(A)` (separate D). Programs are zips of source (or executables).
  Every character of every executed command line is charged. Decompression runs in a fresh
  sandbox without network, with the project directory hidden, and must reproduce
  `image.bin` byte for byte.

## Reference points
`candidates/000_raster_ref` is the best compressor of an earlier, separate project (raster
scan context mixing, no anchors), copied here unchanged and scored with `--ref`. It is the bar
to compare against and never becomes `best.json`. Its model and coder may be reused inside
anchor candidates. `candidates/r00..r16` are that project's other algorithms (store, LZMA,
and every context-mixing version), copied and scored the same way. `006_simprior` (training
the model on simulated images first) is not an anchor method either, so it is scored with
`--ref` too.

## Design families
- **A. Anchors as side information:** a lattice of anchors is stored first; the rest is
  coded in scan order by a model whose stencil also reads ahead into an interpolation of the
  anchors (`001_lattice_lookahead`, built on the reference's model).
- **B. Anchors + hierarchical refinement:** anchors are a coarse lattice and every other
  pixel is coded coarse-to-fine with known neighbours on all sides (`002_quadtree`).
A family's next experiment may branch from that family's best candidate even when the
other family holds `best.json`; `results.tsv` still marks it against the global best.

## Loop
1. Pick one idea (anchor selection, anchor storage, conditional model, coder, golf, ...).
2. Copy the current best candidate to `candidates/NNN_short_name/`, apply the idea.
3. Optionally probe quickly in `lab/` (archive size of a local build) to sweep knobs.
4. Score: `python3 harness/evaluate.py candidates/NNN_short_name -m "what changed"`.
   It appends a row to `results.tsv` and marks it `keep` only if it is verified and beats
   `best.json`; otherwise it is `discard` / `mismatch` / `crash`.
5. Commit the candidate and the new `results.tsv` row (kept or not). Record probes,
   including rejected ones, in `lab/notes.md`.
6. Next experiment branches from the best candidate (`best.json`). Repeat.

Trade-off to keep in mind: with C = D, each zipped source byte costs 2 bytes of S, so a
model feature must save more than twice its compressed source size.
