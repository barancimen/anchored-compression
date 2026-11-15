#!/usr/bin/env python3
"""Where a coder spends its bits on a labelled image (lab tool).

usage: where_bits.py PREFIX COSTFILE [COSTFILE ...]
  PREFIX.npz from lab/gen_labeled.py; COSTFILE = float32 bits per pixel (cost.f32 written by
  lab/costmap.c), one per image (several prefixes: PREFIX1,PREFIX2 with matching cost files).
For each label: share of pixels, bits per pixel, total bytes, and the excess over the ideal
cost of a pixel of plain background (Szego limit from lab/gauss_orders.py at the image's own
scale; for the inside of opaque objects, the ideal of that smoother field). Objects whose
shape is known add almost nothing to the ideal cost, so the excess is roughly what a model
that knew the objects could save, before paying for describing them.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gauss_orders as g  # noqa: E402

NAMES = ["background", "opaque inside", "shape outline", "line", "ring", "texture disc"]
prefixes = sys.argv[1].split(",")
costs = sys.argv[2:]
tot = {}
for pre, cf in zip(prefixes, costs):
    z = np.load(pre + ".npz")
    cat, contrast, scale = z["cat"].ravel(), z["contrast"].ravel(), float(z["scale"])
    q = np.fromfile(cf, np.float32).astype(np.float64)
    bg, inside = g.spectra(scale)
    ideal = {0: g.szego(bg), 1: g.szego(inside)}
    print("%s: scale %.3f, %.0f bytes ideal-bg %.4f bpp, ideal-inside %.4f bpp, coder %.4f bpp" %
          (os.path.basename(pre), scale, q.sum() / 8, ideal[0], ideal[1], q.mean()))
    print("  %-14s %6s %7s %8s %9s" % ("label", "share", "bpp", "bytes", "excess B"))
    for c in range(6):
        sel = cat == c
        n = sel.sum()
        if not n:
            continue
        ex = (q[sel] - ideal.get(c, ideal[0])).sum() / 8
        print("  %-14s %6.3f %7.3f %8.0f %9.0f" % (NAMES[c], n / cat.size, q[sel].mean(), q[sel].sum() / 8, ex))
        t = tot.setdefault(c, [0, 0.0, 0.0])
        t[0] += n
        t[1] += q[sel].sum() / 8
        t[2] += ex
    for c in (2, 3, 4):
        sel = cat == c
        lo, hi = sel & (contrast < 1), sel & (contrast >= 1)
        if sel.any():
            print("    %-12s contrast <1 level: %5d px %.3f bpp | >=1 level: %5d px %.3f bpp" %
                  (NAMES[c], lo.sum(), q[lo].mean() if lo.any() else 0, hi.sum(), q[hi].mean() if hi.any() else 0))
if len(prefixes) > 1:
    n_all = sum(t[0] for t in tot.values())
    print("all images: label, share, bytes per image, excess bytes per image")
    for c, (n, b, ex) in sorted(tot.items()):
        print("  %-14s %6.3f %8.0f %9.0f" % (NAMES[c], n / n_all, b / len(prefixes), ex / len(prefixes)))
