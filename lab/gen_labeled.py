#!/usr/bin/env python3
"""gen_image.py with a seed and per-pixel labels (lab tool; gen_image.py itself stays fixed).

The random draws are the same, in the same order, as in gen_image.py, so the images follow the
same distribution. It also records what each pixel shows, for measuring where a coder spends
its bits (lab/where_bits.py).

usage: gen_labeled.py SEED OUTPREFIX [--background-only]
writes OUTPREFIX.bin (the image.bin format) and OUTPREFIX.npz with
  pix       the 500x500 image (values 0..15)
  cat       label per pixel, the last object drawn there wins:
            0 background, 1 inside an opaque shape (more than 1.5 px from its outline),
            2 within 1.5 px of the outline of a disc, ellipse or rectangle,
            3 line (all of it, up to 1.5 px outside), 4 ring (same), 5 sinusoid texture disc
  contrast  for labels 2..4, the step across the outline in grey levels
  scale     grey levels per field unit (16 / (hi - lo)), and lo, hi
  tmpl      the additive objects' visible part in grey levels (what a decoder could render
            from their parameters); also written as OUTPREFIX.tmpl.f32 for lab/costmap.c -DTMPL
  objects   one row per object: kind index, cx, cy, r, theta, extra1, extra2, opaque, amp
            (kinds: 0 disc, 1 ellipse, 2 rect, 3 ring, 4 line, 5 texture)
"""
import sys

import numpy as np

H = W = 500
seed, out = int(sys.argv[1]), sys.argv[2]
rng = np.random.Generator(np.random.PCG64(seed))
yy, xx = np.mgrid[0:H, 0:W].astype(np.float64)
KINDS = ["disc", "ellipse", "rect", "ring", "line", "texture"]


def grf(beta):
    fy = np.fft.fftfreq(H)[:, None]
    fx = np.fft.rfftfreq(W)[None, :]
    f = np.hypot(fx, fy)
    f[0, 0] = 1.0
    spec = (rng.normal(size=f.shape) + 1j * rng.normal(size=f.shape)) * f ** (-beta / 2)
    spec[0, 0] = 0
    z = np.fft.irfft2(spec, s=(H, W))
    return (z - z.mean()) / z.std()


def soft(d):
    return np.clip(0.5 - d, 0.0, 1.0)


img = 1.00 * grf(3.4) + 0.30 * grf(2.2) + 0.08 * grf(1.0)
img += rng.normal(0, 0.6) * (xx - W / 2) / W + rng.normal(0, 0.6) * (yy - H / 2) / H
cat = np.zeros((H, W), np.uint8)
tmpl = np.zeros((H, W))  # additive objects, field units; opaque objects cover it like the image
step = np.zeros((H, W))  # step across the outline, field units (scaled to grey levels below)
objs = []
nobj = int(rng.integers(30, 45))
for _ in range(0 if "--background-only" in sys.argv else nobj):
    kind = rng.choice(KINDS, p=[0.22, 0.18, 0.22, 0.10, 0.14, 0.14])
    cx, cy = rng.uniform(0, W), rng.uniform(0, H)
    r = rng.uniform(8, 60)
    th = rng.uniform(0, np.pi)
    c, s = np.cos(th), np.sin(th)
    u = (xx - cx) * c + (yy - cy) * s
    v = -(xx - cx) * s + (yy - cy) * c
    e1 = e2 = 0.0
    if kind == "disc":
        d = np.hypot(u, v) - r
    elif kind == "ellipse":
        e1 = rng.uniform(0.3, 0.9)
        a, b = r, r * e1
        d = (np.hypot(u / a, v / b) - 1) * min(a, b)
    elif kind == "rect":
        e1 = rng.uniform(0.3, 1.0)
        a, b = r, r * e1
        d = np.maximum(np.abs(u) - a, np.abs(v) - b)
    elif kind == "ring":
        e1 = rng.uniform(2, 6)
        d = np.abs(np.hypot(u, v) - r) - e1
    elif kind == "line":
        e1, e2 = rng.uniform(0.8, 3.0), rng.uniform(60, 300)
        d = np.maximum(np.abs(v) - e1, np.abs(u) - e2)
    else:
        e1 = rng.uniform(5, 18)
        d = np.hypot(u, v) - r
    m = soft(d)
    if kind == "texture":
        m = m * np.sin(2 * np.pi * u / e1)
    amp = rng.normal(0, 1.1)
    opaque = kind != "texture" and rng.uniform() < 0.5
    if opaque:
        level = rng.normal(0, 1.2) + 0.25 * grf(3.0)
        jump = np.abs(level - img)
        img = img * (1 - m) + level * m
        tmpl *= 1 - m
    else:
        jump = np.full((H, W), abs(amp))
        img += amp * m
        tmpl += amp * m
    near, inside = np.abs(d) < 1.5, d <= -1.5
    if kind == "texture":
        cat[d < 0.5] = 5
    elif kind in ("line", "ring"):
        z = d < 1.5
        cat[z], step[z] = (4 if kind == "ring" else 3), jump[z]
    else:
        if opaque:
            cat[inside] = 1
        cat[near], step[near] = 2, jump[near]
    objs.append([KINDS.index(kind), cx, cy, r, th, e1, e2, float(opaque), amp])

lo, hi = np.percentile(img, [0.5, 99.5])
scale = 16.0 / (hi - lo)
q = (img - lo) * scale + rng.normal(0, 0.28, size=img.shape)
pix = np.clip(np.floor(q), 0, 15).astype(np.uint8)
packed = (pix[:, 0::2] << 4) | pix[:, 1::2]
open(out + ".bin", "wb").write(packed.tobytes())
(tmpl * scale).astype(np.float32).tofile(out + ".tmpl.f32")
np.savez_compressed(out + ".npz", pix=pix, cat=cat, contrast=(step * scale).astype(np.float32),
                    tmpl=(tmpl * scale).astype(np.float32),
                    scale=scale, lo=lo, hi=hi, objects=np.array(objs))
share = np.bincount(cat.ravel(), minlength=6) / cat.size
print("seed %d scale %.3f levels/unit, %d objects, label shares %s" %
      (seed, scale, len(objs), " ".join("%.3f" % x for x in share)))
