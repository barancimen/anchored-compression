#!/usr/bin/env python3
"""Oracle layers for lab/objects/oracle.c (lab tool). Replays lab/gen_labeled.py's random draws
(same seed = same image) and writes float32 layers in grey levels, row-major:
  OUT.add.f32       additive objects (= gen_labeled's tmpl)
  OUT.addK.f32      only the additive objects of kind K (0 disc .. 5 texture); opaque ones still cover
  OUT.full.f32      additive objects plus every opaque object treated as additive, amplitude =
                    the mask-weighted mean of (level - image under it); it covers what is below
  OUT.opq.f32       only the opaque objects treated that way
  OUT.objs.txt      one line per object: kind cx cy r th e1 e2 opaque amp(grey levels, effective)
usage: gen_layers.py SEED OUT
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
add = np.zeros((H, W))
addk = [np.zeros((H, W)) for _ in KINDS]
full = np.zeros((H, W))
opq = np.zeros((H, W))
objs = []
nobj = int(rng.integers(30, 45))
for _ in range(nobj):
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
    ki = KINDS.index(kind)
    if opaque:
        level = rng.normal(0, 1.2) + 0.25 * grf(3.0)
        ae = float((m * (level - img)).sum() / max(m.sum(), 1e-9))
        img = img * (1 - m) + level * m
        add *= 1 - m
        for t in addk:
            t *= 1 - m
        full = full * (1 - m) + ae * m
        opq = opq * (1 - m) + ae * m
        amp = ae
    else:
        img += amp * m
        add += amp * m
        addk[ki] += amp * m
        full += amp * m
        opq *= 1
    objs.append([ki, cx, cy, r, th, e1, e2, float(opaque), amp])

lo, hi = np.percentile(img, [0.5, 99.5])
scale = 16.0 / (hi - lo)
q = (img - lo) * scale + rng.normal(0, 0.28, size=img.shape)
pix = np.clip(np.floor(q), 0, 15).astype(np.uint8)
packed = (pix[:, 0::2] << 4) | pix[:, 1::2]
open(out + ".bin", "wb").write(packed.tobytes())
(add * scale).astype(np.float32).tofile(out + ".add.f32")
for i, t in enumerate(addk):
    (t * scale).astype(np.float32).tofile(out + ".add%d.f32" % i)
(full * scale).astype(np.float32).tofile(out + ".full.f32")
(opq * scale).astype(np.float32).tofile(out + ".opq.f32")
with open(out + ".objs.txt", "w") as f:
    for o in objs:
        f.write("%d %.3f %.3f %.3f %.4f %.3f %.3f %d %.3f\n" % (*o[:8], o[8] * scale))
print("seed %d scale %.3f lo %.4f, %d objects" % (seed, scale, lo, len(objs)))
