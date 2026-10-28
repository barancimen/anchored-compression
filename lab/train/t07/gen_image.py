#!/usr/bin/env python3
"""Generate the benchmark image: 500x500 pixels, 4 bits/pixel (values 0..15).

Content: a multi-scale Gaussian random field (spatially autocorrelated "terrain"),
plus anti-aliased objects (discs, ellipses, rotated rectangles, rings, thick lines),
a few patches of oriented texture, a gentle illumination gradient, and mild
sensor noise, then quantized to 16 levels.

All randomness is drawn from os.urandom and the seed is never stored. This matters
for a Kolmogorov-style contest: if the seed were known, "generator + seed" would be
an unbeatable decompressor. Without it, a compressor has to model the pixels.

Output (canonical benchmark file, the analogue of enwik9):
  image.bin : 125000 bytes, row-major, 2 pixels/byte, left pixel in the HIGH nibble
Viewing aid (not part of the benchmark):
  image.png : 8-bit grayscale (value*17)
"""
import os
import sys

import numpy as np
from PIL import Image

H = W = 500
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "image.bin")
if os.path.exists(OUT) and "--force" not in sys.argv:
    sys.exit("image.bin already exists; refusing to overwrite the benchmark (use --force)")

rng = np.random.Generator(np.random.PCG64(int.from_bytes(os.urandom(32), "little")))
yy, xx = np.mgrid[0:H, 0:W].astype(np.float64)


def grf(beta):
    """Unit-variance Gaussian random field with power spectrum ~ 1/f^beta."""
    fy = np.fft.fftfreq(H)[:, None]
    fx = np.fft.rfftfreq(W)[None, :]
    f = np.hypot(fx, fy)
    f[0, 0] = 1.0
    spec = (rng.normal(size=f.shape) + 1j * rng.normal(size=f.shape)) * f ** (-beta / 2)
    spec[0, 0] = 0
    z = np.fft.irfft2(spec, s=(H, W))
    return (z - z.mean()) / z.std()


def soft(d):
    """Anti-aliased inside-mask from a signed distance d (negative = inside)."""
    return np.clip(0.5 - d, 0.0, 1.0)


img = 1.00 * grf(3.4) + 0.30 * grf(2.2) + 0.08 * grf(1.0)
img += rng.normal(0, 0.6) * (xx - W / 2) / W + rng.normal(0, 0.6) * (yy - H / 2) / H

for _ in range(int(rng.integers(30, 45))):
    kind = rng.choice(["disc", "ellipse", "rect", "ring", "line", "texture"],
                      p=[0.22, 0.18, 0.22, 0.10, 0.14, 0.14])
    cx, cy = rng.uniform(0, W), rng.uniform(0, H)
    r = rng.uniform(8, 60)
    th = rng.uniform(0, np.pi)
    c, s = np.cos(th), np.sin(th)
    u = (xx - cx) * c + (yy - cy) * s          # rotated coordinates
    v = -(xx - cx) * s + (yy - cy) * c
    if kind == "disc":
        m = soft(np.hypot(u, v) - r)
    elif kind == "ellipse":
        a, b = r, r * rng.uniform(0.3, 0.9)
        k = np.hypot(u / a, v / b)
        m = soft((k - 1) * min(a, b))
    elif kind == "rect":
        a, b = r, r * rng.uniform(0.3, 1.0)
        m = soft(np.maximum(np.abs(u) - a, np.abs(v) - b))
    elif kind == "ring":
        w = rng.uniform(2, 6)
        m = soft(np.abs(np.hypot(u, v) - r) - w)
    elif kind == "line":
        w = rng.uniform(0.8, 3.0)
        L = rng.uniform(60, 300)
        m = soft(np.maximum(np.abs(v) - w, np.abs(u) - L))
    else:  # oriented sinusoidal texture inside a disc
        period = rng.uniform(5, 18)
        m = soft(np.hypot(u, v) - r) * np.sin(2 * np.pi * u / period)
    amp = rng.normal(0, 1.1)
    if kind != "texture" and rng.uniform() < 0.5:
        # opaque object: replaces the background with its own lightly shaded level
        level = rng.normal(0, 1.2) + 0.25 * grf(3.0)
        img = img * (1 - m) + level * m
    else:
        img += amp * m

# Map to 16 levels: 0.5/99.5 percentiles -> [0,16), add mild sensor noise, quantize.
lo, hi = np.percentile(img, [0.5, 99.5])
q = (img - lo) / (hi - lo) * 16.0 + rng.normal(0, 0.28, size=img.shape)
pix = np.clip(np.floor(q), 0, 15).astype(np.uint8)

packed = (pix[:, 0::2] << 4) | pix[:, 1::2]
with open(OUT, "wb") as f:
    f.write(packed.tobytes())
Image.fromarray(pix * 17).save(OUT[:-4] + ".png")
print("wrote", OUT, packed.size, "bytes; histogram:", np.bincount(pix.ravel(), minlength=16).tolist())
