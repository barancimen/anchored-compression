#!/usr/bin/env python3
"""Ideal cost of coding orders on the generator's background (lab calculator, no coder runs).

The background of gen_image.py is a stationary Gaussian field: three power-law random fields
(beta 3.4, 2.2, 1.0 with weights 1, .3, .08), times `scale` grey levels per field unit, plus
white sensor noise (sd .28 levels) and the rounding to integers (variance 1/12, treated as
white noise). Inside opaque objects the field is .25 x a beta 3.0 field plus the same noise.
For such a field the best predictor from any set of known pixels is linear, and its error
variance follows from the spectrum. A pixel then costs about 0.5 log2(2 pi e mse) bits (for
pure noise this gives 0.73 bpp, the entropy of the sensor noise after rounding).

So for any coding order (which pixels are known when a pixel is coded) and any stencil radius
this gives the cost an ideal, fully trained linear model would pay: a model-free way to
compare anchor schemes. It knows nothing about objects or about learning costs.

usage: gauss_orders.py [SCALE=2.5]
"""
import sys
from functools import lru_cache

import numpy as np

N = 500
SCALE = float(sys.argv[1]) if __name__ == "__main__" and len(sys.argv) > 1 else 2.5
WHITE = 0.28 ** 2 + 1 / 12
fy = np.fft.fftfreq(N)[:, None]
fx = np.fft.fftfreq(N)[None, :]
F = np.hypot(fx, fy)


def spec(beta):
    s = np.where(F > 0, np.maximum(F, 1e-12) ** -beta, 0.0)
    return s / s.mean()  # unit variance


def spectra(scale):
    """Spectra (grey levels^2) of the background and of an opaque object's inside."""
    bg = scale ** 2 * (spec(3.4) + 0.09 * spec(2.2) + 0.0064 * spec(1.0)) + WHITE
    return bg, scale ** 2 * 0.0625 * spec(3.0) + WHITE


S_BG, S_IN = spectra(SCALE)
R_BG = np.fft.ifft2(S_BG).real
R_IN = np.fft.ifft2(S_IN).real


def bits(mse):
    return 0.5 * np.log2(2 * np.pi * np.e * mse)


def szego(S):
    return bits(np.exp(np.mean(np.log(S))))


def mse_of(R, offs):
    if not offs:
        return R[0, 0]
    o = np.array(offs)
    d = o[:, None, :] - o[None, :, :]
    A = R[d[..., 0] % N, d[..., 1] % N]
    r = R[o[:, 0] % N, o[:, 1] % N]
    w = np.linalg.solve(A + 1e-9 * np.eye(len(o)), r)
    return R[0, 0] - r @ w


def disc(rad):
    k = int(rad)
    return [(dy, dx) for dy in range(-k, k + 1) for dx in range(-k, k + 1)
            if 0 < dy * dy + dx * dx <= rad * rad]


def causal(dy, dx):
    return dy < 0 or (dy == 0 and dx < 0)


def schedule_cost(R, period, pass_of, radius):
    """Mean bits per pixel of a schedule. pass_of(y, x) -> pass index (each pass is coded in
    raster order); radius(pass) -> stencil radius in pixels. Returns (bpp, per-pass bpp)."""
    cache, per = {}, {}
    for y in range(period):
        for x in range(period):
            p = pass_of(y, x)
            st = tuple(d for d in disc(radius(p))
                       if (q := pass_of(y + d[0], x + d[1])) < p or (q == p and causal(*d)))
            if st not in cache:
                cache[st] = bits(mse_of(R, list(st)))
            per.setdefault(p, []).append(cache[st])
    allb = [b for v in per.values() for b in v]
    return np.mean(allb), {p: (np.mean(v), len(v) / len(allb)) for p, v in sorted(per.items())}


def quadtree(T):
    """005-style hierarchy: lattice of step T first, then for t = T..2 centres, then edges."""
    levels = []
    t = T
    while t > 1:
        levels.append(t)
        t //= 2

    def pass_of(y, x):
        y, x = y % T, x % T
        if y == 0 and x == 0:
            return 0
        for i, t in enumerate(levels):
            h = t // 2
            if y % t == h and x % t == h:
                return 1 + 2 * i
            if (y % t == h and x % t == 0) or (y % t == 0 and x % t == h):
                return 2 + 2 * i
        raise AssertionError
    return pass_of


def lattice_then_raster(T):
    return lambda y, x: 0 if (y % T == 0 and x % T == 0) else 1


if __name__ == "__main__":
    print("scale %.2f levels/unit; field variance %.2f levels^2; white noise %.3f" %
          (SCALE, R_BG[0, 0] - WHITE, WHITE))
    print("\nideal raster limit (Szego, whole causal past): background %.4f bpp, opaque interior %.4f bpp"
          % (szego(S_BG), szego(S_IN)))
    print("noise alone: %.4f bpp" % bits(WHITE))
    every = [d for d in disc(15)]
    print("all other pixels known within radius 15 (interpolation floor): background %.4f bpp"
          % bits(mse_of(R_BG, every)))

    print("\nraster, causal stencil of radius r (taps): background bpp / interior bpp")
    for r in (1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15):
        st = [d for d in disc(r) if causal(*d)]
        print("  r=%-4g (%3d taps)  %.4f  %.4f" % (r, len(st), bits(mse_of(R_BG, st)), bits(mse_of(R_IN, st))))

    print("\nquadtree (005 style), anchors on the lattice of step T, radius in px (anchor pass radius in brackets)")
    for T, ra, rf in ((2, 12, 6), (2, 14, 7), (4, 24, 6), (8, 48, 6), (2, 12, 3), (2, 12, 2)):
        b, per = schedule_cost(R_BG, T, quadtree(T), lambda p: ra if p == 0 else rf)
        parts = "  ".join("p%d %.3f x%.3f" % (p, v, w) for p, (v, w) in per.items())
        print("  T=%d r=%g (%g): %.4f bpp   [%s]" % (T, rf, ra, b, parts))

    print("\nlattice anchors of step T, then raster with look-ahead onto the anchors, radius r")
    for T, r in ((16, 7), (8, 7), (4, 7), (4, 4)):
        b, per = schedule_cost(R_BG, T, lattice_then_raster(T), lambda p: 4 * T if p == 0 else r)
        print("  T=%d r=%g: %.4f bpp  (anchors %.3f bpp, rest %.4f bpp)" % (T, r, b, per[0][0], per[1][0]))
