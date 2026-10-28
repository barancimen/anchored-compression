#!/usr/bin/env python3
"""Make the inpainting animation and the anchors-only picture.

Usage: gif.py IMAGE SAMPLEDIR [T=8]
Writes, next to this script:
  animation.gif  one bare 500x500 frame per sample (value v shown as grey 17v, exact
                 16-grey palette), 40 ms per frame = 25 fps, looping forever
  anchors.png    the input: the anchor pixels of IMAGE (x % T == 0 and y % T == 0) on black
Checks that every sample keeps the anchors of IMAGE.
"""
import glob
import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
PALETTE = [17 * v for v in range(16) for _ in range(3)]


def load(fn):
    raw = np.fromfile(fn, np.uint8)
    x = np.empty(2 * raw.size, np.uint8)
    x[0::2], x[1::2] = raw >> 4, raw & 15
    return x.reshape(500, 500)


def frame(x):
    im = Image.fromarray(x, "P")
    im.putpalette(PALETTE)
    return im


def main():
    image, sdir = sys.argv[1], sys.argv[2]
    t = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    orig = load(image)
    files = sorted(glob.glob(os.path.join(sdir, "sample_*.bin")))
    samples = [load(f) for f in files]
    for f, s in zip(files, samples):
        assert (s[::t, ::t] == orig[::t, ::t]).all(), f"{f}: anchors differ from {image}"
    frames = [frame(s) for s in samples]
    frames[0].save(os.path.join(HERE, "animation.gif"), save_all=True, append_images=frames[1:],
                   duration=40, loop=0, optimize=False)
    anchors = np.zeros_like(orig)
    anchors[::t, ::t] = orig[::t, ::t]
    Image.fromarray(anchors * 17, "L").save(os.path.join(HERE, "anchors.png"))
    print(f"{len(frames)} frames -> animation.gif ({len(frames) / 25:.2f} s loop), anchors.png")


if __name__ == "__main__":
    main()
