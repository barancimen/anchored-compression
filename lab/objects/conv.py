#!/usr/bin/env python3
"""Convert gen_layers.py's OUT.objs.txt (generator convention) to render.c's convention (lab tool).
usage: conv.py OBJS.txt > objects.txt   (kinds: 0 ellipse/disc, 1 rect/line, 2 ring, 3 texture)"""
import sys

for line in open(sys.argv[1]):
    k, cx, cy, r, th, e1, e2, q, a = line.split()
    k, cx, cy, r, th, e1, e2, q, a = int(k), float(cx), float(cy), float(r), float(th), float(e1), float(e2), int(q), float(a)
    if k == 0:
        K, R, E = 0, r, r
    elif k == 1:
        K, R, E = 0, r, r * e1
    elif k == 2:
        K, R, E = 1, r, r * e1
    elif k == 3:
        K, R, E = 2, r, e1
    elif k == 4:
        K, R, E = 1, e2, e1
    else:
        K, R, E = 3, r, e1
    print("%d %.4f %.4f %.4f %.6f %.4f %.4f %d" % (K, cx, cy, R, th, E, a, q))
