#!/usr/bin/env python3
"""Local search over the quadtree candidate's 5 counter contexts (lab/qt.c template).
Each context is a tuple of indices into the 12 nearest stencil pixels (or 'mq' = the
(prediction, error bucket) context). Moves: replace one index, add one, drop one.
Score = archive bytes of image.bin (local -O2 build). Usage: qtsearch.py [ROUNDS]"""
import itertools, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, 'build'); os.makedirs(B, exist_ok=True)
IMG = os.path.join(HERE, '..', 'image.bin')
cache = {}
def expr(t):
    if t == 'mq': return '(int)(m*2+.5)<<5|q'
    return '|'.join('r[%d]<<%d' % (i, 4 * (len(t) - 1 - j)) for j, i in enumerate(t))
def score(ctx):
    key = tuple(ctx)
    if key not in cache:
        defs = ['-DCX%d=%s' % (i, expr(t)) for i, t in enumerate(ctx)]
        exe = os.path.join(B, 'qts')
        subprocess.run(['cc', '-O2', '-o', exe, os.path.join(HERE, 'qt.c'), '-lm'] + defs, check=True,
                       stderr=subprocess.DEVNULL)
        cache[key] = len(subprocess.run([exe], stdin=open(IMG, 'rb'), capture_output=True).stdout)
    return cache[key]
cur = [(0, 1, 2, 3), (0, 1), (2, 3), 'mq', (0, 1, 4, 5)]
best = score(cur); print('start', cur, best, flush=True)
for rnd in range(int(sys.argv[1]) if len(sys.argv) > 1 else 2):
    improved = False
    for s in range(5):
        if cur[s] == 'mq': continue
        t = cur[s]; moves = []
        for j in range(len(t)):
            for v in range(12):
                if v not in t: moves.append(t[:j] + (v,) + t[j + 1:])
            if len(t) > 1: moves.append(t[:j] + t[j + 1:])
        if len(t) < 4: moves += [t + (v,) for v in range(12) if v not in t]
        for m in moves:
            c = cur[:s] + [m] + cur[s + 1:]
            sc = score(c)
            if sc < best:
                best, cur, improved = sc, c, True
                print('  slot %d -> %s: %d' % (s, m, sc), flush=True)
    print('round', rnd, cur, best, flush=True)
    if not improved: break
print('final', cur, best)
print(' '.join("-DCX%d='%s'" % (i, expr(t)) for i, t in enumerate(cur)))
