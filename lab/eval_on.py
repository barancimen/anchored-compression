#!/usr/bin/env python3
"""Score candidates on a different image with the unchanged harness (same sandbox, same
Hutter accounting). Usage: eval_on.py IMAGE CAND_DIR [CAND_DIR ...]  -> lab/generalization.tsv"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'harness'))
import evaluate
img = os.path.abspath(sys.argv[1])
evaluate.IMAGE = img  # the only change: which file is image.bin
log = os.path.join(HERE, 'generalization.tsv')
new = not os.path.exists(log)
for c in sys.argv[2:]:
    try:
        info = evaluate.evaluate(os.path.abspath(c)); st = 'ok' if info['verified'] else 'MISMATCH'
    except Exception as e:
        info, st = {}, 'CRASH: %s' % str(e).splitlines()[0]
    row = [os.path.relpath(img, os.path.join(HERE, '..')), os.path.basename(c.rstrip('/')), st] + \
          [str(info.get(k, '-')) for k in ('S', 'code', 'L_A', 'charges')] + ['%.4f' % info['bpp'] if 'bpp' in info else '-']
    with open(log, 'a') as f:
        if new: f.write('image\tcand\tstatus\tS\tcode\tL_A\tcharges\tbpp\n'); new = False
        f.write('\t'.join(row) + '\n')
    print('\t'.join(row), flush=True)
