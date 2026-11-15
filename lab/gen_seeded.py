#!/usr/bin/env python3
"""Run the unchanged gen_image.py with a known 32-byte seed in place of os.urandom (lab tool).
Usage: gen_seeded.py SEEDHEX OUTDIR  ->  OUTDIR/image.bin, OUTDIR/image.png, OUTDIR/seed
The script is copied to OUTDIR and run from there, because it writes next to itself and
refuses to overwrite an existing image.bin."""
import os, runpy, shutil, sys
here = os.path.dirname(os.path.abspath(__file__))
seed = bytes.fromhex(sys.argv[1])
assert len(seed) == 32
out = os.path.abspath(sys.argv[2])
os.makedirs(out, exist_ok=True)
shutil.copy(os.path.join(here, '..', 'gen_image.py'), out)
open(os.path.join(out, 'seed'), 'wb').write(seed)
calls, real = [], os.urandom
def urandom(n):
    calls.append(n)
    return seed if n == 32 else real(n)
os.urandom = urandom
sys.argv = ['gen_image.py']
runpy.run_path(os.path.join(out, 'gen_image.py'), run_name='__main__')
assert calls.count(32) == 1, calls  # the generator's one draw of entropy got our seed
