#!/usr/bin/env python3
"""Check and measure a seed-replay candidate (016) on images whose seed is known.

image.bin's own seed was never stored, so harness/evaluate.py cannot score such a candidate.
This runs the unchanged harness on images that the unchanged gen_image.py makes from known
seeds (lab/gen_seeded.py), with one substitution: the compress command is not run. Its
archive is written down from the seed instead, because no program can find it from the image.
Packaging, build, sandbox, decompression with the project hidden, the byte-for-byte
comparison, S and the charges (the compress command line included) are the harness's own.

cand.json names what the archive holds in "seed_archive":
  "pcg64_state": the 16-byte state and 16-byte increment of numpy's PCG64 after seeding
  "urandom":     the 32 bytes gen_image.py read from os.urandom

Usage: seed_check.py CAND_DIR SEED [SEED ...]   -> appends to lab/seed_replay.tsv
SEED is 64 hex digits, or a number n standing for sha256("seed replay test n")."""
import hashlib, json, os, subprocess, sys, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'harness'))
import evaluate


def archive_from(kind, seed):
    if kind == 'urandom':
        return seed
    st = np.random.PCG64(int.from_bytes(seed, 'little')).state['state']
    return st['state'].to_bytes(16, 'little') + st['inc'].to_bytes(16, 'little')


cdir = os.path.abspath(sys.argv[1])
cand = json.load(open(os.path.join(cdir, 'cand.json')))
assert len(cand['archive']) == 1
log = os.path.join(HERE, 'seed_replay.tsv')
new = not os.path.exists(log)
os.makedirs(evaluate.SCRATCH, exist_ok=True)
harness_run = evaluate.run
for arg in sys.argv[2:]:
    seed = bytes.fromhex(arg) if len(arg) == 64 else hashlib.sha256(b'seed replay test ' + arg.encode()).digest()
    d = tempfile.mkdtemp(prefix='seed-', dir=evaluate.SCRATCH)
    subprocess.run([sys.executable, os.path.join(HERE, 'gen_seeded.py'), seed.hex(), d], check=True,
                   stdout=subprocess.DEVNULL)
    evaluate.IMAGE = os.path.join(d, 'image.bin')
    archive = archive_from(cand['seed_archive'], seed)

    def run(cmd, cwd, hide_project):
        if cmd == cand['comp']:  # the one substitution: the archive comes from the seed
            open(os.path.join(cwd, cand['archive'][0]), 'wb').write(archive)
            return 0.0, 0.0, 0.0
        return harness_run(cmd, cwd, hide_project)

    evaluate.run = run
    try:
        info = evaluate.evaluate(cdir)
        st = 'ok' if info['verified'] else 'MISMATCH'
    except Exception as e:
        info, st = {}, 'CRASH: %s' % str(e).splitlines()[0]
    row = [seed.hex()[:16], os.path.basename(cdir.rstrip('/')), st] + \
          [str(info.get(k, '-')) for k in ('S', 'code', 'L_A', 'charges')] + \
          ['%.1f' % info['t_dec'] if 't_dec' in info else '-']
    with open(log, 'a') as f:
        if new:
            f.write('seed\tcand\tstatus\tS\tcode\tL_A\tcharges\tt_dec\n')
            new = False
        f.write('\t'.join(row) + '\n')
    print('\t'.join(row), flush=True)
    subprocess.run(['rm', '-rf', d])
