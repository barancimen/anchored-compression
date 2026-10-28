#!/usr/bin/env python3
"""Hutter-Prize-style scorer for compressors of image.bin.  THIS FILE IS FIXED: the
research loop may not change it (it is the analogue of the contest rules).

Rules implemented (http://prize.hutter1.net/hrules.htm, hfaq.htm):
  * Programs are Linux x86 executables, or zip files of C/C++/Python/asm source
    (+ makefile).  L(X) = byte length of the executable or of the zip.
  * Forms and total size S (+ "charges", see below):
      sfx      : comp C outputs a self-extracting executable DA.   S = L(C) + L(DA)
      same     : one program P is both comp and decomp (C = D).     S = L(P) + L(P) + L(A)
      separate : comp C and decomp D submitted separately.          S = L(C) + 2 L(D) + L(A)
  * "If command-line options for execution or compilation are necessary, their length
    is added to S."  We are stricter: EVERY character of EVERY command line we execute
    (build lines, the compress line including its input redirection, the decompress
    line including its output redirection) is charged.
  * Decompression must use no outside information: it runs in a fresh directory that
    holds only the decompressor (rebuilt from its zip) and the archive, with no network,
    and with this project directory (where the original lives) hidden behind an empty
    tmpfs.  The output must be byte-identical to image.bin.
  * The input file is always called image.bin (like enwik9); all other names are free.
  * Resource limits: single core, <=10 GB address space, and a 30 min wall clock per
    step (far inside the contest's ~50 h budget).

Usage: evaluate.py CANDIDATE_DIR [-m "description"] [--ref]
  --ref: score a reference point (status "ref"); it never becomes best.json.
Work files go to $HP_SCRATCH (default: <tmp>/anchored-compression-eval), which must lie
outside this project directory because decompression runs with the project hidden.
CANDIDATE_DIR/cand.json, e.g. for the C = D form:
  {"form": "same",
   "P": {"files": ["p.c"], "minify": "c", "build": ["make p"]},
   "comp": "./p<image.bin>a", "archive": ["a"],
   "decomp": "./p<a>o", "output": "o"}
"""
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile, time, resource, zlib, struct
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IMAGE = os.path.join(ROOT, "image.bin")
RESULTS = os.path.join(ROOT, "results.tsv")
BEST = os.path.join(ROOT, "best.json")
SCRATCH = os.path.abspath(os.environ.get("HP_SCRATCH", os.path.join(tempfile.gettempdir(), "anchored-compression-eval")))
assert not (SCRATCH + os.sep).startswith(ROOT + os.sep), "HP_SCRATCH must be outside the project"
TIMEOUT = 1800
MEM = 10 * 10**9
NPIX = 500 * 500
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import minify  # noqa: E402

try:
    import zopfli.zlib as _zopfli
except ImportError:
    _zopfli = None


def deflate_raw(data):
    """Best raw-deflate stream we can make (plain deflate, method 8: any unzip reads it)."""
    cands = []
    c = zlib.compressobj(9, zlib.DEFLATED, -15, 9)
    cands.append(c.compress(data) + c.flush())
    if _zopfli is not None:
        z = _zopfli.compress(data, numiterations=1000)  # zlib-wrapped: strip 2-byte header + adler32
        cands.append(z[2:-4])
    return min(cands, key=len)


def make_zip(entries):
    """Minimal standard zip (no extra fields, no comments, fixed timestamps)."""
    out, central = bytearray(), bytearray()
    for name, data in entries:
        nb = name.encode()
        crc = zlib.crc32(data) & 0xFFFFFFFF
        comp = deflate_raw(data)
        method, payload = (8, comp) if len(comp) < len(data) else (0, data)
        off = len(out)
        out += struct.pack("<IHHHHHIIIHH", 0x04034B50, 20, 0, method, 0, 0x21, crc,
                           len(payload), len(data), len(nb), 0) + nb + payload
        central += struct.pack("<IHHHHHHIIIHHHHHII", 0x02014B50, 20, 20, 0, method, 0, 0x21, crc,
                               len(payload), len(data), len(nb), 0, 0, 0, 0, 0, off) + nb
    cd_off = len(out)
    out += central
    out += struct.pack("<IHHHHIIH", 0x06054B50, 0, 0, len(entries), len(entries), len(central), cd_off, 0)
    return bytes(out)


def package(cdir, spec, outdir, tag):
    """Return (path, length) of the measured artefact for one program."""
    if "exe" in spec:  # a ready executable is measured as-is
        dst = os.path.join(outdir, tag)
        shutil.copy(os.path.join(cdir, spec["exe"]), dst)
        os.chmod(dst, 0o755)
        return dst, os.path.getsize(dst)
    entries = []
    for f in spec["files"]:
        data = open(os.path.join(cdir, f), "rb").read()
        if spec.get("minify") == "c" and f.endswith((".c", ".h", ".cpp")):
            data = minify.minify_c(data.decode()).encode()
        entries.append((os.path.basename(f), data))
    z = make_zip(entries)
    dst = os.path.join(outdir, tag + ".zip")
    open(dst, "wb").write(z)
    with zipfile.ZipFile(dst) as zf:  # must be a valid zip that round-trips
        assert zf.testzip() is None
        for name, data in entries:
            assert zf.read(name) == data
    return dst, len(z)


def install(artefact, d):
    if artefact.endswith(".zip"):
        subprocess.run(["unzip", "-q", artefact, "-d", d], check=True)
    else:
        shutil.copy(artefact, d)


def run(cmd, cwd, hide_project):
    """Run one command line in a sandbox (no network; optionally project dir hidden)."""
    inner = f"cd {cwd} && exec prlimit --as={MEM} timeout {TIMEOUT} sh -c '{cmd}'"
    if hide_project:
        inner = f"mount -t tmpfs -o size=1m tmpfs {ROOT} && " + inner
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.perf_counter()
    p = subprocess.run(["unshare", "-n", "-m", "sh", "-c", inner], capture_output=True)
    dt = time.perf_counter() - t0
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (r1.ru_utime + r1.ru_stime) - (r0.ru_utime + r0.ru_stime)
    if p.returncode != 0:
        raise RuntimeError(f"command failed ({p.returncode}): {cmd}\n{p.stderr.decode()[-2000:]}")
    return dt, cpu, r1.ru_maxrss / 1024


def evaluate(cdir):
    cand = json.load(open(os.path.join(cdir, "cand.json")))
    form = cand["form"]
    os.makedirs(SCRATCH, exist_ok=True)
    work = tempfile.mkdtemp(prefix="ev-", dir=SCRATCH)
    sub = os.path.join(cdir, "submission")
    shutil.rmtree(sub, ignore_errors=True)
    os.makedirs(sub)
    charges, info = 0, {}
    try:
        if form == "same":
            P, LP = package(cdir, cand["P"], sub, "P")
            progC = progD = (P, cand["P"])
            code = 2 * LP
            info["L_P"] = LP
        elif form == "separate":
            C, LC = package(cdir, cand["C"], sub, "C")
            D, LD = package(cdir, cand["D"], sub, "D")
            progC, progD = (C, cand["C"]), (D, cand["D"])
            code = LC + 2 * LD
            info.update(L_C=LC, L_D=LD)
        elif form == "sfx":
            C, LC = package(cdir, cand["C"], sub, "C")
            progC, progD = (C, cand["C"]), None
            code = LC
            info["L_C"] = LC
        else:
            raise ValueError(form)

        # ---- compress: program + image.bin only
        dc = os.path.join(work, "c")
        os.makedirs(dc)
        install(progC[0], dc)
        shutil.copy(IMAGE, os.path.join(dc, "image.bin"))
        for b in progC[1].get("build", []):
            run(b, dc, False)
            charges += len(b)
        t_c, cpu_c, rss_c = run(cand["comp"], dc, False)
        charges += len(cand["comp"])
        arch = os.path.join(work, "archive")
        os.makedirs(arch)
        LA = 0
        for a in cand["archive"]:
            shutil.copy(os.path.join(dc, a), os.path.join(arch, a))
            LA += os.path.getsize(os.path.join(dc, a))
        shutil.rmtree(dc)  # nothing from the compress run survives except the archive

        # ---- decompress: fresh dir, rebuilt decompressor + archive, project hidden
        dd = os.path.join(work, "d")
        os.makedirs(dd)
        if progD is not None:
            install(progD[0], dd)
        for a in cand["archive"]:
            shutil.copy(os.path.join(arch, a), os.path.join(dd, a))
            if form == "sfx":
                os.chmod(os.path.join(dd, a), 0o644)  # as received; chmod must be in the command
        if progD is not None:
            for b in progD[1].get("build", []):
                run(b, dd, True)
                charges += len(b)
        t_d, cpu_d, rss_d = run(cand["decomp"], dd, True)
        charges += len(cand["decomp"])
        out = open(os.path.join(dd, cand["output"]), "rb").read()
        ok = hashlib.sha256(out).digest() == hashlib.sha256(open(IMAGE, "rb").read()).digest()
        S = code + LA + charges
        info.update(S=S, code=code, L_A=LA, charges=charges, bpp=LA * 8 / NPIX,
                    t_comp=t_c, t_dec=t_d, cpu_comp=cpu_c, cpu_dec=cpu_d,
                    rss_mb=max(rss_c, rss_d), verified=ok)
        return info
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cand")
    ap.add_argument("-m", "--desc", default="")
    ap.add_argument("--ref", action="store_true", help="reference point: never becomes best")
    args = ap.parse_args()
    cdir = os.path.abspath(args.cand)
    name = os.path.basename(cdir)
    best = json.load(open(BEST)) if os.path.exists(BEST) else {"S": float("inf"), "exp": None}
    try:
        info = evaluate(cdir)
        if not info["verified"]:
            status = "mismatch"
        elif args.ref:
            status = "ref"
        elif info["S"] < best["S"]:
            status = "keep"
        else:
            status = "discard"
    except Exception as e:  # noqa: BLE001
        info, status = {}, "crash"
        print(e, file=sys.stderr)
    new = not os.path.exists(RESULTS)
    with open(RESULTS, "a") as f:
        if new:
            f.write("exp\tstatus\tS\tcode\tL_A\tcharges\tbpp\tt_comp\tt_dec\trss_mb\tdescription\n")
        g = lambda k, fmt: (fmt % info[k]) if k in info else "-"  # noqa: E731
        f.write("\t".join([name, status, g("S", "%d"), g("code", "%d"), g("L_A", "%d"), g("charges", "%d"),
                           g("bpp", "%.4f"), g("t_comp", "%.1f"), g("t_dec", "%.1f"), g("rss_mb", "%.0f"),
                           args.desc]) + "\n")
    if status == "keep":
        json.dump({"exp": name, **info}, open(BEST, "w"), indent=1)
    print(json.dumps({"exp": name, "status": status, **info}, indent=1))


if __name__ == "__main__":
    main()
