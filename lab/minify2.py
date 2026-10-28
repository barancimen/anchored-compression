"""Tighter C minifier for submissions (lab tool; the harness keeps its own conservative one).
Keeps a space only between two word chars or where the two chars could fuse into a different
C token (multi-char punctuators, comment openers, number/dot adjacency)."""
import re, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'harness'))
import minify as base
_W = re.compile(r"[A-Za-z0-9_]")
FUSE = {"->", "++", "--", "<<", ">>", "<=", ">=", "==", "!=", "&&", "||", "*=", "/=", "%=", "+=",
        "-=", "&=", "^=", "|=", "##", "<:", ":>", "<%", "%>", "%:", "/*", "//", ".."}
def need(a, b):
    if _W.match(a) and _W.match(b): return True
    if a + b in FUSE: return True
    if (a == "." and b.isdigit()) or (a.isdigit() and b == "."): return True
    if a in "+-" and b in "+-": return True
    return False
def squeeze(line):
    toks, i, n = [], 0, len(line)
    while i < n:
        c = line[i]
        if c in "\"'":
            j = i + 1
            while line[j] != c: j += 2 if line[j] == "\\" else 1
            toks.append(line[i:j + 1]); i = j + 1
        elif c.isspace():
            while i < n and line[i].isspace(): i += 1
            toks.append(" ")
        else:
            toks.append(c); i += 1
    out = []
    for k, t in enumerate(toks):
        if t == " ":
            a = out[-1][-1] if out else ""; b = toks[k + 1][0] if k + 1 < len(toks) else ""
            if a and b and need(a, b): out.append(" ")
        else: out.append(t)
    return "".join(out)
def minify_c(src):
    src = base._strip_comments(src)
    out, body = [], []
    for line in src.split("\n"):
        s = line.strip()
        if s.startswith("#"):
            if body: out.append(squeeze(" ".join(body))); body = []
            m = re.match(r"#\s*define\s+(\w+(?:\([^)]*\))?)\s*(.*)", s)
            out.append("#define " + m.group(1).replace(" ", "") + " " + squeeze(m.group(2)) if m else re.sub(r"\s+", " ", s).replace("#include <", "#include<"))
        elif s: body.append(line)
    if body: out.append(squeeze(" ".join(body)))
    return "\n".join(x for x in out if x) + "\n"
if __name__ == "__main__":
    sys.stdout.write(minify_c(open(sys.argv[1]).read()))
