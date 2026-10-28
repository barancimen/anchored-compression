"""Conservative C minifier: drop comments, collapse whitespace where tokens cannot merge.
Preprocessor lines keep their own line and internal spacing (collapsed to one space)."""
import re

_WORD = re.compile(r"[A-Za-z0-9_]")
_OPS = set("+-*/%&|^<>=!.:#")


def _strip_comments(s):
    out, i, n = [], 0, len(s)
    while i < n:
        c = s[i]
        if c in "\"'":
            j = i + 1
            while s[j] != c:
                j += 2 if s[j] == "\\" else 1
            out.append(s[i:j + 1])
            i = j + 1
        elif s.startswith("//", i):
            while i < n and s[i] != "\n":
                i += 1
        elif s.startswith("/*", i):
            j = s.index("*/", i + 2)
            out.append("\n" if "\n" in s[i:j] else " ")
            i = j + 2
        else:
            out.append(c)
            i += 1
    return "".join(out)


def _squeeze(line):
    """Collapse whitespace outside literals; keep one space only where needed."""
    toks, i, n = [], 0, len(line)
    while i < n:
        c = line[i]
        if c in "\"'":
            j = i + 1
            while line[j] != c:
                j += 2 if line[j] == "\\" else 1
            toks.append(line[i:j + 1])
            i = j + 1
        elif c.isspace():
            while i < n and line[i].isspace():
                i += 1
            toks.append(" ")
        else:
            toks.append(c)
            i += 1
    out = []
    for k, t in enumerate(toks):
        if t == " ":
            a = out[-1][-1] if out else ""
            b = toks[k + 1][0] if k + 1 < len(toks) else ""
            if a and b and ((_WORD.match(a) and _WORD.match(b)) or (a in _OPS and b in _OPS)):
                out.append(" ")
        else:
            out.append(t)
    return "".join(out)


def minify_c(src):
    src = _strip_comments(src)
    out, body = [], []
    for line in src.split("\n"):
        if line.strip().startswith("#"):
            if body:
                out.append(_squeeze(" ".join(body)))
                body = []
            out.append(re.sub(r"\s+", " ", line.strip()))
        elif line.strip():
            body.append(line)
    if body:
        out.append(_squeeze(" ".join(body)))
    return "\n".join(x for x in out if x) + "\n"


if __name__ == "__main__":
    import sys
    sys.stdout.write(minify_c(open(sys.argv[1]).read()))
