"""Measure our prose against the prose of the papers we cite.

The complaint is that the manuscript reads as machine-written. Rewriting it
against my own taste would risk reproducing the same voice, so this builds the
target from the corpus instead: the same measurements, run on the full text of
papers from the venues this work is aimed at.

Sources are whatever is present in refs/ -- PubMed Central XML for the Pecho
papers, PDFs for the rest.
"""
import os
import re
import sys
import glob

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PATTERNS = {
    "not-X-but-Y":       r"\bis not\b[^.]{0,80}\bbut\b|\bnot\b[^.]{0,60}\brather\b",
    "That/This opener":  r"^(That|This) is\b",
    "colon-expansion":   r"^[^:]{15,90}:\s+[a-z]",
    "em-dash":           r"—|---| - ",
    "short (<=8 words)": None,
}


def sentences(text):
    t = re.sub(r"\s+", " ", text)
    return [x.strip() for x in re.split(r"(?<=[.!?])\s+", t)
            if 15 < len(x.strip()) < 600]


def from_xml(path):
    raw = open(path, encoding="utf8", errors="replace").read()
    body = re.search(r"<body>(.*?)</body>", raw, re.S)
    if not body:
        return ""
    t = body.group(1)
    t = re.sub(r"<(table-wrap|fig|xref|inline-formula|disp-formula)[^>]*>.*?</\1>",
               " ", t, flags=re.S)
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", t)


def from_pdf(path):
    try:
        from pypdf import PdfReader
    except ImportError:
        try:
            from PyPDF2 import PdfReader
        except ImportError:
            return ""
    try:
        r = PdfReader(path)
        return " ".join((p.extract_text() or "") for p in r.pages)
    except Exception:
        return ""


def profile(name, text):
    sents = sentences(text)
    if len(sents) < 40:
        return None
    n = len(sents)
    out = {"n": n, "mean words": sum(len(s.split()) for s in sents) / n}
    for k, p in PATTERNS.items():
        if p is None:
            hits = sum(1 for s in sents if len(s.split()) <= 8)
        else:
            hits = sum(1 for s in sents if re.search(p, s, re.I))
        out[k] = hits / n * 100
    return name, out


def main():
    rows = []
    refs = os.path.join(ROOT, "refs")
    for p in sorted(glob.glob(os.path.join(refs, "*.xml"))):
        if "efetch" not in p:
            continue
        r = profile(os.path.basename(p)[:28], from_xml(p))
        if r:
            rows.append(r)
    for p in sorted(glob.glob(os.path.join(refs, "*.pdf"))):
        r = profile(os.path.basename(p)[:28], from_pdf(p))
        if r:
            rows.append(r)

    tex = open(os.path.join(ROOT, "out", "writeup", "manuscript.tex"),
               encoding="utf8").read()
    body = tex.split(r"\section{Introduction}")[1].split(r"\begin{thebibliography}")[0]
    body = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?", " ", body)
    body = re.sub(r"[{}$&%]", " ", body)
    ours = profile("*** THIS MANUSCRIPT ***", body)
    if ours:
        rows.append(ours)

    if not rows:
        print("no usable sources found in refs/")
        return 1
    cols = ["n", "mean words"] + list(PATTERNS)
    w = max(len(r[0]) for r in rows) + 2
    print(f"{'source':<{w}}" + "".join(f"{c:>18}" for c in cols))
    print("-" * (w + 18 * len(cols)))
    for name, d in rows:
        line = f"{name:<{w}}"
        for c in cols:
            v = d[c]
            line += f"{v:>18.1f}" if isinstance(v, float) else f"{v:>18d}"
        print(line)
    print("\npercentages are share of sentences; 'mean words' is sentence length")
    return 0


if __name__ == "__main__":
    sys.exit(main())
