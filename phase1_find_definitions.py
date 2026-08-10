"""PHASE 1f — locate the papers' own DEFINITIONS of the quantities we must match.

We need to know exactly what the authors mean by:
  * P, the "percolation factor"  (is it connected-to-both-faces, or
    connected-to-one-inlet-plane?  This changes what we must compute.)
  * TPB "total" vs "active"
  * how TPB length is measured (voxel edges? centroid/skeleton?)
Because our Phase-2/4/5 gates compare against these numbers, a mismatch in
DEFINITION would look like a mismatch in RESULT.
"""

from __future__ import annotations

import io
import os
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "refs")
OUT = os.path.join(HERE, "out", "phase1")
os.makedirs(OUT, exist_ok=True)

DOCS = [
    ("TRANSPORT ma8095265", "ma8095265_transport_PMC5512617.epmc.xml"),
    ("TPB       ma8105370", "ma8105370_TPB_PMC5455394.epmc.xml"),
]

KEYWORDS = [
    "percolation factor",
    "percolat",
    "connectivity check",
    "active TPB",
    "TPB length",
    "centroid",
    "segmentation",
    "threshold",
]

WINDOW = 500


def plain(path):
    x = open(path, encoding="utf-8", errors="replace").read()
    t = re.sub(r"<[^>]+>", " ", x)
    return re.sub(r"\s+", " ", t)


def main():
    lines = []
    for tag, fname in DOCS:
        p = os.path.join(REFS, fname)
        if not os.path.exists(p):
            print("MISSING", p)
            continue
        t = plain(p)
        lines.append("=" * 78)
        lines.append(tag)
        lines.append("=" * 78)
        seen = set()
        for kw in KEYWORDS:
            hits = list(re.finditer(kw, t, re.I))
            if not hits:
                continue
            lines.append(f"\n### keyword: {kw!r}  ({len(hits)} hits)")
            for m in hits[:6]:
                a = max(0, m.start() - WINDOW)
                b = min(len(t), m.start() + WINDOW)
                key = t[a:b][:120]
                if key in seen:
                    continue
                seen.add(key)
                lines.append("   ..." + t[a:b] + "...")
                lines.append("")
    out = "\n".join(lines)
    print(out)
    dest = os.path.join(OUT, "phase1_definitions_context.txt")
    with open(dest, "w", encoding="utf-8") as f:
        f.write(out)
    print(f"\n[saved] {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
