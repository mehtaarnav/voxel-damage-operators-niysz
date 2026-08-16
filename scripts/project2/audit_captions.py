"""Audit figure captions: does each one define its own dispersion marks?

A caption must be readable alone. Where a panel shows spread -- error bars,
scattered points, a class-mean rule -- the caption has to say what the mark
means, or the reader is guessing.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEX = os.path.join(ROOT, "out", "writeup", "manuscript.tex")

s = open(TEX, encoding="utf8").read()
pat = re.compile(r"\\caption\{(.*?)\}\s*\n\s*\\label\{fig:(\w+)\}", re.S)
defines = re.compile(r"error bar|standard deviation|one s\.d\.|individual "
                     r"bisection|horizontal bar|class mean|points are",
                     re.I)

print(f"{'figure':<16}{'defines spread':>16}")
print("-" * 34)
for m in pat.finditer(s):
    cap = re.sub(r"\s+", " ", m.group(1))
    lab = m.group(2)
    print(f"{lab:<16}{str(bool(defines.search(cap))):>16}")
    print(f"    {cap[:160]}")
    print()
