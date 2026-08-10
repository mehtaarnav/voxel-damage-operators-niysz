"""PHASE 1d — dump every table (and every figure caption) from the two papers,
verbatim, from the PMC JATS XML.

No LLM summarisation: cells are read straight out of the XML so the ground-truth
numbers used for the Phase-2 / Phase-4 gates are exactly what is printed in the
papers.
"""

from __future__ import annotations

import os
import re
import sys
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "refs")
OUT = os.path.join(HERE, "out", "phase1")
os.makedirs(OUT, exist_ok=True)

FILES = [
    ("ma8105370  Materials 2015, 8(10), 7129-7147  (TPB paper)",
     "ma8105370_TPB_PMC5455394.efetch.xml"),
    ("ma8095265  Materials 2015, 8(9), 5554-5585  (transport paper)",
     "ma8095265_transport_PMC5512617.efetch.xml"),
]


def txt(el):
    """All descendant text of an element, whitespace-normalised."""
    if el is None:
        return ""
    s = "".join(el.itertext())
    return re.sub(r"\s+", " ", s).strip()


def render_table(tw):
    """Render a <table-wrap> as aligned plain text."""
    lines = []
    label = txt(tw.find(".//label"))
    caption = txt(tw.find(".//caption"))
    lines.append(f"### {label}  {caption}")

    tbl = tw.find(".//table")
    if tbl is None:
        # some tables are supplied only as images
        gr = tw.find(".//graphic")
        lines.append("   [table supplied as an image only"
                     + (f": {gr.attrib}" if gr is not None else "") + "]")
        return lines

    rows = []
    for tr in tbl.iter("tr"):
        cells = []
        for td in tr:
            if td.tag in ("td", "th"):
                span = int(td.attrib.get("colspan", "1"))
                cells.append(txt(td))
                for _ in range(span - 1):
                    cells.append("")
        rows.append(cells)

    if not rows:
        return lines
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]
    widths = [max(len(r[c]) for r in rows) for c in range(ncol)]
    widths = [min(w, 30) for w in widths]
    for r in rows:
        lines.append("   " + " | ".join(r[c][:30].ljust(widths[c])
                                        for c in range(ncol)).rstrip())
    # footnotes
    for fn in tw.iter("table-wrap-foot"):
        t = txt(fn)
        if t:
            lines.append(f"   [foot] {t}")
    return lines


def main():
    all_lines = []
    for title, fname in FILES:
        path = os.path.join(REFS, fname)
        if not os.path.exists(path):
            print("MISSING:", path)
            continue
        tree = ET.parse(path)
        root = tree.getroot()

        block = []
        block.append("=" * 78)
        block.append(title)
        block.append("=" * 78)

        art_title = txt(root.find(".//article-title"))
        block.append(f"Title: {art_title}")
        doi = None
        for aid in root.iter("article-id"):
            if aid.attrib.get("pub-id-type") == "doi":
                doi = aid.text
        block.append(f"DOI:   {doi}")

        block.append("\n--- TABLES ---")
        tws = list(root.iter("table-wrap"))
        if not tws:
            block.append("   (no <table-wrap> elements found)")
        for tw in tws:
            block.extend(render_table(tw))
            block.append("")

        block.append("--- FIGURE CAPTIONS ---")
        for fig in root.iter("fig"):
            lab = txt(fig.find(".//label"))
            cap = txt(fig.find(".//caption"))
            block.append(f"   {lab}: {cap}")
        block.append("")

        block.append("--- SUPPLEMENTARY MATERIAL POINTERS ---")
        for sm in root.iter("supplementary-material"):
            block.append("   " + txt(sm) + "  " + str(sm.attrib))
        block.append("")

        print("\n".join(block))
        all_lines.extend(block)

    dest = os.path.join(OUT, "phase1_paper_tables_verbatim.txt")
    with open(dest, "w", encoding="utf-8") as f:
        f.write("\n".join(all_lines))
    print(f"\n[saved] {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
