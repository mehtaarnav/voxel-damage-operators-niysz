"""PHASE 1b — dump the dataset's own metadata files verbatim.

Reads 2_3D_Data_Info.xlsx (every sheet, every non-empty cell) and 1_Read_Me.docx
(every paragraph + table) so that voxel size, stack dimensions and phase-label
convention come from the source, not from assumption.
"""

from __future__ import annotations

import os
import sys
import zipfile
import re

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def dump_xlsx(path):
    print("=" * 78)
    print("2_3D_Data_Info.xlsx")
    print("=" * 78)
    wb = openpyxl.load_workbook(path, data_only=True)
    for ws in wb.worksheets:
        print(f"\n--- sheet: {ws.title!r}   dims={ws.dimensions} ---")
        for row in ws.iter_rows():
            vals = []
            for c in row:
                if c.value is not None and str(c.value).strip() != "":
                    vals.append(f"{c.coordinate}={c.value!r}")
            if vals:
                print("   " + " | ".join(vals))


def dump_docx(path):
    """Minimal docx text extraction (no python-docx dependency)."""
    print("\n" + "=" * 78)
    print("1_Read_Me.docx")
    print("=" * 78)
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8", "replace")
    # split into paragraphs, strip tags, keep text runs
    paras = re.split(r"</w:p>", xml)
    for p in paras:
        texts = re.findall(r"<w:t[^>]*>(.*?)</w:t>", p, flags=re.S)
        line = "".join(texts)
        line = re.sub(r"&amp;", "&", line)
        line = re.sub(r"&lt;", "<", line)
        line = re.sub(r"&gt;", ">", line)
        line = line.strip()
        if line:
            print("   " + line)


def main():
    x = os.path.join(DATA, "2_3D_Data_Info.xlsx")
    d = os.path.join(DATA, "1_Read_Me.docx")
    if os.path.exists(x):
        dump_xlsx(x)
    else:
        print("MISSING:", x)
    if os.path.exists(d):
        dump_docx(d)
    else:
        print("MISSING:", d)
    return 0


if __name__ == "__main__":
    sys.exit(main())
