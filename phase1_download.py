"""
PHASE 1a — Download the Zenodo dataset (metadata + segmented stacks only).

Zenodo record 4056538:
  "FIB-tomography data of Ni-YSZ anodes for Solid Oxide Fuel Cells (SOFC):
   Comparison of pristine and degraded materials (before/after redox cycling)"
  Holzer, Pecho, Hocker, Iwanschitz, Mai.

DELIBERATE SCOPE CHOICE: only the *_Segmented.zip archives (~140 MB total) and
the two metadata files are fetched.  The *_RawData.zip archives total ~11 GB and
the *_Visualization.zip archives ~660 MB; neither is needed, because the whole
study operates on the already phase-segmented volumes.  If the segmentation
turns out to be ambiguous we can revisit and pull raw data for one sample.

Usage:  python phase1_download.py [--meta-only]
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import requests

RECORD = "4056538"
BASE = f"https://zenodo.org/records/{RECORD}/files/"

META_FILES = ["1_Read_Me.docx", "2_3D_Data_Info.xlsx"]
SEG_FILES = [
    "3_Rx36_Segmented.zip",
    "4_Rx37_Segmented.zip",
    "5_Rx38_Segmented.zip",
    "6_Rx41-1_Segmented.zip",
    "7_Rx41-2_Segmented.zip",
    "8_Rx41-3_Segmented.zip",
]

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def download(fname: str, dest_dir: str) -> str:
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, fname)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  [skip] {fname}  already present ({os.path.getsize(dest)/1e6:.1f} MB)")
        return dest
    url = BASE + fname + "?download=1"
    print(f"  [get ] {fname}")
    t0 = time.time()
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        got = 0
        tmp = dest + ".part"
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                got += len(chunk)
        os.replace(tmp, dest)
    dt = time.time() - t0
    print(f"         {got/1e6:8.1f} MB in {dt:5.1f} s"
          + (f"  (expected {total/1e6:.1f} MB)" if total else ""))
    return dest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta-only", action="store_true")
    args = ap.parse_args()

    print("=" * 78)
    print(f"PHASE 1a — downloading from Zenodo record {RECORD}")
    print("=" * 78)

    print("\nMetadata:")
    for f in META_FILES:
        download(f, DATA)

    if args.meta_only:
        print("\n--meta-only: stopping before the segmented stacks.")
        return 0

    print("\nSegmented stacks (~140 MB total):")
    for f in SEG_FILES:
        download(f, DATA)

    print("\nDone. Files in:", DATA)
    return 0


if __name__ == "__main__":
    sys.exit(main())
