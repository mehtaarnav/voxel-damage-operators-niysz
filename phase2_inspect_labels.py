"""
PHASE 2a — inspect the segmented stacks BEFORE interpreting any label.

Several filenames say "4Phases" while the metadata says three phases, and one
stack (5_Rx38) is documented with Ni and YSZ swapped.  So: read the actual label
values and their exact voxel counts first, assume nothing.

Streams slice-by-slice; never holds a full volume in RAM.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, stack_shape  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase2")
os.makedirs(OUT, exist_ok=True)


def main():
    print("=" * 78)
    print("PHASE 2a — raw label inventory of the six segmented stacks")
    print("=" * 78)

    rows = []
    for key, folder, grain, state, nx, ny, nz, vx, vy, vz in SAMPLES:
        print(f"\n--- {key}  ({folder}) ---")
        shp = stack_shape(folder)
        print(f"  stack shape (nz, ny, nx) from files : {shp}")
        print(f"  metadata (nz, ny, nx)               : ({nz}, {ny}, {nx})")
        match = (shp == (nz, ny, nx))
        print(f"  shape matches metadata              : {match}")
        if not match:
            print("  *** SHAPE MISMATCH — investigate before trusting anything ***")

        h = label_histogram(folder)
        c = h["counts"]
        present = np.nonzero(c)[0]
        tot = h["total"]
        print(f"  dtype(s)                            : {h['dtypes']}")
        print(f"  in-plane shape(s)                   : {h['shapes']}")
        print(f"  total voxels                        : {tot:,}")
        print(f"  metadata label note                 : {ZENODO_LABEL_NOTE[key]}")
        print(f"  DISTINCT LABEL VALUES PRESENT       : {list(present)}")
        for v in present:
            print(f"      value {v:3d} : {c[v]:14,d}  ({100.0*c[v]/tot:6.3f} %)")

        row = {"sample": key, "folder": folder, "grain": grain, "state": state,
               "shape_files": str(shp), "shape_meta": f"({nz}, {ny}, {nx})",
               "shape_match": match, "dtype": ",".join(h["dtypes"]),
               "total_voxels": tot, "n_labels": len(present),
               "labels": ",".join(str(v) for v in present),
               "zenodo_note": ZENODO_LABEL_NOTE[key]}
        for i, v in enumerate(present):
            row[f"label{i}_value"] = int(v)
            row[f"label{i}_count"] = int(c[v])
            row[f"label{i}_frac"] = float(c[v] / tot)
        rows.append(row)

    df = pd.DataFrame(rows)
    dest = os.path.join(OUT, "phase2_label_inventory.csv")
    df.to_csv(dest, index=False)
    print(f"\n[saved] {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
