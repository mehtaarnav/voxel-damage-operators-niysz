"""26-connectivity sensitivity check.

The paper evaluates connectivity at 6-connectivity throughout, which is the
conservative convention: face contacts only, no face-diagonal or corner
contacts. Section 4.2 states the concern that this is not neutral across the
series -- the finest anode has the most surface per unit volume and therefore
the most opportunity for diagonal contacts, so its connectivity should be the
most understated -- and records that the check was not run.

This runs it. For each anode, the pristine spanning fraction is computed under
both conventions on the full ROI, and the difference is compared across the
series. If the concern is real, the fine anode gains most.

Reported per ROI so the spread is visible rather than averaged away.
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import tifffile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE      # noqa: E402
from cmlib.io import label_histogram, slice_paths              # noqa: E402
from cmlib.percolation import percolating_mask                 # noqa: E402
from cmlib.phases import assign_labels                         # noqa: E402
from cmlib.roi import tile_rois                                # noqa: E402

OUT = os.path.join(ROOT, "out", "project2")
AXIS = 2
SIDE = {"fine": 8.0, "medium": 8.0, "coarse": 12.0}
NROI = 3


def load(folder, mapping, r):
    ps = slice_paths(folder)
    sh = (r["z1"] - r["z0"], r["y1"] - r["y0"], r["x1"] - r["x0"])
    ni = np.empty(sh, bool)
    for i, z in enumerate(range(r["z0"], r["z1"])):
        ni[i] = tifffile.imread(ps[z])[r["y0"]:r["y1"],
                                       r["x0"]:r["x1"]] == mapping["Ni"]
    return ni


def main():
    sample = {s[2]: s for s in SAMPLES if s[3] == "pristine"}
    rows = []
    for grain in ("fine", "medium", "coarse"):
        k = sample[grain]
        mapping = assign_labels(label_histogram(k[1])["counts"],
                                ZENODO_LABEL_NOTE[k[0]])
        for ri, r in enumerate(tile_rois(k[6], k[5], k[4], k[9], k[8], k[7],
                                         SIDE[grain], max_rois=NROI)):
            t0 = time.time()
            ni = load(k[1], mapping, r)
            n0 = int(ni.sum())
            p6 = int(percolating_mask(ni, axis=AXIS, connectivity=6).sum()) / n0
            p26 = int(percolating_mask(ni, axis=AXIS, connectivity=26).sum()) / n0
            rows.append(dict(anode=grain, roi=ri, p_span_6=p6, p_span_26=p26,
                             gain=p26 - p6,
                             disconnected_6=1 - p6, disconnected_26=1 - p26))
            print(f"  {grain:7s} roi{ri}  6-conn {p6:.4f}  26-conn {p26:.4f}"
                  f"  gain {p26-p6:+.4f}   [{time.time()-t0:.0f}s]",
                  flush=True)
            del ni
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/o7_connectivity_check.csv", index=False)

    g = df.groupby("anode")[["p_span_6", "p_span_26", "gain"]].agg(
        ["mean", "std"]).round(4)
    print("\nclass means")
    print(g.to_string())

    m = df.groupby("anode").gain.mean().sort_values(ascending=False)
    print(f"\ngain from 26-connectivity, largest first: "
          f"{', '.join(f'{a} {v:+.4f}' for a, v in m.items())}")
    print(f"the stated concern predicts fine gains most -> "
          f"{'CONFIRMED' if m.index[0] == 'fine' else 'NOT CONFIRMED'}")
    dis6 = df.groupby("anode").disconnected_6.mean()
    dis26 = df.groupby("anode").disconnected_26.mean()
    print(f"\ndisconnected Ni fraction, 6-conn:  "
          f"{', '.join(f'{a} {v*100:.1f}%' for a, v in dis6.items())}")
    print(f"disconnected Ni fraction, 26-conn: "
          f"{', '.join(f'{a} {v*100:.1f}%' for a, v in dis26.items())}")


if __name__ == "__main__":
    main()
