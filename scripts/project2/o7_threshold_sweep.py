"""Is the erosion sign reversal robust to the choice of R_Ni threshold?

The pre-registered bisection used two thresholds, 0.50 and 0.10, and found that
the fine anode is the LAST to fail under uniform surface erosion -- the reverse
of the measured ordering, in which fine retains percolation worst. A referee
can reasonably ask whether that reversal is an artifact of where the threshold
was placed, particularly since the measured fine retention is 0.680, above both
pre-registered thresholds.

This sweeps the threshold instead of assuming it.

METHOD. cmlib.damage2.apply_o6 erodes round by round from a single RNG stream
and prunes to the largest component once at the end. Stepping one round at a
time with one stream therefore reproduces apply_o6(n) exactly for every n in a
single pass, at a twentieth of the cost of re-running it per threshold. The
R_Ni curve is recorded for n = 1..20, and every threshold is then read off the
same curve, so all thresholds see identical damage realisations.

R_Ni is the frozen metric: spanning-cluster voxels divided by PRISTINE voxels
(out/project2/PREREG_RNI_METRIC.md), which is invariant to the pruning step.

Design matches the pre-registered bisection: 3 ROIs per anode, damage seeds
300/301/302, p_erode = 0.35, n in [1, 20].
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import tifffile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.damage2 import STRUCT6                              # noqa: E402
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE      # noqa: E402
from cmlib.io import label_histogram, slice_paths              # noqa: E402
from cmlib.percolation import percolating_mask                 # noqa: E402
from cmlib.phases import assign_labels                         # noqa: E402
from cmlib.roi import tile_rois                                # noqa: E402

OUT = os.path.join(ROOT, "out", "project2")
AXIS, CONN = 2, 6
SIDE = {"fine": 8.0, "medium": 8.0, "coarse": 12.0}
SEEDS = (300, 301, 302)
NROI = 3
P_ERODE = 0.35
NMAX = 20
THRESHOLDS = (0.80, 0.70, 0.68, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)


def load(folder, mapping, r):
    ps = slice_paths(folder)
    sh = (r["z1"] - r["z0"], r["y1"] - r["y0"], r["x1"] - r["x0"])
    ni = np.empty(sh, bool)
    for i, z in enumerate(range(r["z0"], r["z1"])):
        ni[i] = tifffile.imread(ps[z])[r["y0"]:r["y1"],
                                       r["x0"]:r["x1"]] == mapping["Ni"]
    return ni


def curve(ni0, seed):
    """R_Ni after each of NMAX erosion rounds, semantics of apply_o6."""
    n0 = int(ni0.sum())
    rng = np.random.default_rng(seed)
    cur = ni0.copy()
    out = []
    for n in range(1, NMAX + 1):
        eroded = ndi.binary_erosion(cur, structure=STRUCT6)
        boundary = cur & ~eroded
        cur &= ~(boundary & (rng.random(cur.shape) < P_ERODE))
        lab, k = ndi.label(cur, structure=STRUCT6)
        if k == 0:
            out.append(0.0)
            break
        c = np.bincount(lab.ravel()); c[0] = 0
        final = lab == int(np.argmax(c))
        out.append(int(percolating_mask(final, axis=AXIS,
                                        connectivity=CONN).sum()) / n0)
        del lab, final
    while len(out) < NMAX:
        out.append(0.0)
    return out


def main():
    sample = {s[2]: s for s in SAMPLES if s[3] == "pristine"}
    rows = []
    for grain in ("fine", "medium", "coarse"):
        k = sample[grain]
        mapping = assign_labels(label_histogram(k[1])["counts"],
                                ZENODO_LABEL_NOTE[k[0]])
        rois = tile_rois(k[6], k[5], k[4], k[9], k[8], k[7],
                         SIDE[grain], max_rois=NROI)
        for ri, r in enumerate(rois):
            t0 = time.time()
            ni0 = load(k[1], mapping, r)
            for seed in SEEDS:
                cv = curve(ni0, seed)
                for n, v in enumerate(cv, 1):
                    rows.append(dict(anode=grain, roi=ri, seed=seed,
                                     n_rounds=n, r_ni=v))
            print(f"  {grain:7s} roi{ri}  {time.time()-t0:.0f}s", flush=True)
            del ni0
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/o7_threshold_curves.csv", index=False)
    print("curves written")

    # transitions: first n at which R_Ni falls below each threshold
    tr = []
    for (a, ri, sd), g in df.groupby(["anode", "roi", "seed"]):
        g = g.sort_values("n_rounds")
        for T in THRESHOLDS:
            below = g[g.r_ni < T]
            if below.empty:
                n, flag = NMAX + 0.5, "ceiling"
            else:
                n, flag = float(below.n_rounds.iloc[0]) - 0.5, "ok"
            tr.append(dict(anode=a, roi=ri, seed=sd, threshold=T,
                           transition=n, flag=flag))
    t = pd.DataFrame(tr)
    t.to_csv(f"{OUT}/o7_threshold_transitions.csv", index=False)

    print(f"\n{'thresh':>7} {'fine':>7} {'medium':>7} {'coarse':>7} "
          f"{'ordering (latest to fail first)':>34} {'fine last?':>11}")
    print("-" * 78)
    summary = []
    for T in THRESHOLDS:
        s = t[t.threshold == T].groupby("anode").transition.mean()
        order = s.sort_values(ascending=False)
        fine_last = order.index[0] == "fine"
        summary.append(dict(threshold=T, fine=s["fine"], medium=s["medium"],
                            coarse=s["coarse"],
                            ordering=" > ".join(order.index),
                            fine_most_resistant=bool(fine_last)))
        print(f"{T:7.2f} {s['fine']:7.2f} {s['medium']:7.2f} {s['coarse']:7.2f} "
              f"{' > '.join(order.index):>34} {str(fine_last):>11}")
    pd.DataFrame(summary).to_csv(f"{OUT}/o7_threshold_summary.csv", index=False)
    n_ok = sum(r["fine_most_resistant"] for r in summary)
    print(f"\nfine is the most erosion-resistant at {n_ok}/{len(summary)} "
          f"thresholds")
    print("measured ordering is the opposite: fine retains percolation WORST")


if __name__ == "__main__":
    main()
