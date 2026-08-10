"""
PHASE 3c — sensitivity of the connectivity-margin metrics to the SNOW watershed
marker parameter r_max.

r_max (the radius of the maximum filter that merges nearby distance-transform
peaks into one marker) is the dominant over/under-segmentation control, and it
is specified in VOXELS, so the same value means different physical distances on
samples with different voxel sizes (19.53, 24.41 and 29.14 nm here).  If the
ranking of the three anodes depends on r_max, the ranking is an artifact of the
parameter rather than of the microstructure.  This sweep tests exactly that.

One ROI per anode (the first), r_max in {2, 4, 6, 8}.
"""

from __future__ import annotations

import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, load_subvolume  # noqa: E402
from cmlib.metrics import summarise_network  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.pnm import extract_ni_network, largest_component  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase3")
os.makedirs(OUT, exist_ok=True)

PRISTINE = ["fine_pre", "medium_pre", "coarse_pre"]
RMAX = [2, 4, 6, 8]


def main():
    rows = []
    print("=" * 78)
    print("PHASE 3c — SNOW r_max sensitivity")
    print("=" * 78)
    for key in PRISTINE:
        _, folder, grain, state, nx_, ny_, nz_, vx, vy, vz = \
            [s for s in SAMPLES if s[0] == key][0]
        spacing = (vz, vy, vx)
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        r = tile_rois(nz_, ny_, nx_, vz, vy, vx, 8.0, 1)[0]
        vol = load_subvolume(folder, r["z0"], r["z1"], r["y0"], r["y1"],
                             r["x0"], r["x1"])
        ni = vol == mapping["Ni"]
        del vol
        print(f"\n--- {key} (voxel {vx:.2f} nm) ---")
        for rm in RMAX:
            G, diag, extras = extract_ni_network(ni, spacing, axis=2, r_max=rm)
            if G is None or G.number_of_edges() == 0:
                print(f"  r_max={rm}: no network")
                continue
            Gc = largest_component(G)
            Gc.graph["face_lo"] = [n for n in extras["face_lo"] if n in Gc]
            Gc.graph["face_hi"] = [n for n in extras["face_hi"] if n in Gc]
            m = summarise_network(Gc)
            m.update(sample=key, grain=grain, r_max=rm,
                     r_max_nm=rm * float(np.prod(spacing) ** (1 / 3)))
            rows.append(m)
            print(f"  r_max={rm} ({m['r_max_nm']:.0f} nm): nodes={m['n_nodes']:5d} "
                  f"edges={m['n_edges']:5d}  l2={m['lambda2_raw']:9.3f}  "
                  f"mincut={m['mincut']:9.1f}  g_eff={m['g_eff']:8.1f}  "
                  f"neck_p10={m['neck_p10_nm']:6.0f} nm")
        del ni

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase3c_rmax_sensitivity.csv"), index=False)

    print("\n" + "=" * 78)
    print("DOES THE ANODE RANKING SURVIVE THE PARAMETER CHOICE?")
    print("=" * 78)
    for metric, hib in (("neck_p10_nm", True), ("lambda2_raw", True),
                        ("mincut", True), ("g_eff", True), ("n_nodes", True)):
        print(f"\n  {metric}:")
        for rm in RMAX:
            d = df[df.r_max == rm]
            if d.empty:
                continue
            s = d.set_index("grain")[metric].replace([np.inf, -np.inf], np.nan)
            order = list(s.sort_values(ascending=not hib).dropna().index)
            print(f"     r_max={rm}: {order}")

    metrics = ["neck_p10_nm", "lambda2_raw", "mincut", "g_eff"]
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.3))
    for ax, m in zip(axes, metrics):
        for key in PRISTINE:
            d = df[df["sample"] == key].sort_values("r_max")
            ax.plot(d.r_max, d[m], "o-", label=key)
        ax.set_xlabel("SNOW r_max (voxels)")
        ax.set_ylabel(m)
        ax.set_title(m, fontsize=10)
        ax.grid(alpha=0.25)
        if m in ("lambda2_raw", "mincut", "g_eff"):
            ax.set_yscale("log")
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("Phase 3c: sensitivity of the metrics to the watershed marker "
                 "parameter r_max (one ROI per anode)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase3c_rmax_sensitivity.png"), dpi=145)
    plt.close(fig)
    print(f"\n[saved] {os.path.join(OUT, 'phase3c_rmax_sensitivity.csv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
