"""
PHASE 4b — TPB length density on the FULL stacks (all six samples).

Runs at native resolution over every voxel of every stack (0.48-1.11 gigavoxel)
by streaming two z-slices at a time, so no sub-sampling and no REV compromise
is involved here.  This is the measurement that the Phase-4 gate compares
against the published Figure 7 values.

Convention: voxel-edge counting; see cmlib/tpb.py for the full statement and
phase4a_validate_tpb.py for its validation (exact on axis-aligned test cases;
staircase over-estimate measured at 1.713 for the worst-case (1,1,1) orientation
against a theoretical bound of sqrt(3) = 1.732).
"""

from __future__ import annotations

import os
import sys
import time

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE, ground_truth_frame  # noqa: E402
from cmlib.io import iter_slices, label_histogram  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.tpb import tpb_density_streaming  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase4")
os.makedirs(OUT, exist_ok=True)


def main():
    print("=" * 78)
    print("PHASE 4b — TPB density, full stacks, voxel-edge convention")
    print("=" * 78)

    gt = ground_truth_frame().set_index("sample")
    rows = []
    for key, folder, grain, state, nx_, ny_, nz_, vx, vy, vz in SAMPLES:
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        lab = {"Ni": mapping["Ni"], "YSZ": mapping["YSZ"],
               "pore": mapping["pore"]}
        spacing = (vz, vy, vx)
        t0 = time.time()
        r = tpb_density_streaming(iter_slices(folder), lab, spacing)
        dt = time.time() - t0

        pub = float(gt.loc[key, "TPB_total_um-2__P-F7_digitized"])
        mine = r["tpb_density_um-2"]
        ratio = mine / pub
        print(f"\n--- {key} ({grain}, {state}) ---")
        print(f"  voxels {r['nz']}x{r['ny']}x{r['nx']}   "
              f"volume {r['volume_um3']:.1f} um^3   [{dt:.0f} s]")
        print(f"  TPB edges (z,y,x) = ({r['tpb_edges_z']:,}, "
              f"{r['tpb_edges_y']:,}, {r['tpb_edges_x']:,})")
        print(f"  TPB length        = {r['tpb_length_um']:,.0f} um")
        print(f"  TPB density       = {mine:.4f} um^-2   "
              f"(published total {pub:.2f})   ratio = {ratio:.3f}")

        row = dict(sample=key, grain=grain, state=state, **r,
                   tpb_published_total=pub,
                   tpb_published_active=float(
                       gt.loc[key, "TPB_active_um-2__P-F7_digitized"]),
                   ratio_mine_over_published=ratio,
                   seconds=round(dt, 1))
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase4b_tpb_full_stacks.csv"), index=False)

    print("\n" + "=" * 78)
    print("PHASE 4 TPB GATE")
    print("=" * 78)
    pri = df[df.state == "pristine"]
    print(df[["sample", "grain", "state", "tpb_density_um-2",
              "tpb_published_total", "ratio_mine_over_published"]]
          .to_string(index=False))

    within2 = ((df["ratio_mine_over_published"] > 0.5) &
               (df["ratio_mine_over_published"] < 2.0)).all()
    # ordering must be fine > medium > coarse for the pristine samples
    order = pri.sort_values("tpb_density_um-2", ascending=False)["grain"].tolist()
    order_ok = order == ["fine", "medium", "coarse"]
    rmin = df["ratio_mine_over_published"].min()
    rmax = df["ratio_mine_over_published"].max()
    print(f"\n  ratio range over all six samples: {rmin:.3f} - {rmax:.3f}")
    print(f"  all within a factor of 2          : {within2}")
    print(f"  pristine ordering by TPB density  : {order}  "
          f"-> {'PASS' if order_ok else 'FAIL'}")
    print(f"  ratio spread (max/min)            : {rmax/rmin:.3f}   "
          "(a CONSTANT ratio would indicate a pure convention offset)")
    print(f"\n  RESULT: {'PASS' if (within2 and order_ok) else 'FAIL'}")

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    w = 0.35
    xs = np.arange(len(df))
    axes[0].bar(xs - w/2, df["tpb_density_um-2"], w, label="this work (voxel-edge)")
    axes[0].bar(xs + w/2, df["tpb_published_total"], w,
                label="published TPB$_{total}$ (Fig. 7)")
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels([f"{g}\n{s}" for g, s in zip(df.grain, df.state)],
                            fontsize=8)
    axes[0].set_ylabel("TPB density (um$^{-2}$)")
    axes[0].set_title("TPB density: this work vs published")
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].grid(alpha=0.25, axis="y")

    axes[1].bar(xs, df["ratio_mine_over_published"], color="C2")
    axes[1].axhline(1.0, color="k", lw=1)
    axes[1].axhline(1.713, color="crimson", ls="--", lw=1.3,
                    label="measured worst-case staircase factor 1.713")
    axes[1].axhline(df["ratio_mine_over_published"].mean(), color="C0", ls=":",
                    lw=1.5,
                    label=f"mean ratio = {df['ratio_mine_over_published'].mean():.2f}")
    axes[1].set_xticks(xs)
    axes[1].set_xticklabels([f"{g}\n{s}" for g, s in zip(df.grain, df.state)],
                            fontsize=8)
    axes[1].set_ylabel("this work / published")
    axes[1].set_title("Ratio — a constant value means a pure convention offset")
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].grid(alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase4b_tpb.png"), dpi=145)
    plt.close(fig)
    print(f"\n[saved] {os.path.join(OUT, 'phase4b_tpb.png')}")
    return 0 if (within2 and order_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
