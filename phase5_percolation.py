"""
PHASE 5 — degraded-state outcomes: does Ni still percolate, and how much of it?

Run on the FULL stacks (0.48-1.11 gigavoxel) at native resolution, one sample at
a time, using the Phase-0-validated connected-component + face-spanning code.
No sub-volume, no REV compromise: this is the outcome variable the hypothesis
must predict, so it is measured on everything.

THREE DIFFERENT "PERCOLATING FRACTIONS" ARE REPORTED, because the papers'
P is not defined the same way as a spanning fraction:

  P_span    fraction of Ni voxels in cluster(s) touching BOTH opposite faces of
            the transport axis.  Strictest; this is our Phase-0-validated
            definition.
  P_reach   fraction of Ni voxels in cluster(s) touching AT LEAST ONE of those
            two faces.  This is the closest analogue to the papers' MIP-PSD
            derived P, which measures the fraction of a phase reachable from a
            boundary by simulated intrusion.  P_reach >= P_span always.
  P_largest largest single cluster / total Ni.

Comparing our numbers against the published P therefore uses P_reach as the
like-for-like column, with P_span reported alongside as the stricter bound.

Transport axis = x = array axis 2 (through-thickness; see cmlib/graph.py note 2).
Connectivity = 6 (face-sharing), as validated in Phase 0.
"""

from __future__ import annotations

import argparse
import gc
import os
import sys
import time

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tifffile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE, ground_truth_frame  # noqa: E402
from cmlib.io import label_histogram, slice_paths  # noqa: E402
from cmlib.percolation import percolation_summary  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase5")
os.makedirs(OUT, exist_ok=True)

AXIS = 2          # x, through-thickness
CONN = 6


def load_ni_mask(folder, ni_label):
    """Full-stack boolean Ni mask, built slice by slice (1 byte/voxel)."""
    ps = slice_paths(folder)
    a0 = tifffile.imread(ps[0])
    out = np.empty((len(ps), a0.shape[0], a0.shape[1]), dtype=bool)
    out[0] = a0 == ni_label
    for i, p in enumerate(ps[1:], start=1):
        out[i] = tifffile.imread(p) == ni_label
    return out


def analyse(mask, axis=AXIS, conn=CONN):
    """Thin wrapper over cmlib.percolation.percolation_summary.

    The P_span/P_reach/P_largest/percolation logic previously lived here
    (moved to the library, unchanged in definition, 2026-08-10 -- see
    out/next/EXECUTION_SPEC.md Phase 0). This function now only renames keys
    to the column names phase5_percolation.csv has always used, so a re-run
    reproduces the exact same CSV schema.
    """
    s = percolation_summary(mask, axis=axis, connectivity=conn,
                            check_other_axes=True)
    res = {
        "n_ni_voxels": s["n_phase_voxels"],
        "n_clusters": s["n_clusters"],
        "percolates_x": s["percolates"],
        "n_spanning_clusters": s["n_spanning_clusters"],
        "P_span": s["P_span"],
        "P_reach": s["P_reach"],
        "P_largest": s["P_largest"],
        "percolates_z": s["percolates_axis0"],
        "percolates_y": s["percolates_axis1"],
    }
    gc.collect()
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="*",
                    default=[s[0] for s in SAMPLES])
    args = ap.parse_args()

    print("=" * 78)
    print("PHASE 5 — full-stack Ni percolation (6-connectivity, x axis)")
    print("=" * 78)

    gt = ground_truth_frame().set_index("sample")
    rows = []
    for key in args.samples:
        _, folder, grain, state, nx_, ny_, nz_, vx, vy, vz = \
            [s for s in SAMPLES if s[0] == key][0]
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        print(f"\n--- {key} ({grain}, {state}) ---")
        t0 = time.time()
        mask = load_ni_mask(folder, mapping["Ni"])
        print(f"  mask {mask.shape} = {mask.size/1e6:.0f} Mvoxel, "
              f"{mask.nbytes/1e9:.2f} GB   [{time.time()-t0:.0f} s]")
        t0 = time.time()
        r = analyse(mask)
        dt = time.time() - t0
        del mask
        gc.collect()

        pub_P = float(gt.loc[key, "Ni_P__T-S4"])
        pub_phi = float(gt.loc[key, "Ni_Phi__T-S4"])
        print(f"  clusters {r['n_clusters']:,}   spans x: {r['percolates_x']}"
              f"  (y: {r['percolates_y']}, z: {r['percolates_z']})")
        print(f"  P_span   = {r['P_span']:.4f}")
        print(f"  P_reach  = {r['P_reach']:.4f}   <- like-for-like vs published")
        print(f"  P_largest= {r['P_largest']:.4f}")
        print(f"  published Ni P = {pub_P:.3f}   "
              f"(P_reach/published = {r['P_reach']/pub_P:.3f})   [{dt:.0f} s]")

        rows.append(dict(sample=key, grain=grain, state=state,
                         phi_ni_published=pub_phi, P_published=pub_P, **r,
                         P_reach_over_published=r["P_reach"] / pub_P,
                         seconds=round(dt, 1)))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase5_percolation.csv"), index=False)

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    cols = ["sample", "grain", "state", "n_clusters", "percolates_x",
            "P_span", "P_reach", "P_largest", "P_published",
            "P_reach_over_published"]
    print(df[cols].to_string(index=False))

    # retention: degraded / pristine, per grain
    print("\nRetention (degraded / pristine), same grain:")
    ret = []
    for g in ("fine", "medium", "coarse"):
        a = df[(df.grain == g) & (df.state == "pristine")]
        b = df[(df.grain == g) & (df.state == "degraded")]
        if len(a) and len(b):
            a, b = a.iloc[0], b.iloc[0]
            ret.append(dict(
                grain=g,
                P_span_pre=a.P_span, P_span_post=b.P_span,
                P_span_retained=b.P_span / a.P_span if a.P_span else np.nan,
                P_reach_pre=a.P_reach, P_reach_post=b.P_reach,
                P_reach_retained=b.P_reach / a.P_reach if a.P_reach else np.nan,
                P_pub_pre=a.P_published, P_pub_post=b.P_published,
                P_pub_retained=b.P_published / a.P_published))
    rdf = pd.DataFrame(ret)
    rdf.to_csv(os.path.join(OUT, "phase5_retention.csv"), index=False)
    print(rdf.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    xs = np.arange(len(df))
    w = 0.27
    axes[0].bar(xs - w, df.P_span, w, label="P_span (both faces)")
    axes[0].bar(xs, df.P_reach, w, label="P_reach (either face)")
    axes[0].bar(xs + w, df.P_published, w, label="published P (MIP-PSD)")
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels([f"{g}\n{s}" for g, s in zip(df.grain, df.state)],
                            fontsize=8)
    axes[0].set_ylabel("percolating fraction of Ni")
    axes[0].set_title("Ni percolation, full stacks")
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].grid(alpha=0.25, axis="y")

    if len(rdf):
        xs2 = np.arange(len(rdf))
        axes[1].bar(xs2 - 0.2, rdf.P_reach_retained, 0.4,
                    label="this work (P_reach)")
        axes[1].bar(xs2 + 0.2, rdf.P_pub_retained, 0.4, label="published P")
        axes[1].set_xticks(xs2)
        axes[1].set_xticklabels(rdf.grain)
        axes[1].axhline(1.0, color="k", lw=1)
        axes[1].set_ylabel("retained fraction (degraded / pristine)")
        axes[1].set_title("Retention of Ni percolation after redox cycling")
        axes[1].legend(fontsize=8, frameon=False)
        axes[1].grid(alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase5_percolation.png"), dpi=145)
    plt.close(fig)
    print(f"\n[saved] {os.path.join(OUT, 'phase5_percolation.png')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
