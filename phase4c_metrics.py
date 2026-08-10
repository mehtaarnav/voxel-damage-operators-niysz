"""
PHASE 4c — connectivity-margin metrics on every extracted Ni network.

Loads the per-ROI networks saved by phase3_extract_network.py and computes, for
each: weighted algebraic connectivity (raw and normalised), face-to-face
min-cut, effective conductance, and neck-width quantiles.  Then aggregates to
one row per anode as mean +/- standard deviation across ROIs, which is the
error bar the multi-ROI strategy exists to provide.

See cmlib/metrics.py for the size-dependence caveat on raw lambda_2.
"""

from __future__ import annotations

import argparse
import glob
import os
import pickle
import sys
import time

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.metrics import summarise_network  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase4")
NETS = os.path.join(HERE, "out", "networks")
os.makedirs(OUT, exist_ok=True)

ORDER = ["fine_pre", "medium_pre", "coarse_pre",
         "fine_post", "medium_post", "coarse_post"]

METRICS = ["lambda2_raw", "lambda2_norm", "mincut", "g_eff",
           "neck_p10_nm", "neck_p25_nm", "neck_p50_nm",
           "n_nodes", "n_edges", "mean_degree",
           "chamber_equiv_diam_mean_nm"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roi-um", type=float, default=8.0)
    args = ap.parse_args()

    print("=" * 78)
    print("PHASE 4c — connectivity-margin metrics per ROI")
    print("=" * 78)

    files = sorted(glob.glob(os.path.join(NETS, f"*__{args.roi_um}um.pkl")))
    if not files:
        print(f"No networks found in {NETS} for ROI {args.roi_um} um.")
        return 1
    print(f"{len(files)} networks found\n")

    rows = []
    for f in files:
        base = os.path.basename(f)
        sample, roi, _ = base.split("__")
        with open(f, "rb") as fh:
            G = pickle.load(fh)
        t0 = time.time()
        m = summarise_network(G)
        m.update(sample=sample, roi=roi, roi_um=args.roi_um,
                 seconds=round(time.time() - t0, 1))
        rows.append(m)
        print(f"  {sample:12s} {roi:8s} n={m['n_nodes']:5d} e={m['n_edges']:5d}  "
              f"l2={m['lambda2_raw']:10.2f}  l2n={m['lambda2_norm']:.4f}  "
              f"mincut={m['mincut']:10.2f}  g_eff={m['g_eff']:10.2f}  "
              f"p10={m['neck_p10_nm']:6.0f}nm  [{m['seconds']}s]")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, f"phase4c_metrics_per_roi_{args.roi_um}um.csv"),
              index=False)

    # ---- aggregate per anode -------------------------------------------
    print("\n" + "=" * 78)
    print("Aggregated per anode (mean +/- sd across ROIs)")
    print("=" * 78)
    agg_rows = []
    for s in [x for x in ORDER if x in set(df["sample"])]:
        d = df[df["sample"] == s]
        r = {"sample": s, "n_roi": len(d)}
        for m in METRICS:
            if m in d.columns:
                vals = pd.to_numeric(d[m], errors="coerce")
                vals = vals.replace([np.inf, -np.inf], np.nan).dropna()
                r[f"{m}_mean"] = float(vals.mean()) if len(vals) else np.nan
                r[f"{m}_sd"] = float(vals.std(ddof=1)) if len(vals) > 1 else np.nan
                r[f"{m}_n"] = int(len(vals))
        agg_rows.append(r)
    agg = pd.DataFrame(agg_rows)
    agg.to_csv(os.path.join(OUT, f"phase4c_metrics_per_anode_{args.roi_um}um.csv"),
               index=False)

    show = ["sample", "n_roi"]
    for m in ("lambda2_raw", "lambda2_norm", "mincut", "g_eff", "neck_p10_nm",
              "n_nodes"):
        show += [f"{m}_mean", f"{m}_sd"]
    with pd.option_context("display.width", 250, "display.max_columns", 60):
        print(agg[show].to_string(index=False))

    # ---- figure ---------------------------------------------------------
    plot_metrics = ["lambda2_raw", "lambda2_norm", "mincut", "g_eff",
                    "neck_p10_nm", "n_nodes"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    colours = {"fine_pre": "C0", "medium_pre": "C1", "coarse_pre": "C2"}
    for ax, m in zip(axes.ravel(), plot_metrics):
        for i, s in enumerate(agg["sample"]):
            d = df[df["sample"] == s]
            vals = pd.to_numeric(d[m], errors="coerce")
            vals = vals.replace([np.inf, -np.inf], np.nan).dropna()
            c = colours.get(s, f"C{i}")
            ax.scatter(np.full(len(vals), i) + np.random.default_rng(0)
                       .uniform(-0.09, 0.09, len(vals)), vals, s=26,
                       color=c, alpha=0.75, zorder=3)
            if len(vals):
                ax.errorbar(i, vals.mean(),
                            yerr=(vals.std(ddof=1) if len(vals) > 1 else 0),
                            fmt="_", ms=26, color="k", lw=1.6, capsize=6,
                            zorder=4)
        ax.set_xticks(range(len(agg)))
        ax.set_xticklabels(agg["sample"], rotation=15, fontsize=8)
        ax.set_title(m, fontsize=10)
        ax.grid(alpha=0.25, axis="y")
        if m in ("lambda2_raw", "mincut", "g_eff"):
            ax.set_yscale("log")
    fig.suptitle(f"Phase 4 connectivity-margin metrics, {args.roi_um} um ROIs "
                 "(points = individual ROIs, bars = mean +/- sd)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"phase4c_metrics_{args.roi_um}um.png"), dpi=145)
    plt.close(fig)
    print(f"\n[saved] {os.path.join(OUT, f'phase4c_metrics_{args.roi_um}um.png')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
