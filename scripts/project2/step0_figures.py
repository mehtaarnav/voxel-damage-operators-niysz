"""Regenerate the Phase 5 figure (now complete, both phases) from the CSVs.

Separated from step0_percolation.py so the figure can be rebuilt without
re-reading 2.1 GB of TIFF stacks. The previously committed
out/phase5/phase5_percolation.png was drawn from a single sample and is
replaced here.
"""
from __future__ import annotations

import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT5 = os.path.join(ROOT, "out", "phase5")
OUT2 = os.path.join(ROOT, "out", "project2")

ORDER = ["fine", "medium", "coarse"]


def main():
    ni = pd.read_csv(os.path.join(OUT5, "phase5_percolation.csv"))
    ysz = pd.read_csv(os.path.join(OUT2, "step0_ysz_percolation.csv"))
    nir = pd.read_csv(os.path.join(OUT5, "phase5_retention.csv"))
    yr = pd.read_csv(os.path.join(OUT2, "step0_ysz_retention.csv"))
    nir = nir.set_index("grain").loc[ORDER].reset_index()
    yr = yr.set_index("grain").loc[ORDER].reset_index()

    # --- Phase 5 figure: Ni, as before, but complete -----------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    xs, w = np.arange(len(ni)), 0.27
    axes[0].bar(xs - w, ni.P_span, w, label="P_span (both faces)")
    axes[0].bar(xs, ni.P_reach, w, label="P_reach (either face)")
    axes[0].bar(xs + w, ni.P_published, w, label="published P (MIP-PSD)")
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels([f"{g}\n{s}" for g, s in zip(ni.grain, ni.state)],
                            fontsize=8)
    axes[0].set_ylabel("percolating fraction of Ni")
    axes[0].set_title("Ni percolation, full stacks")
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].grid(alpha=0.25, axis="y")

    xs2 = np.arange(len(nir))
    axes[1].bar(xs2 - 0.2, nir.P_reach_retained, 0.4, label="this work (P_reach)")
    axes[1].bar(xs2 + 0.2, nir.P_pub_retained, 0.4, label="published P")
    axes[1].set_xticks(xs2)
    axes[1].set_xticklabels(nir.grain)
    axes[1].axhline(1.0, color="k", lw=1)
    axes[1].set_ylabel("retained fraction (degraded / pristine)")
    axes[1].set_title("Retention of Ni percolation after redox cycling")
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].grid(alpha=0.25, axis="y")
    fig.tight_layout()
    p = os.path.join(OUT5, "phase5_percolation.png")
    fig.savefig(p, dpi=145)
    plt.close(fig)
    print("[saved]", p)

    # --- Step 0 figure: the divergence ------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    xs = np.arange(len(ysz))
    axes[0].bar(xs - 0.2, ysz.P_span, 0.4, color="#4C72B0", label="YSZ P_span")
    axes[0].bar(xs + 0.2, ysz.P_largest, 0.4, color="#A7BEDE",
                label="YSZ P_largest")
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels([f"{g}\n{s}" for g, s in zip(ysz.grain, ysz.state)],
                            fontsize=8)
    axes[0].set_ylabel("fraction of YSZ phase")
    axes[0].set_title("YSZ percolation (Step 0, first measurement)")
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].grid(alpha=0.25, axis="y")
    axes[0].annotate("no spanning\ncluster", xy=(5, 0.02), xytext=(4.3, 0.35),
                     fontsize=8, ha="center",
                     arrowprops=dict(arrowstyle="->", lw=1))

    xs2 = np.arange(len(ORDER))
    axes[1].bar(xs2 - 0.2, nir.P_span_retained, 0.4, color="#C44E52",
                label="Ni P_span retained")
    axes[1].bar(xs2 + 0.2, yr.P_span_retained, 0.4, color="#4C72B0",
                label="YSZ P_span retained")
    axes[1].set_xticks(xs2)
    axes[1].set_xticklabels(ORDER)
    axes[1].axhline(1.0, color="k", lw=1)
    axes[1].set_ylabel("retained (degraded / pristine)")
    axes[1].set_title("The divergence: Ni worst in fine, YSZ worst in coarse")
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].grid(alpha=0.25, axis="y")

    post = ysz[ysz.state == "degraded"].set_index("grain").loc[ORDER]
    axes[2].scatter(post.volume_fraction, post.P_span, s=90, color="#4C72B0",
                    zorder=3)
    for g, r in post.iterrows():
        axes[2].annotate(g, (r.volume_fraction, r.P_span),
                         textcoords="offset points", xytext=(8, 4), fontsize=9)
    axes[2].axvline(0.3116, color="k", ls="--", lw=1)
    axes[2].annotate("random-site $p_c$ = 0.3116", xy=(0.3116, 0.5),
                     xytext=(4, 0), textcoords="offset points", rotation=90,
                     fontsize=7, va="center")
    axes[2].set_xlabel(r"degraded $\Phi_{YSZ}$")
    axes[2].set_ylabel("YSZ P_span (degraded)")
    axes[2].set_title("YSZ collapse is not a volume-fraction effect")
    axes[2].grid(alpha=0.25)
    fig.tight_layout()
    p = os.path.join(OUT2, "step0_percolation.png")
    fig.savefig(p, dpi=145)
    plt.close(fig)
    print("[saved]", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
