"""Threshold sensitivity of the erosion transition.

Reads out/project2/o7_threshold_transitions.csv (27 bisection-equivalent
transitions per threshold: 3 ROIs x 3 seeds x 3 anodes, read off a single
damage realisation per ROI and seed so that all thresholds see identical
damage).

The point of the panel is not the ordering, which the spread swallows, but
that the fine anode never fails first at any threshold, while the measurement
says it is the first to lose percolation.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotstyle import ANODES, COLOR, MARKER, COL_W, tidy, use

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P2 = os.path.join(ROOT, "out", "project2")
FIG = os.path.join(ROOT, "out", "writeup", "figs")
use()


def main():
    t = pd.read_csv(f"{P2}/o7_threshold_transitions.csv")
    g = t.groupby(["threshold", "anode"]).transition.agg(["mean", "std"])

    fig, ax = plt.subplots(figsize=(COL_W, 2.5))
    for a in ANODES:
        s = g.xs(a, level="anode").sort_index()
        ax.errorbar(s.index, s["mean"], yerr=s["std"], color=COLOR[a],
                    marker=MARKER[a], ms=3.0, lw=1.1, capsize=1.8,
                    elinewidth=0.6, label=a, zorder=3)
    ax.set_xlabel(r"$\Rni$ threshold".replace(r"\Rni",
                                              r"$R_{\mathrm{Ni}}$"))
    ax.set_xlabel(r"$R_{\mathrm{Ni}}$ threshold")
    ax.set_ylabel("damage rounds to transition")
    ax.invert_xaxis()
    ax.legend(loc="upper left")
    ax.set_title("the fine anode never fails first")
    tidy(ax)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"fig6_threshold_sweep.{ext}"))
    plt.close(fig)
    print("  wrote fig6_threshold_sweep")

    # the statement the figure is asserting, checked rather than assumed
    order = g["mean"].unstack()
    first = order.idxmin(axis=1)
    print("  first to fail, by threshold:")
    for thr, who in first.items():
        print(f"    {thr:.2f}  {who}")
    print(f"  fine is first to fail at "
          f"{int((first == 'fine').sum())}/{len(first)} thresholds")


if __name__ == "__main__":
    main()
