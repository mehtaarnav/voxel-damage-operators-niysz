"""Manuscript figures 7 and 8, in the print style of the other six.

Every value is read from a committed CSV written by o7_export_tables.py.

Figure 7 (double column, two panels) carries the central result: the gated
quantity falls while the quantity that matters climbs, and the reason, which is
that most accepted moves cost nothing under the gate.

Figure 8 (single column) isolates where the triple-line inflation comes from by
re-applying the two groups of moves separately.

Neither panel uses a second y-axis. Where two quantities of different scale must
be compared they are divided by their pristine values and share one axis.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import NullFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotstyle import ACCENT, COLOR, COL_W, FULL_W, panel, tidy, use

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P2 = os.path.join(ROOT, "out", "project2")
FIG = os.path.join(ROOT, "out", "writeup", "figs")
use()

NI_C = COLOR["fine"]        # the gated quantity, in the anode ramp
TPB_C = ACCENT              # reserved in this paper for violations and targets


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig)
    print("  wrote", name)


def figure7():
    tr = pd.read_csv(f"{P2}/o7_trajectory.csv")
    hs = pd.read_csv(f"{P2}/o7_da_histogram.csv")

    fig, ax = plt.subplots(1, 2, figsize=(FULL_W, 2.5))

    a = ax[0]
    a.axhline(1.0, color="0.55", lw=0.6, zorder=1)
    a.plot(tr.accepted / 1e3, tr.tpb_ratio, color=TPB_C, marker="o",
           ms=2.6, zorder=3, label="TPB density")
    a.plot(tr.accepted / 1e3, tr.s_ratio, color=NI_C, marker="s",
           ms=2.6, zorder=3, label=r"$S_{\mathrm{spec}}$ (gated)")
    a.set_yscale("log")
    a.set_yticks([1, 2, 3, 5])
    a.set_yticklabels(["1", "2", "3", "5"])
    a.yaxis.set_minor_formatter(NullFormatter())
    a.set_xlabel(r"accepted moves ($\times 10^{3}$)")
    a.set_ylabel("value / pristine value")
    a.set_title("gate satisfied, triple line manufactured")
    a.legend(loc="center right")
    a.annotate(f"$\\times{tr.tpb_ratio.iloc[-1]:.2f}$",
               xy=(tr.accepted.iloc[-1] / 1e3, tr.tpb_ratio.iloc[-1]),
               xytext=(-4, 3), textcoords="offset points",
               ha="right", color=TPB_C, fontsize=7.5)
    a.annotate(f"${(tr.s_ratio.iloc[-1]-1)*100:+.1f}\\%$",
               xy=(tr.accepted.iloc[-1] / 1e3, tr.s_ratio.iloc[-1]),
               xytext=(-4, -9), textcoords="offset points",
               ha="right", color=NI_C, fontsize=7.5)
    panel(a, "a")

    b = ax[1]
    cols = [TPB_C if v == 0 else NI_C for v in hs.dA]
    b.bar(hs.dA, hs.share * 100, width=1.35, color=cols, zorder=3)
    for v, sh in zip(hs.dA, hs.share):
        if sh > 0.04:
            b.text(v, sh * 100 + 1.6, f"{sh*100:.0f}", ha="center",
                   fontsize=7, color="0.25")
    b.set_xticks(sorted(hs.dA))
    b.set_xlabel(r"$\Delta A$ (exposed faces)")
    b.set_ylabel("share of accepted moves (\\%)")
    b.set_title("four in five moves cost nothing")
    b.set_ylim(0, 95)
    tidy(b, minor_x=False)
    panel(b, "b")

    fig.tight_layout(w_pad=2.0)
    save(fig, "fig8_gate_vs_tpb")


def figure8():
    cf = pd.read_csv(f"{P2}/o7_counterfactual.csv").set_index("case")
    order = ["all_moves", "contact_adjacent_only", "away_from_contact_only"]
    labels = ["all\nmoves", "only at the\nNi--YSZ contact",
              "only away from\nthe contact"]
    vals = [cf.loc[k, "tpb_ratio"] for k in order]

    fig, ax = plt.subplots(figsize=(COL_W, 2.15))
    ax.axhline(1.0, color="0.55", lw=0.6, zorder=1)
    bars = ax.bar(range(3), vals, 0.56,
                  color=["0.55", TPB_C, COLOR["coarse"]], zorder=3)
    for r, v in zip(bars, vals):
        ax.text(r.get_x() + r.get_width() / 2, v + 0.12, f"{v:.2f}",
                ha="center", fontsize=7.5)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("TPB / pristine TPB")
    ax.set_ylim(0, 6.1)
    ax.set_title("where the inflation comes from")
    tidy(ax, minor_x=False)
    fig.tight_layout()
    save(fig, "fig9_counterfactual")


if __name__ == "__main__":
    figure7()
    figure8()
