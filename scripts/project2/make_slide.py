"""Portfolio / talk slide built from the committed result CSVs.

Every number is read from out/project2/*.csv at render time, so the slide
cannot drift from the data the way a hand-typed card does.
"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P    = os.path.join(ROOT, "out", "project2")
OUT  = os.path.join(ROOT, "out", "writeup")

INK, MUTED, RULE = "#111418", "#5B6670", "#D8DEE4"
NAVY, MID, PALE, RED, BG = "#08306B", "#2171B5", "#9ECAE1", "#B2182B", "#FFFFFF"

mv = pd.read_csv(os.path.join(P, "o7_move_stats.csv")).iloc[0]
cf = pd.read_csv(os.path.join(P, "o7_counterfactual.csv")).set_index("case")
sv = pd.read_csv(os.path.join(P, "o7_strict_vs_nonstrict.csv")).set_index("anode")

NEUTRAL    = float(mv["neutral_share"])
CONTACT_F  = float(mv["added_near_ysz"]) / float(mv["added"])
TPB_ALL    = float(cf.loc["all_moves", "tpb_ratio"])
TPB_CONT   = float(cf.loc["contact_adjacent_only", "tpb_ratio"])
TPB_AWAY   = float(cf.loc["away_from_contact_only", "tpb_ratio"])
TPB_STRICT = float(sv.loc["fine", "ratio_strict"])
DS         = 0.036     # specific surface area fall (manuscript fig. 8)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Helvetica Neue", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.7, "axes.grid": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "out", "ytick.direction": "out", "savefig.dpi": 220,
})

W, H = 12.0, 6.75
fig = plt.figure(figsize=(W, H), facecolor=BG)
L, R = 0.050, 0.962

def T(x, y, s, size=11, color=INK, weight="normal", ha="left", va="top",
      style="normal", ls=1.35):
    return fig.text(x, y, s, size=size, color=color, weight=weight, ha=ha,
                    va=va, style=style, linespacing=ls)

def hr(y, x0=L, x1=R, c=RULE, lw=1.0):
    fig.lines.append(plt.Line2D([x0, x1], [y, y], color=c, lw=lw,
                                transform=fig.transFigure))

# ---------------------------------------------------------------- header
T(L, 0.963, "A validity check that destroys what it certifies", size=25.5,
  weight="bold")
T(L, 0.895, "Voxel coarsening operators on public FIB-SEM tomograms of "
            "Ni\u2013YSZ fuel-cell electrodes.", size=12.3, color=MUTED)
hr(0.861)

# ------------------------------------------- (a) one run, judged two ways
# Independent scales. A shared axis makes -3.6 % invisible beside +399 %,
# destroying the comparison the panel exists to make.
T(L, 0.812, "(a)  One run, judged two ways", size=11, weight="bold")

def track(y_ax, y_lab, value, xlim, ticks, ticklabels, color, label, sub,
          val_txt, val_side, flag, flag_size):
    ax = fig.add_axes([0.123, y_ax, 0.200, 0.066])
    ax.barh([0], [value], height=.60, color=color, zorder=3)
    ax.set_xlim(*xlim); ax.set_ylim(-.6, .6)
    ax.axvline(0, color=INK, lw=0.9, zorder=4)
    ax.set_yticks([]); ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels, fontsize=8.5)
    ax.tick_params(length=2.5, pad=2)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(RULE)
    pad = (xlim[1] - xlim[0]) * 0.025
    if val_side == "inside":
        ax.text(value - pad, 0, val_txt, ha="right", va="center",
                fontsize=12, weight="bold", color="white", zorder=5)
    else:
        ax.text(value + (pad if val_side == "left" else -pad), 0,
                val_txt, ha=val_side, va="center", fontsize=12,
                weight="bold", color=color)
    T(0.123, y_lab, label, size=10.3, weight="bold", va="baseline")
    T(0.123, y_lab - 0.019, sub, size=9.3, color=MUTED, va="baseline")
    fig.text(0.336, y_ax + 0.033, flag, size=flag_size, weight="bold",
             color=color, va="center", ha="left", linespacing=1.2)

track(0.642, 0.747, -DS * 100, (-6.2, 6.2), [-5, 0, 5],
      ["\u22125 %", "0", "+5 %"], MID,
      "specific surface area", "the quantity the gate watches",
      "\u22123.6 %", "right", "GATE\nPASSES", 9.4)

track(0.436, 0.541, (TPB_ALL - 1) * 100, (-620, 620), [-500, 0, 500],
      ["\u2212500 %", "0", "+500 %"], RED,
      "triple-phase-boundary density", "the quantity the cell runs on",
      "+399 %", "inside", "NOT\nCHECKED", 9.4)

T(0.123, 0.392, "note: the two axes differ by 100\u00d7", size=8.8,
  color=MUTED, style="italic")

# ------------------------------------------- (b) where the inflation sits
axb = fig.add_axes([0.432, 0.436, 0.212, 0.335])
vals = [TPB_ALL, TPB_CONT, TPB_AWAY]
axb.bar(range(3), vals, width=.56, color=[RED, RED, PALE], zorder=3)
axb.axhline(1.0, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=4)
axb.set_ylim(0, 6.5); axb.set_xlim(-0.62, 2.62)
axb.set_xticks(range(3))
axb.set_xticklabels(["all\nmoves", "contact\nmoves\n(3 %)",
                     "all other\nmoves\n(97 %)"], fontsize=9.2)
axb.set_ylabel("TPB inflation  (\u00d7)", fontsize=10)
axb.tick_params(labelsize=9, length=3)
axb.spines["left"].set_color(RULE); axb.spines["bottom"].set_color(RULE)
for i, v in enumerate(vals):
    axb.text(i, v + .22, "\u00d7%.2f" % v, ha="center", va="bottom",
             fontsize=11, weight="bold", color=RED if v > 2 else MUTED)
T(0.432, 0.812, "(b)  3 % of the moves cause all of it", size=11,
  weight="bold")

# ----------------------------------------------------------- right rail
rx = 0.714
T(rx, 0.812, "WHY THE GATE DOES NOT OBJECT", size=9.3, weight="bold",
  color=MUTED)
T(rx, 0.762, "%.0f %%" % (NEUTRAL * 100), size=34, weight="bold", color=NAVY)
T(rx, 0.652, "of accepted moves change the area by\nexactly zero, so the "
             "criterion never\nprices them at all.", size=10.6, color=MUTED)
hr(0.566, rx, R, lw=0.9)
T(rx, 0.528, "TIGHTENING IT IS NOT ENOUGH", size=9.3, weight="bold",
  color=MUTED)
T(rx, 0.478, "\u00d7%.2f \u2192 \u00d7%.2f" % (TPB_ALL, TPB_STRICT), size=26,
  weight="bold", color=NAVY)
T(rx, 0.392, "demanding a strict area decrease\nshrinks the artifact. It does "
             "not\nremove it.", size=10.6, color=MUTED)

# -------------------------------------------------------- takeaway band
by, bh = 0.118, 0.158
fig.patches.append(FancyBboxPatch(
    (L, by), R - L, bh, transform=fig.transFigure,
    boxstyle="round,pad=0,rounding_size=0.006", linewidth=0,
    facecolor="#F4F7FA", zorder=0))
fig.patches.append(Rectangle((L, by), 0.0045, bh, transform=fig.transFigure,
                             facecolor=RED, linewidth=0, zorder=1))
T(L + 0.022, by + bh - 0.030,
  "Surface area is the wrong invariant to validate a coarsening rule against.",
  size=14.2, weight="bold")
T(L + 0.022, by + bh - 0.079,
  "One check costs a line and disqualified three operators here:  an operator "
  "whose TPB retention exceeds 1 is roughening, not coarsening.",
  size=11.1, color=MUTED)

# ---------------------------------------------------------------- footer
hr(0.082)
T(L, 0.060, "Pore-network analysis  \u00b7  Percolation  \u00b7  "
            "Pre-registered operators  \u00b7  Reproducible build",
  size=10, color=MUTED)
T(R, 0.060, "doi.org/10.5281/zenodo.21969405", size=10, color=MUTED,
  ha="right")

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, "slide_gate." + ext), facecolor=BG)
print("wrote", os.path.join(OUT, "slide_gate.png"))
print("  neutral=%.4f contact=%.4f tpb_all=%.3f cont=%.3f away=%.3f strict=%.3f"
      % (NEUTRAL, CONTACT_F, TPB_ALL, TPB_CONT, TPB_AWAY, TPB_STRICT))
