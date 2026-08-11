"""Manuscript figures for the Ni-YSZ voxel-operator paper.

Every value plotted is read from a committed CSV under out/project2 or
out/phase6. Nothing is hard-coded except axis limits and annotation text.
Run `python scripts/project2/o5v2_transcribe.py` first if
out/project2/o5v2_area_barrier.csv is absent.

Outputs six PDFs (vector, for the LaTeX build) and six PNGs (for the markdown
preview) into out/writeup/figs.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotstyle import (ANODES, COLOR, MARKER, SYNTH, ACCENT, COL_W, FULL_W,
                       use, panel, tidy)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P2 = os.path.join(ROOT, "out", "project2")
P6 = os.path.join(ROOT, "out", "phase6")
FIG = os.path.join(ROOT, "out", "writeup", "figs")
os.makedirs(FIG, exist_ok=True)
use()


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig)
    print("  wrote", name)


# --------------------------------------------------------------------------
# Figure 1 -- the measured signature the operators were built to reproduce.
# --------------------------------------------------------------------------
t = pd.read_csv(os.path.join(P6, "phase6_comparison_table.csv"), index_col=0)
x = np.arange(3)
w = 0.26

fig, ax = plt.subplots(1, 2, figsize=(FULL_W, 2.35))

series_ni = [("P_span_retained", "$P_{\\mathrm{span}}$ (this work)", "#08306B"),
             ("P_reach_retained", "$P_{\\mathrm{reach}}$ (this work)", "#4292C6"),
             ("P_published_retained", "published $P$", "#C6DBEF")]
for k, (col, lab, c) in enumerate(series_ni):
    ax[0].bar(x + (k - 1) * w, t.loc[ANODES, col], w, color=c,
              edgecolor="k", linewidth=0.4, label=lab)
ax[0].set_ylim(0, 1.05)
ax[0].set_ylabel("Ni percolation retained")
ax[0].set_title("Fine retains Ni percolation worst")
ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3,
             columnspacing=1.2)

series_tpb = [("tpb_retained", "this work", "#08306B"),
              ("tpb_retained_published", "published", "#C6DBEF")]
for k, (col, lab, c) in enumerate(series_tpb):
    ax[1].bar(x + (k - 0.5) * w, t.loc[ANODES, col], w, color=c,
              edgecolor="k", linewidth=0.4, label=lab)
ax[1].set_ylim(0, 1.05)
ax[1].set_ylabel("TPB density retained")
ax[1].set_title("Fine retains TPB best")
ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2,
             columnspacing=1.2)
# bracket the non-monotone pair in the published series
xa, xb, yb = 1 + 0.5 * w, 2 + 0.5 * w, 0.72
ax[1].plot([xa, xa, xb, xb], [0.62, yb, yb, 0.64], color=ACCENT, lw=0.7,
           clip_on=False)
ax[1].text((xa + xb) / 2, yb + 0.02, "non-monotone", fontsize=6.5,
           color=ACCENT, ha="center", va="bottom")

for a in ax:
    a.set_xticks(x)
    a.set_xticklabels(ANODES)
    tidy(a, minor_x=False)
panel(ax[0], "a", dx=-0.13)
panel(ax[1], "b", dx=-0.13)
fig.tight_layout(w_pad=2.0)
save(fig, "fig1_measured_signature")


# --------------------------------------------------------------------------
# Figure 2 -- the outcome metric rewrites its own denominator.
# --------------------------------------------------------------------------
val = pd.read_csv(os.path.join(P2, "c1real_o6_validity.csv"))
rni = pd.read_csv(os.path.join(P2, "c1real_rni_gate.csv"))

fig, ax = plt.subplots(1, 2, figsize=(FULL_W, 2.4))

for a in ANODES:
    s = val[val.anode == a].sort_values("n_rounds")
    n = [0] + list(s.n_rounds)
    p = [s.pristine_P_span.iloc[0]] + list(s.P_span)
    ax[0].plot(n, p, MARKER[a] + "-", color=COLOR[a], label=a, clip_on=False)
ax[0].axhline(1.0, color="k", lw=0.6, ls=(0, (1, 2)))
ax[0].set_ylim(0.86, 1.012)
ax[0].set_xlim(0, 5.2)
ax[0].set_xlabel("damage rounds $n$")
ax[0].set_ylabel("$P_{\\mathrm{span}}$")
ax[0].set_title("$P_{\\mathrm{span}}$: pinned at unity by pruning")
ax[0].text(2.5, 0.995, "$P_{\\mathrm{span}}=1$ identically after pruning",
           fontsize=6.5, ha="center", va="top", color=ACCENT)
ax[0].legend(loc="lower right")

for a in ANODES:
    s = rni[rni.anode == a].sort_values("n_rounds")
    n = [0] + list(s.n_rounds)
    r = [s.R_Ni_0.iloc[0]] + list(s.R_Ni)
    ax[1].plot(n, r, MARKER[a] + "-", color=COLOR[a], label=a, clip_on=False)
ax[1].set_ylim(0.5, 1.02)
ax[1].set_xlim(0, 8.3)
ax[1].set_xlabel("damage rounds $n$")
ax[1].set_ylabel("$R_{\\mathrm{Ni}}$")
ax[1].set_title("$R_{\\mathrm{Ni}}$: monotone, pruning-invariant")

for a in ax:
    tidy(a)
panel(ax[0], "a", dx=-0.14)
panel(ax[1], "b", dx=-0.14)
fig.tight_layout(w_pad=2.0)
save(fig, "fig2_metric")


# --------------------------------------------------------------------------
# Figure 3 -- voxel erosion manufactures TPB before destroying it.
# --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(COL_W, 2.5))
for a in ANODES:
    s = rni[rni.anode == a].sort_values("n_rounds")
    ax.plot([0] + list(s.n_rounds), [1.0] + list(s.tpb_um2 / s.tpb_pristine),
            MARKER[a] + "-", color=COLOR[a], label=a)
ax.axhline(1.0, color="k", lw=0.6, ls=(0, (1, 2)))
ax.set_yscale("log")
ax.set_xlim(0, 8.3)
ax.set_xlabel("damage rounds $n$")
ax.set_ylabel("TPB / TPB$_{\\mathrm{pristine}}$")
ax.set_ylim(3e-2, 40)
ax.annotate("manufactured", xy=(0.55, 4.0), xytext=(1.5, 25),
            fontsize=6.5, color=ACCENT,
            arrowprops=dict(arrowstyle="->", lw=0.6, color=ACCENT))
ax.annotate("destroyed", xy=(7.2, 0.25), xytext=(4.9, 0.09),
            fontsize=6.5, color=ACCENT,
            arrowprops=dict(arrowstyle="->", lw=0.6, color=ACCENT))
ax.legend(loc="lower left", bbox_to_anchor=(-0.02, -0.02))
tidy(ax, minor_y=False)
fig.tight_layout()
save(fig, "fig3_tpb")


# --------------------------------------------------------------------------
# Figure 4 -- lattice cuts at a full cross-section; real networks do not.
# --------------------------------------------------------------------------
lat = pd.read_csv(os.path.join(P2, "audit_ni_vulnerability.csv"))
real = pd.read_csv(os.path.join(P2, "c1real_results.csv")).drop_duplicates(
    ["anode", "roi"])

fig, ax = plt.subplots(figsize=(COL_W, 2.5))
rng = np.random.default_rng(0)
for i, a in enumerate(ANODES):
    l = lat[lat.analog == a]
    r = real[real.anode == a].mincut_frac.dropna().values
    # all five lattice seeds return the identical fraction, so one marker
    assert l.frac_mincut.nunique() == 1
    ax.scatter([i - 0.19], [l.frac_mincut.iloc[0]], s=26, facecolor="none",
               edgecolor=SYNTH, linewidth=0.9, marker="s",
               label="synthetic lattice" if i == 0 else None)
    ax.annotate(f"${int(round(l.mincut_edges.iloc[0] ** 0.5))}^2$",
                (i - 0.19, l.frac_mincut.iloc[0]),
                textcoords="offset points", xytext=(-14, -3), fontsize=6.5,
                color=SYNTH)
    ax.scatter(i + 0.19 + rng.uniform(-.05, .05, len(r)), r, s=22,
               color=COLOR[a], marker=MARKER[a],
               label="real ROIs" if i == 0 else None)
ax.set_xticks(range(3))
ax.set_xticklabels(ANODES)
ax.set_xlim(-0.6, 2.6)
ax.set_yscale("log")
ax.set_ylim(3e-3, 0.13)
ax.set_ylabel("minimum cut, fraction of throats")
ax.text(0.5, 0.088, "5 seeds, zero variance", fontsize=6.5, color=SYNTH)
ax.legend(loc="lower left")
fig.tight_layout()
save(fig, "fig4_mincut")


# --------------------------------------------------------------------------
# Figure 5 -- the simulation inverts the measurement, and neither candidate
# predictor accounts for it.
# --------------------------------------------------------------------------
res = pd.read_csv(os.path.join(P2, "c1real_results.csv"))
tr = res[res.threshold == 0.50]
means = tr.groupby("anode").transition.mean().loc[ANODES]
per_roi = tr.groupby(["anode", "roi"]).transition.mean()

fig, ax = plt.subplots(1, 3, figsize=(FULL_W, 2.25))

ax[0].bar(range(3), t.loc[ANODES, "P_span_retained"], 0.55,
          color=[COLOR[a] for a in ANODES], edgecolor="k", linewidth=0.4)
ax[0].set_ylim(0, 1.05)
ax[0].set_ylabel("Ni percolation retained")
ax[0].set_title("measurement: fine worst")

for i, a in enumerate(ANODES):
    y = tr[tr.anode == a].transition.values
    ax[1].scatter(i + rng.uniform(-.13, .13, len(y)), y, s=11,
                  color=COLOR[a], marker=MARKER[a], alpha=0.75, zorder=3)
    ax[1].plot([i - 0.28, i + 0.28], [means[a]] * 2, color="k", lw=1.3,
               zorder=4)
ax[1].set_ylim(7.1, 10.9)
ax[1].set_ylabel("erosion rounds to $R_{\\mathrm{Ni}}<0.5$")
ax[1].set_title("simulation: fine last")

for a in ANODES:
    s = tr[tr.anode == a].groupby("roi").agg(
        {"spec_surface": "first", "transition": "mean"})
    ax[2].scatter(s.spec_surface, s.transition, s=22, color=COLOR[a],
                  marker=MARKER[a], label=a)
ax[2].set_xlabel("specific Ni surface area (vox$^{-1}$)")
ax[2].set_ylabel("erosion rounds")
ax[2].set_title("no rate dependence")
ax[2].set_xlim(0.108, 0.172)
ax[2].set_ylim(8.9, 10.15)
ax[2].text(0.04, 0.93, r"$\rho=+0.018$", transform=ax[2].transAxes,
           ha="left", va="top", fontsize=7)
ax[2].legend(loc="lower right", ncol=3, columnspacing=0.7,
             handletextpad=0.3, bbox_to_anchor=(1.03, -0.03))

for a in ax[:2]:
    a.set_xticks(range(3))
    a.set_xticklabels(ANODES)
for a in ax:
    tidy(a, minor_x=(a is ax[2]))
for i, L in enumerate("abc"):
    panel(ax[i], L, dx=-0.24)
fig.tight_layout(w_pad=1.6)
save(fig, "fig5_reversal")


# --------------------------------------------------------------------------
# Figure 6 -- the area barrier that brackets every single-swap operator.
# --------------------------------------------------------------------------
bar = pd.read_csv(os.path.join(P2, "o5v2_area_barrier.csv"))
cur6 = bar[(bar.operator == "curvature_ranked") & (bar.stencil == 6)].sort_values("n_rounds")
cur26 = bar[(bar.operator == "curvature_ranked") & (bar.stencil == 26)].sort_values("n_rounds")
grd = bar[bar.operator == "greedy_area"].sort_values("n_rounds")
S0 = float(grd.S_spec.iloc[0])

fig, ax = plt.subplots(1, 3, figsize=(FULL_W, 2.35))

# (a) schematic of the exact swap identity on a 6-connected lattice
sch = ax[0]
occ = {(i, j) for i in range(1, 7) for j in (1, 2)} - {(5, 2)}
occ.add((2, 3))
for (i, j) in occ:
    sch.add_patch(Rectangle((i, j), 1, 1, facecolor="#DEEBF7",
                            edgecolor="#9ECAE1", linewidth=0.5))
sch.add_patch(Rectangle((2, 3), 1, 1, facecolor="#FDDBC7",
                        edgecolor=ACCENT, linewidth=1.1))
sch.add_patch(Rectangle((5, 2), 1, 1, facecolor="none",
                        edgecolor=ACCENT, linewidth=1.1,
                        linestyle=(0, (2, 1.5))))
sch.annotate("", xy=(5.4, 2.9), xytext=(2.6, 3.9),
             arrowprops=dict(arrowstyle="->", lw=0.9, color=ACCENT,
                             connectionstyle="arc3,rad=0.28"))
sch.text(2.5, 4.15, "$a$: nb $=1$", ha="center", fontsize=7.5, color=ACCENT)
sch.text(5.5, 3.25, "$b$: nb $=3$", ha="center", fontsize=7.5, color=ACCENT)
sch.text(3.5, 0.25,
         r"$\Delta A = 2\,[\,\mathrm{nb}(a)-\mathrm{nb}(b)\,] = -4$",
         ha="center", fontsize=8)
sch.set_xlim(0.6, 7.4)
sch.set_ylim(-0.1, 4.7)
sch.set_aspect("equal")
sch.axis("off")

ax[1].plot(cur26.n_rounds, cur26.neck_voxels, "o-", color="#08306B",
           label="curvature-ranked, 26-conn")
ax[1].plot(cur6.n_rounds, cur6.neck_voxels, "s--", color="#6BAED6",
           label="curvature-ranked, 6-conn")
ax[1].plot(grd.n_rounds, grd.neck_voxels, "^-", color=ACCENT,
           label=r"greedy $\Delta A\leq0$")
ax[1].set_xlabel("damage rounds $n$")
ax[1].set_ylabel("neck volume (voxels)")
ax[1].set_ylim(-3, 70)
ax[1].set_title("neck thinning")
ax[1].legend(loc="lower left", bbox_to_anchor=(-0.02, -0.03))

ax[2].axhspan(S0, 0.4545, color=ACCENT, alpha=0.10, lw=0)
ax[2].axhline(S0, color="k", lw=0.6, ls=(0, (1, 2)))
ax[2].plot(cur26.n_rounds, cur26.S_spec, "o-", color="#08306B")
ax[2].plot(cur6.n_rounds, cur6.S_spec, "s--", color="#6BAED6")
ax[2].plot(grd.n_rounds, grd.S_spec, "^-", color=ACCENT)
ax[2].set_xlim(-0.2, 5.4)
ax[2].set_ylim(0.4405, 0.4545)
ax[2].set_xlabel("damage rounds $n$")
ax[2].set_ylabel("$S_{\\mathrm{spec}}$")
ax[2].set_title("surface area")
ax[2].text(2.7, 0.45325, "forbidden by the validity gate", fontsize=6.5,
           color=ACCENT, ha="center")
ax[2].annotate("greedy: no move accepted", xy=(4.0, S0), xytext=(2.15, 0.4478),
               fontsize=6.5, color=ACCENT,
               arrowprops=dict(arrowstyle="->", lw=0.6, color=ACCENT))

for a in ax[1:]:
    tidy(a)
for i, L in enumerate("abc"):
    panel(ax[i], L, dx=-0.22 if i else 0.0, dy=1.09 if i else 1.02)
fig.tight_layout(w_pad=1.8)
save(fig, "fig6_area_barrier")

print("done ->", FIG)
