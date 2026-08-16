"""Every figure in the primer, derived from the actual FIB-SEM tomograms.

No schematics drawn from imagination and no reproduced third-party figures:
each panel is either a real slice out of the segmented volumes, or a plot of
numbers measured on them, or (for the athermal barrier) an arithmetic curve
evaluated at this dataset's own voxel sizes.

Palette follows the dataviz reference instance. Charts are single-axis by
construction -- where two quantities of different scale must be compared they
are indexed to their pristine values and drawn as dimensionless ratios.

Output: out/writeup/figs_primer/*.png
"""
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import scipy.ndimage as ndi
import tifffile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.damage2 import tpb_density_um2                      # noqa: E402
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE      # noqa: E402
from cmlib.io import label_histogram, slice_paths              # noqa: E402
from cmlib.percolation import percolating_mask                 # noqa: E402
from cmlib.phases import assign_labels                         # noqa: E402
from cmlib.roi import tile_rois                                # noqa: E402
from cmlib.seqgreedy import SeqGreedy, ST6, ST7                # noqa: E402

OUT = os.path.join(ROOT, "out", "writeup", "figs_primer")
os.makedirs(OUT, exist_ok=True)

# ---- palette (dataviz reference instance, light surface) -------------------
SURF = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
NI = "#2a78d6"        # slot 1 blue
YSZ = "#eb6834"       # slot 2 orange
PORE = "#eceae4"      # recessive, reads as empty
JUNCT = "#e34948"     # slot 8 red -- reserved here for the triple junction
AQUA = "#1baf7a"      # slot 3
VIOLET = "#4a3aa7"    # slot 7

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.facecolor": SURF,
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False,
})

SIDE = {"fine": 8.0, "medium": 8.0, "coarse": 12.0}
AXIS_P, CONN = 2, 6


def tidy(ax):
    ax.set_axisbelow(True)
    ax.tick_params(length=0)


def load_roi(grain, side=None):
    k = {s[2]: s for s in SAMPLES if s[3] == "pristine"}[grain]
    mapping = assign_labels(label_histogram(k[1])["counts"], ZENODO_LABEL_NOTE[k[0]])
    r = tile_rois(k[6], k[5], k[4], k[9], k[8], k[7],
                  side or SIDE[grain], max_rois=1)[0]
    ps = slice_paths(k[1])
    sh = (r["z1"] - r["z0"], r["y1"] - r["y0"], r["x1"] - r["x0"])
    ni = np.empty(sh, bool)
    ysz = np.empty(sh, bool)
    for i, z in enumerate(range(r["z0"], r["z1"])):
        a = tifffile.imread(ps[z])[r["y0"]:r["y1"], r["x0"]:r["x1"]]
        ni[i] = a == mapping["Ni"]
        ysz[i] = a == mapping["YSZ"]
    vox = float((k[9] * k[8] * k[7]) ** (1 / 3))
    return ni, ysz, vox


def rgb_slice(ni2d, ysz2d):
    """Three-phase false colour of one 2D slice."""
    img = np.zeros(ni2d.shape + (3,), float)
    for mask, hexc in ((ni2d, NI), (ysz2d, YSZ),
                       (~(ni2d | ysz2d), PORE)):
        c = np.array([int(hexc[i:i + 2], 16) / 255 for i in (1, 3, 5)])
        img[mask] = c
    return img


def scalebar(ax, vox_nm, npx, length_um=2.0, y=0.94):
    px = length_um * 1000.0 / vox_nm
    x0 = npx * 0.05
    ax.plot([x0, x0 + px], [npx * y, npx * y], color=INK, lw=2.5,
            solid_capstyle="butt")
    ax.text(x0 + px / 2, npx * y - npx * 0.025, f"{length_um:g} µm",
            ha="center", va="bottom", color=INK, fontsize=8)


def spec_surf(m):
    s = 0
    for ax_ in range(3):
        for sh in (1, -1):
            s += int((m & ~np.roll(m, sh, axis=ax_)).sum())
    return s / max(int(m.sum()), 1)


# =========================================================== FIGURE 1
def fig_series(cache):
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.5))
    for ax, grain in zip(axes, ("fine", "medium", "coarse")):
        ni, ysz, vox = cache[grain]
        n = 300
        z = ni.shape[0] // 2
        sl_ni, sl_ysz = ni[z, :n, :n], ysz[z, :n, :n]
        ax.imshow(rgb_slice(sl_ni, sl_ysz), interpolation="nearest")
        scalebar(ax, vox, n)
        ax.set_title(f"{grain}   ({vox:.0f} nm voxel)", color=INK)
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color(AXIS)
    fig.legend(handles=[Patch(facecolor=NI, label="nickel"),
                        Patch(facecolor=YSZ, label="YSZ"),
                        Patch(facecolor=PORE, label="pore")],
               loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("The same electrode, made three ways: real FIB-SEM slices",
                 color=INK, y=1.0)
    fig.tight_layout(rect=[0, 0.05, 1, 0.97])
    fig.savefig(f"{OUT}/fig1_series.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig1_series")


# =========================================================== FIGURE 2
def fig_what_is_tpb(cache):
    """Real slice, zoomed, with the triple-junction voxels marked."""
    ni, ysz, vox = cache["fine"]
    z = ni.shape[0] // 2
    # a 2D proxy for the junction: Ni voxels touching both YSZ and pore
    pore = ~(ni | ysz)
    st = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)
    a, b = 120, 200
    sl_ni = ni[z, a:b, a:b]; sl_ysz = ysz[z, a:b, a:b]; sl_po = pore[z, a:b, a:b]
    junc = (sl_ni & ndi.binary_dilation(sl_ysz, st) &
            ndi.binary_dilation(sl_po, st))

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.3))
    for ax, show in zip(axes, (False, True)):
        ax.imshow(rgb_slice(sl_ni, sl_ysz), interpolation="nearest")
        if show:
            ys, xs = np.nonzero(junc)
            ax.scatter(xs, ys, s=7, color=JUNCT, linewidths=0)
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color(AXIS)
    axes[0].set_title("Ni, YSZ, pore", color=INK)
    axes[1].set_title(f"triple junction marked  ({int(junc.sum())} sites "
                      f"in this slice)", color=INK)
    scalebar(axes[0], vox, b - a, length_um=0.5)
    fig.suptitle("Where the reaction happens: only on the red line",
                 color=INK)
    fig.text(0.5, -0.02, "The electrochemistry runs only where all three "
             "phases meet. In 3D this set is a curve; a slice cuts it in "
             "points.", ha="center", color=INK2, fontsize=8.5)
    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    fig.savefig(f"{OUT}/fig2_tpb.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig2_tpb")


# =========================================================== FIGURE 3
def fig_percolation(cache):
    """Spanning nickel vs already-disconnected nickel, in real data."""
    ni, ysz, vox = cache["coarse"]          # coarse has the most isolated Ni
    span = percolating_mask(ni, axis=AXIS_P, connectivity=CONN)
    iso = ni & ~span
    frac = iso.sum() / ni.sum()
    z = ni.shape[0] // 2
    n = 320
    sl_span, sl_iso = span[z, :n, :n], iso[z, :n, :n]

    img = np.ones(sl_span.shape + (3,))
    img[:] = np.array([int(PORE[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    img[ysz[z, :n, :n]] = [0.93, 0.91, 0.87]
    img[sl_span] = [int(NI[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    img[sl_iso] = [int(JUNCT[i:i + 2], 16) / 255 for i in (1, 3, 5)]

    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    ax.imshow(img, interpolation="nearest")
    scalebar(ax, vox, n)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(AXIS)
    ax.legend(handles=[Patch(facecolor=NI, label="nickel in the spanning network"),
                       Patch(facecolor=JUNCT, label="nickel already disconnected")],
              loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=1)
    ax.set_title(f"Real electrodes are not fully connected\n"
                 f"coarse anode: {frac*100:.1f}% of nickel is isolated "
                 f"before any degradation", color=INK, fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig3_percolation.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  fig3_percolation (isolated {frac*100:.1f}%)")
    return frac


# =========================================================== FIGURE 4
def fig_signature():
    """The measured anticorrelation -- and the fact that it is two-level."""
    anodes = ["fine", "medium", "coarse"]
    pspan = [0.680, 0.855, 0.947]
    preach = [0.857, 0.979, 0.942]
    tpb_ours = [0.799, 0.746, 0.590]
    tpb_pub = [0.7434, 0.5862, 0.6075]
    x = np.arange(3)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    ax = axes[0]
    ax.bar(x - 0.19, pspan, 0.34, color=NI, label="$P_{span}$")
    ax.bar(x + 0.19, preach, 0.34, color=AQUA, label="$P_{reach}$")
    for xi, (a, b) in enumerate(zip(pspan, preach)):
        ax.text(xi - 0.19, a + .012, f"{a:.3f}", ha="center", fontsize=8, color=INK2)
        ax.text(xi + 0.19, b + .012, f"{b:.3f}", ha="center", fontsize=8, color=INK2)
    ax.set_ylim(0, 1.13); ax.set_xticks(x); ax.set_xticklabels(anodes)
    ax.set_ylabel("nickel percolation retained")
    ax.set_title("Ni connectivity — the ordering flips\nwith the definition",
                 color=INK, fontsize=10)
    ax.legend(loc="lower right"); tidy(ax)

    ax = axes[1]
    ax.bar(x - 0.19, tpb_ours, 0.34, color=NI, label="this work")
    ax.bar(x + 0.19, tpb_pub, 0.34, color=YSZ, label="published")
    for xi, (a, b) in enumerate(zip(tpb_ours, tpb_pub)):
        ax.text(xi - 0.19, a + .012, f"{a:.3f}", ha="center", fontsize=8, color=INK2)
        ax.text(xi + 0.19, b + .012, f"{b:.3f}", ha="center", fontsize=8, color=INK2)
    ax.set_ylim(0, 1.0); ax.set_xticks(x); ax.set_xticklabels(anodes)
    ax.set_ylabel("TPB retained")
    ax.set_title("TPB — the published series is\nnon-monotone (coarse > medium)",
                 color=INK, fontsize=10)
    ax.legend(loc="lower right"); tidy(ax)
    fig.suptitle("Fine retains connectivity worst and reaction sites best. "
                 "Only that much is defensible.", color=INK, fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{OUT}/fig4_signature.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig4_signature")


# =========================================================== FIGURE 5
def fig_swap_identity(cache):
    """A real swap, taken out of the fine tomogram, with dA computed."""
    ni, ysz, vox = cache["fine"]
    nb = ndi.convolve(ni.astype(np.int8), ST6.astype(np.int8),
                      mode="constant", cval=0)
    interior = np.zeros_like(ni)
    interior[1:-1, 1:-1, 1:-1] = True
    surf = ni & (nb < 6) & interior
    front = (~ni) & (~ysz) & (nb >= 1) & interior
    rng = np.random.default_rng(0)
    # a convex source (few Ni neighbours) and a concave sink (many)
    src = np.argwhere(surf & (nb == 1))
    dst = np.argwhere(front & (nb == 4))
    a = src[rng.integers(len(src))]
    b = dst[rng.integers(len(dst))]
    nA, nB = int(nb[tuple(a)]), int(nb[tuple(b)])
    dA = 2 * (nA - nB)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2))
    for ax, c, lab, nn in ((axes[0], a, "source $a$ — convex", nA),
                           (axes[1], b, "sink $b$ — concave", nB)):
        z, y, x = c
        h = 11
        sl_ni = ni[z, y - h:y + h + 1, x - h:x + h + 1]
        sl_ysz = ysz[z, y - h:y + h + 1, x - h:x + h + 1]
        ax.imshow(rgb_slice(sl_ni, sl_ysz), interpolation="nearest")
        ax.add_patch(Rectangle((h - .5, h - .5), 1, 1, fill=False,
                               edgecolor=JUNCT, lw=2.2))
        ax.set_title(f"{lab}\n$n_N$ = {nn}", color=INK, fontsize=10)
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color(AXIS)
    fig.suptitle(f"A real swap:  $\\Delta A = 2({nA} - {nB}) = {dA}$ faces  "
                 f"— accepted", color=INK, fontsize=11)
    fig.text(0.5, -0.03,
             "Both voxels are taken from the fine tomogram. Moving nickel from "
             "the convex site to the concave one\nremoves "
             f"{abs(dA)} exposed faces. This is the whole of the "
             "surface-area criterion.",
             ha="center", color=INK2, fontsize=8.5)
    fig.tight_layout(rect=[0, 0.02, 1, 0.92])
    fig.savefig(f"{OUT}/fig5_swap.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  fig5_swap (dA={dA})")


# =========================================================== FIGURES 6+7+8
def fig_run(cache):
    """One instrumented run on the fine ROI: dA histogram, indexed
    trajectories, and the before/after slice."""
    ni0, ysz, vox = cache["fine"]
    s0 = spec_surf(ni0)
    t0 = tpb_density_um2(ni0, ysz, vox)
    n0 = int(ni0.sum())
    r0 = int(percolating_mask(ni0, axis=AXIS_P, connectivity=CONN).sum()) / n0

    op = SeqGreedy(ni0, ysz, seed=300)
    op.dA_log = []
    total = int(round(0.03 * op.n_surf0)) * 5
    checkpoints = np.unique(np.linspace(0, total, 13).astype(int))[1:]
    xs, s_r, t_r, r_r = [0], [1.0], [1.0], [1.0]
    done = 0
    tic = time.time()
    for cp in checkpoints:
        op.run(cp - done); done = cp
        s_r.append(spec_surf(op.ni) / s0)
        t_r.append(tpb_density_um2(op.ni, ysz, vox) / t0)
        r_r.append((int(percolating_mask(op.ni, axis=AXIS_P,
                                         connectivity=CONN).sum()) / n0) / r0)
        xs.append(op.accepted)
        print(f"    cp {op.accepted:>7}  S {s_r[-1]:.4f}  TPB {t_r[-1]:.3f} "
              f"  [{time.time()-tic:.0f}s]")
    dA = np.array(op.dA_log)

    # ---- FIGURE 6: the dA histogram
    fig, ax = plt.subplots(figsize=(6.6, 3.9))
    vals, counts = np.unique(dA, return_counts=True)
    cols = [JUNCT if v == 0 else NI for v in vals]
    ax.bar(vals, counts / counts.sum() * 100, width=1.4, color=cols)
    for v, c in zip(vals, counts):
        if c / counts.sum() > 0.03:
            ax.text(v, c / counts.sum() * 100 + 1.2,
                    f"{c/counts.sum()*100:.0f}%", ha="center",
                    fontsize=8.5, color=INK2)
    ax.set_xlabel("$\\Delta A$  (change in exposed nickel faces)")
    ax.set_ylabel("share of accepted moves  (%)")
    ax.set_title("Four in five accepted moves cost exactly nothing",
                 color=INK, fontsize=10.5)
    ax.set_xticks(sorted(vals))
    tidy(ax)
    ax.legend(handles=[Patch(facecolor=JUNCT, label="$\\Delta A = 0$ — unpriced by the gate"),
                       Patch(facecolor=NI, label="$\\Delta A < 0$ — area actually removed")],
              loc="upper left")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig6_plateau.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  fig6_plateau (neutral {100*(dA==0).mean():.1f}%)")

    # ---- FIGURE 7: the money figure (single axis, indexed)
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.axhline(1.0, color=AXIS, lw=1.2, zorder=1)
    ax.plot(xs, t_r, color=JUNCT, lw=2, marker="o", ms=4.5, zorder=3,
            label="TPB density")
    ax.plot(xs, s_r, color=NI, lw=2, marker="o", ms=4.5, zorder=3,
            label="specific surface area  (the gated quantity)")
    ax.set_yscale("log")
    ax.set_yticks([0.9, 1, 2, 3, 5])
    ax.set_yticklabels(["0.9×", "1×", "2×", "3×", "5×"])
    ax.set_xlabel("accepted moves")
    ax.set_ylabel("value ÷ pristine value")
    ax.annotate(f"TPB ×{t_r[-1]:.1f}", xy=(xs[-1], t_r[-1]),
                xytext=(-6, 8), textcoords="offset points",
                ha="right", color=JUNCT, fontsize=10, fontweight="bold")
    ax.annotate(f"area {(s_r[-1]-1)*100:+.1f}%  — gate satisfied",
                xy=(xs[-1], s_r[-1]), xytext=(-6, -16),
                textcoords="offset points", ha="right", color=NI,
                fontsize=10, fontweight="bold")
    ax.set_title("The validity check passes while the physics is destroyed\n"
                 "fine anode, 67 million voxels, volume conserved exactly",
                 color=INK, fontsize=11)
    ax.legend(loc="center right")
    tidy(ax)
    fig.text(0.5, -0.04, "Both quantities are divided by their pristine value, "
             "so they share one axis. The operator was required to reduce "
             "surface area.\nIt did. Nothing required it to preserve the "
             "three-phase junction, and it did not.",
             ha="center", color=INK2, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig7_money.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig7_money")

    # ---- FIGURE 8: before / after, real slice
    z = ni0.shape[0] // 2
    a, b = 150, 230
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.5))
    for ax, m, lab in ((axes[0], ni0, "pristine"),
                       (axes[1], op.ni, f"after {op.accepted:,} swaps")):
        ax.imshow(rgb_slice(m[z, a:b, a:b], ysz[z, a:b, a:b]),
                  interpolation="nearest")
        ax.set_title(lab, color=INK)
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color(AXIS)
    scalebar(axes[0], vox, b - a, length_um=0.5)
    fig.suptitle("What the operator actually does to the nickel",
                 color=INK, fontsize=11)
    fig.text(0.5, -0.03,
             "Nickel is not retracting into compact particles. It is "
             "speckling along the interface — and every new fleck that touches "
             "zirconia\nand pore adds triple line. That is the 5× rise, and it "
             "is the opposite of coarsening.",
             ha="center", color=INK2, fontsize=8.5)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    fig.savefig(f"{OUT}/fig8_before_after.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig8_before_after")
    del op


# =========================================================== FIGURE 9
def fig_strict():
    # Read from the committed tables rather than hard-coding: values baked into
    # a figure script go stale silently when the underlying run is repeated.
    import pandas as pd
    P2 = os.path.join(ROOT, "out", "project2")
    anodes = ["fine", "medium", "coarse"]
    ns = pd.read_csv(f"{P2}/o7_strict_vs_nonstrict.csv").set_index("anode")
    nonstrict = [float(ns.loc[a, "ratio_nonstrict"]) for a in anodes]
    strict = [float(ns.loc[a, "ratio_strict"]) for a in anodes]
    x = np.arange(3)
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    ax.axhline(1.0, color=AXIS, lw=1.2)
    ax.bar(x - 0.19, nonstrict, 0.34, color=JUNCT, label="$\\Delta A \\leq 0$  (neutral moves allowed)")
    ax.bar(x + 0.19, strict, 0.34, color=NI, label="$\\Delta A < 0$  (strictly area-lowering)")
    for xi, (a, b) in enumerate(zip(nonstrict, strict)):
        ax.text(xi - 0.19, a + .08, f"{a:.2f}×", ha="center", fontsize=8.5, color=INK2)
        ax.text(xi + 0.19, b + .08, f"{b:.2f}×", ha="center", fontsize=8.5, color=INK2)
    ax.text(2.44, 1.04, "no change", fontsize=8, color=MUTED, va="bottom", ha="right")
    ax.set_xticks(x); ax.set_xticklabels(anodes)
    ax.set_ylabel("TPB after ÷ TPB pristine")
    ax.set_ylim(0, 5.9)
    ax.set_title("Closing the loophole helps, and is not enough",
                 color=INK, fontsize=10.5)
    ax.legend(loc="upper right"); tidy(ax)
    fig.text(0.5, -0.04, "Forbidding the unpriced moves removes most of the "
             "artifact. What remains — 1.33× to 1.53× — is manufactured by "
             "moves that\nstrictly lower surface area. Area is the wrong "
             "invariant, not a leaky one.",
             ha="center", color=INK2, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig9_strict.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig9_strict")


# =========================================================== FIGURE 10
def fig_athermal():
    gamma = 2.0
    kT = 1.380649e-23 * (950 + 273.15)
    edges = np.linspace(0.3, 35, 400)
    ratio = gamma * (edges * 1e-9) ** 2 / kT
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    ax.plot(edges, ratio, color=NI, lw=2, zorder=3)
    ax.axhline(1.0, color=JUNCT, lw=1.6, ls="--", zorder=2)
    ax.text(0.42, 1.5, "$J_{NP} = k_BT$   — thermal fluctuation can act here",
            color=JUNCT, fontsize=8.5, va="bottom")
    ax.set_yscale("log"); ax.set_xscale("log")
    for e, lab in ((17.9, ""), (20.0, "this dataset"), (29.1, "")):
        r = gamma * (e * 1e-9) ** 2 / kT
        ax.plot([e], [r], "o", color=INK, ms=6, zorder=4)
        if lab:
            ax.annotate(f"{lab}\n{r:.1e} $k_BT$", xy=(e, r),
                        xytext=(-14, -46), textcoords="offset points",
                        ha="center", color=INK, fontsize=9,
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    ax.set_xlabel("voxel edge  (nm)")
    ax.set_ylabel("energy of one voxel face  ($k_BT$ at 950 °C)")
    ax.set_title("Why finite temperature cannot rescue the method",
                 color=INK, fontsize=10.5)
    tidy(ax)
    fig.text(0.5, -0.04,
             "A voxel face at tomographic resolution carries ~10⁴–10⁵ $k_BT$. "
             "The Boltzmann factor for an area-raising move underflows to "
             "zero.\nThermal activation works at the atomic scale, which this "
             "lattice has discarded.",
             ha="center", color=INK2, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig10_athermal.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig10_athermal")


def main():
    print("loading tomograms…")
    cache = {}
    for g in ("fine", "medium", "coarse"):
        t = time.time()
        cache[g] = load_roi(g)
        print(f"  {g}: {cache[g][0].shape}  {time.time()-t:.0f}s")
    print("rendering…")
    fig_series(cache)
    fig_what_is_tpb(cache)
    fig_percolation(cache)
    fig_signature()
    fig_swap_identity(cache)
    fig_strict()
    fig_athermal()
    del cache["medium"], cache["coarse"]
    fig_run(cache)
    print("done ->", OUT)


if __name__ == "__main__":
    main()
