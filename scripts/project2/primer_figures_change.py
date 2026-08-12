"""Rework of the before/after panel, plus a tick fix on the trajectory figure.

The first attempt zoomed on an arbitrary window and the two panels looked
identical, because 2.6e5 swaps touch only ~0.4% of a 67 Mvoxel volume. Rather
than assert a description of the change, this script measures it: where the
changed voxels are, whether they are isolated or contiguous, and what they are
adjacent to. The caption is then written from the measurement.
"""
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Patch
import scipy.ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.damage2 import tpb_density_um2                      # noqa: E402
from cmlib.percolation import percolating_mask                 # noqa: E402
from cmlib.seqgreedy import SeqGreedy, ST7                     # noqa: E402

import primer_figures as pf                                    # noqa: E402

OUT = pf.OUT
NI, YSZ, PORE, JUNCT = pf.NI, pf.YSZ, pf.PORE, pf.JUNCT
INK, INK2, MUTED, AXIS = pf.INK, pf.INK2, pf.MUTED, pf.AXIS
CACHE = os.path.join(ROOT, "out", "project2", "_cache", "primer_run.npz")


def main():
    ni0, ysz, vox = pf.load_roi("fine")
    s0 = pf.spec_surf(ni0)
    t0 = tpb_density_um2(ni0, ysz, vox)

    op = SeqGreedy(ni0, ysz, seed=300)
    total = int(round(0.03 * op.n_surf0)) * 5
    checkpoints = np.unique(np.linspace(0, total, 13).astype(int))[1:]
    xs, s_r, t_r = [0], [1.0], [1.0]
    done = 0
    tic = time.time()
    for cp in checkpoints:
        op.run(cp - done); done = cp
        s_r.append(pf.spec_surf(op.ni) / s0)
        t_r.append(tpb_density_um2(op.ni, ysz, vox) / t0)
        xs.append(op.accepted)
    print(f"trajectory done in {time.time()-tic:.0f}s; "
          f"S {s_r[-1]:.4f}  TPB {t_r[-1]:.3f}")

    ni1 = op.ni
    added = ni1 & ~ni0
    removed = ni0 & ~ni1
    changed = added | removed

    # ---- measure, do not assume ------------------------------------------
    lab, n_cl = ndi.label(added, structure=ST7)
    sizes = np.bincount(lab.ravel())[1:]
    pore0 = ~(ni0 | ysz)
    add_near_ysz = int((added & ndi.binary_dilation(ysz, ST7)).sum())
    print(f"  added {added.sum():,}  removed {removed.sum():,}")
    print(f"  added clusters: {n_cl:,}   median size {np.median(sizes):.0f} "
          f"  singletons {(sizes==1).sum():,} ({(sizes==1).mean()*100:.0f}%)")
    print(f"  added voxels adjacent to YSZ: {add_near_ysz:,} "
          f"({add_near_ysz/added.sum()*100:.0f}%)")

    # pick the slice and window with the most change
    per_slice = changed.sum(axis=(1, 2))
    z = int(np.argmax(per_slice))
    cs = changed[z]
    k = 60
    ker = np.ones((k, k))
    dens = ndi.convolve(cs.astype(float), ker, mode="constant")
    yc, xc = np.unravel_index(int(np.argmax(dens)), dens.shape)
    y0 = int(np.clip(yc - k // 2, 0, cs.shape[0] - k))
    x0 = int(np.clip(xc - k // 2, 0, cs.shape[1] - k))
    print(f"  busiest slice z={z} ({per_slice[z]} changes), window "
          f"({y0},{x0}) {k}x{k} with {cs[y0:y0+k, x0:x0+k].sum()} changes")

    np.savez_compressed(
        CACHE, xs=xs, s_r=s_r, t_r=t_r,
        slice_ni0=ni0[z, y0:y0+k, x0:x0+k], slice_ni1=ni1[z, y0:y0+k, x0:x0+k],
        slice_ysz=ysz[z, y0:y0+k, x0:x0+k], vox=vox, accepted=op.accepted,
        singleton_frac=(sizes == 1).mean(), near_ysz=add_near_ysz/added.sum())

    # ---- FIGURE 7 redo, minor ticks suppressed ---------------------------
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.axhline(1.0, color=AXIS, lw=1.2, zorder=1)
    ax.plot(xs, t_r, color=JUNCT, lw=2, marker="o", ms=4.5, zorder=3,
            label="TPB density")
    ax.plot(xs, s_r, color=NI, lw=2, marker="o", ms=4.5, zorder=3,
            label="specific surface area  (the gated quantity)")
    ax.set_yscale("log")
    ax.set_yticks([0.9, 1, 2, 3, 4, 5])
    ax.set_yticklabels(["0.9×", "1×", "2×", "3×", "4×", "5×"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("accepted moves")
    ax.set_ylabel("value ÷ pristine value")
    ax.annotate(f"TPB ×{t_r[-1]:.1f}", xy=(xs[-1], t_r[-1]),
                xytext=(-8, 6), textcoords="offset points", ha="right",
                color=JUNCT, fontsize=10.5, fontweight="bold")
    ax.annotate(f"area {(s_r[-1]-1)*100:+.1f}%  — gate satisfied",
                xy=(xs[-1], s_r[-1]), xytext=(-8, 10),
                textcoords="offset points", ha="right", color=NI,
                fontsize=10.5, fontweight="bold")
    ax.set_title("The validity check passes while the physics is destroyed\n"
                 "fine anode, 67 million voxels, volume conserved exactly",
                 color=INK, fontsize=11)
    ax.legend(loc="center right")
    pf.tidy(ax)
    fig.text(0.5, -0.04, "Both quantities are divided by their pristine value, "
             "so they share one axis. The operator was required to reduce "
             "surface area.\nIt did. Nothing required it to preserve the "
             "three-phase junction, and it did not.",
             ha="center", color=INK2, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig7_money.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig7_money (fixed)")

    # ---- FIGURE 8 redo: three panels, busiest window ---------------------
    a0, a1 = ni0[z, y0:y0+k, x0:x0+k], ni1[z, y0:y0+k, x0:x0+k]
    yz = ysz[z, y0:y0+k, x0:x0+k]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.4))
    axes[0].imshow(pf.rgb_slice(a0, yz), interpolation="nearest")
    axes[0].set_title("pristine", color=INK)
    axes[1].imshow(pf.rgb_slice(a1, yz), interpolation="nearest")
    axes[1].set_title(f"after {op.accepted:,} swaps", color=INK)

    diff = np.ones(a0.shape + (3,))
    diff[:] = [0.97, 0.97, 0.96]
    diff[yz] = [0.91, 0.89, 0.85]
    diff[a0 & a1] = [0.80, 0.85, 0.92]
    diff[a1 & ~a0] = [int(JUNCT[i:i+2], 16)/255 for i in (1, 3, 5)]
    diff[a0 & ~a1] = [int(NI[i:i+2], 16)/255 for i in (1, 3, 5)]
    axes[2].imshow(diff, interpolation="nearest")
    axes[2].set_title("what moved", color=INK)
    axes[2].legend(handles=[
        Patch(facecolor=JUNCT, label="nickel added"),
        Patch(facecolor=NI, label="nickel removed"),
        Patch(facecolor=[0.80, 0.85, 0.92], label="nickel unchanged")],
        loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=3, fontsize=8)

    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color(AXIS)
    pf.scalebar(axes[0], vox, k, length_um=0.5)
    fig.suptitle("What the operator actually does to the nickel",
                 color=INK, fontsize=11)
    fig.tight_layout(rect=[0, 0.04, 1, 0.94])
    fig.savefig(f"{OUT}/fig8_before_after.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig8_before_after (rebuilt on the busiest window)")


if __name__ == "__main__":
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main()
