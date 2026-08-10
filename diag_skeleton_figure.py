"""Figure for the Phase-3 gate failure: skeleton is a CURVE for the fine anode
but a SHEET for the medium and coarse anodes."""
from __future__ import annotations

import os
import sys

import matplotlib
import numpy as np
from scipy import ndimage as ndi
from skan import csr
from skimage.morphology import skeletonize

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, load_subvolume  # noqa: E402
from cmlib.percolation import percolating_mask  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "phase3")
NEIGH = ndi.generate_binary_structure(3, 3)
NEIGH[1, 1, 1] = False

fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
keys = ["fine_pre", "medium_pre", "coarse_pre"]
cols = ["C0", "C1", "C2"]

for j, (key, col) in enumerate(zip(keys, cols)):
    _, folder, *_ , vx, vy, vz = [s for s in SAMPLES if s[0] == key][0]
    _, folder, grain, state, nx_, ny_, nz_, vx, vy, vz = \
        [s for s in SAMPLES if s[0] == key][0]
    spacing = (vz, vy, vx)
    counts = label_histogram(folder)["counts"]
    mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
    r = tile_rois(nz_, ny_, nx_, vz, vy, vx, 8.0, 1)[0]
    vol = load_subvolume(folder, r["z0"], r["z1"], r["y0"], r["y1"],
                         r["x0"], r["x1"])
    ni = vol == mapping["Ni"]
    del vol
    mask = percolating_mask(ni, axis=2, connectivity=6)
    edt = ndi.distance_transform_edt(mask, sampling=spacing)
    skel = skeletonize(mask)
    deg = ndi.convolve(skel.astype(np.uint8), NEIGH.astype(np.uint8),
                       mode="constant", cval=0)[skel]
    S = csr.Skeleton(skel.astype(np.float64) * edt, spacing=spacing)
    summ = csr.summarize(S, separator="-")
    blen_vox = summ["branch-distance"].to_numpy() / min(spacing)

    h = np.bincount(deg, minlength=15)[:15] / len(deg)
    ax = axes[0, j]
    ax.bar(np.arange(15), h * 100, color=col)
    ax.axvspan(1.5, 2.5, color="green", alpha=0.13)
    ax.set_title(f"{key}\nskeleton-voxel degree (26-conn)", fontsize=10)
    ax.set_xlabel("degree")
    ax.set_ylabel("% of skeleton voxels")
    ax.set_ylim(0, 100)
    ax.text(0.97, 0.93,
            f"deg2 = {h[2]*100:.1f}%\ndeg>=4 = {h[4:].sum()*100:.1f}%\n"
            f"mean = {deg.mean():.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(fc="white", ec="0.7"))
    ax.grid(alpha=0.25, axis="y")

    ax = axes[1, j]
    bins = np.arange(0, 41)
    ax.hist(np.clip(blen_vox, 0, 40), bins=bins, color=col)
    ax.axvline(np.median(blen_vox), color="k", ls="--", lw=1.4,
               label=f"median = {np.median(blen_vox):.1f} vox")
    ax.axvline(5, color="crimson", ls=":", lw=1.6, label="5 voxels")
    ax.set_xlabel("branch length (voxels)")
    ax.set_ylabel("count")
    ax.set_title(f"n_branches = {len(blen_vox):,}; "
                 f"{100*(blen_vox<5).mean():.0f}% shorter than 5 vox",
                 fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25, axis="y")
    del ni, mask, edt, skel

fig.suptitle("PHASE 3 GATE FAILURE — skeletonize() returns a CURVE skeleton for "
             "the fine anode (deg2 = 83%) but a MEDIAL SHEET for medium and "
             "coarse (deg>=4 = 79% and 94%).\nA curve skeleton has mean degree "
             "~2; medium/coarse have ~6.7 and ~7.7, the signature of a 2D "
             "surface. Branch statistics are therefore not comparable between "
             "anodes.", fontsize=10.5)
fig.tight_layout(rect=[0, 0, 1, 0.93])
dest = os.path.join(OUT, "phase3_GATE_FAILURE_skeleton_dimensionality.png")
fig.savefig(dest, dpi=145)
print("[saved]", dest)
