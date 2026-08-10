"""
DIAGNOSTIC — is the 3D skeleton a clean CURVE skeleton, or a sheet/hairy mess?

Symptom that triggered this: at 8 um ROIs the mean branch length is 11 voxels
for fine_pre but only 3.7 voxels for coarse_pre, and coarse_pre yields 4x MORE
graph edges than fine_pre from a 3x SMALLER volume.  That ordering is backwards
for a coarser microstructure and suggests the skeleton quality differs
systematically between anodes -- which would bias every graph metric.

Tests:
  A. Reference shapes with known skeletons (cylinder -> line, sphere -> point,
     slab -> sheet), to confirm skeletonize()'s 3D behaviour in this version.
  B. On real ROIs: the distribution of skeleton-voxel DEGREE (number of
     26-connected skeleton neighbours).  A curve skeleton is dominated by
     degree 2.  A SHEET skeleton has a large fraction with degree >= 4.
  C. Branch length and branch-type distributions, and how much of the graph is
     terminal spurs (branch-type 1) vs real junction-to-junction links (type 2).
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skan import csr
from skimage.morphology import skeletonize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, load_subvolume  # noqa: E402
from cmlib.percolation import percolating_mask  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402

NEIGH = ndi.generate_binary_structure(3, 3)   # 26-connectivity
NEIGH[1, 1, 1] = False


def degree_profile(skel):
    """Fraction of skeleton voxels with each number of 26-conn skeleton neighbours."""
    deg = ndi.convolve(skel.astype(np.uint8), NEIGH.astype(np.uint8),
                       mode="constant", cval=0)
    d = deg[skel]
    hist = np.bincount(d, minlength=15)[:15]
    return hist / max(hist.sum(), 1), d


def report_shape(name, vol):
    skel = skeletonize(vol)
    prof, d = degree_profile(skel)
    print(f"  {name:22s} vox={vol.sum():8d}  skel={skel.sum():6d}  "
          f"deg2={prof[2]*100:5.1f}%  deg>=4={prof[4:].sum()*100:5.1f}%  "
          f"meandeg={d.mean():.2f}")


def main():
    print("=" * 78)
    print("A. Reference shapes (does skeletonize() behave as expected in 3D?)")
    print("=" * 78)
    N = 80
    zz, yy, xx = np.ogrid[:N, :N, :N]

    cyl = ((yy - 40) ** 2 + (xx - 40) ** 2) < 8 ** 2
    cyl = np.broadcast_to(cyl, (N, N, N)).copy()
    cyl[:5] = False; cyl[-5:] = False
    report_shape("long cylinder", cyl)

    sph = ((zz - 40) ** 2 + (yy - 40) ** 2 + (xx - 40) ** 2) < 20 ** 2
    report_shape("solid sphere", sph)

    slab = np.zeros((N, N, N), bool)
    slab[30:50, 10:70, 10:70] = True
    report_shape("thick slab (plate)", slab)

    blob = ndi.gaussian_filter(
        np.random.default_rng(0).random((N, N, N)).astype(np.float32), 6) > 0.5
    report_shape("smooth random blob", blob)

    print("\n  EXPECTATION: cylinder -> almost all degree 2 (a curve).")
    print("               slab     -> large degree>=4 fraction (a SHEET).")

    print("\n" + "=" * 78)
    print("B/C. Real ROIs (one 8 um ROI per pristine anode)")
    print("=" * 78)
    rows = []
    for key in ("fine_pre", "medium_pre", "coarse_pre"):
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
        prof, d = degree_profile(skel)

        S = csr.Skeleton(skel.astype(np.float64) * edt, spacing=spacing)
        summ = csr.summarize(S, separator="-")
        blen = summ["branch-distance"].to_numpy()
        btype = summ["branch-type"].to_numpy()
        vox_len = blen / min(spacing)

        print(f"\n--- {key} (voxel {vx:.2f} nm) ---")
        print(f"  mask voxels          {mask.sum():,}")
        print(f"  skeleton voxels      {skel.sum():,}  "
              f"({100*skel.sum()/max(mask.sum(),1):.3f} % of mask)")
        print(f"  skeleton-voxel degree: deg1={prof[1]*100:5.1f}%  "
              f"deg2={prof[2]*100:5.1f}%  deg3={prof[3]*100:5.1f}%  "
              f"deg>=4={prof[4:].sum()*100:5.1f}%   mean={d.mean():.2f}")
        print(f"  n_branches           {len(blen):,}")
        print(f"  branch length (voxels): median={np.median(vox_len):6.1f}  "
              f"p10={np.percentile(vox_len,10):5.1f}  "
              f"p90={np.percentile(vox_len,90):6.1f}")
        for t, nm in ((0, "endpt-endpt"), (1, "junction-endpt (SPUR)"),
                      (2, "junction-junction"), (3, "cycle")):
            n = int((btype == t).sum())
            print(f"    type {t} {nm:24s}: {n:7,d}  ({100*n/len(btype):5.1f} %)")
        short = vox_len < 5
        print(f"  branches shorter than 5 voxels: {short.sum():,} "
              f"({100*short.mean():.1f} %)")

        rows.append(dict(sample=key, mask_vox=int(mask.sum()),
                         skel_vox=int(skel.sum()),
                         skel_pct=100*skel.sum()/max(mask.sum(), 1),
                         deg2_pct=prof[2]*100, deg_ge4_pct=prof[4:].sum()*100,
                         mean_deg=float(d.mean()), n_branches=len(blen),
                         med_len_vox=float(np.median(vox_len)),
                         spur_pct=100*float((btype == 1).mean()),
                         short_pct=100*float(short.mean())))
        del ni, mask, edt, skel

    df = pd.DataFrame(rows)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "out", "phase3", "diag_skeleton.csv")
    df.to_csv(out, index=False)
    print("\n" + df.to_string(index=False))
    print("\n[saved]", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
