"""
PHASE 4d — mean Ni particle size by watershed, using the explicitly requested
libraries: scipy.ndimage distance transform + skimage.feature.peak_local_max
for markers + skimage.segmentation.watershed.

CONVENTIONS (all stated; watershed marker parameters are exactly the kind of
silent choice that would invalidate a comparison)
------------------------------------------------------------------------------
  * phase          : ALL Ni voxels (not only the percolating cluster) -- this is
                     a conventional descriptor of the microstructure, and the
                     published particle sizes are not percolation-restricted.
  * distance map   : scipy.ndimage.distance_transform_edt with `sampling` set to
                     the true physical voxel size, so the map is in nm and
                     anisotropy is handled.
  * marker blur    : Gaussian, sigma = 0.4 VOXELS.  Suppresses single-voxel
                     noise maxima that would over-segment.  Same value as SNOW.
  * markers        : skimage.feature.peak_local_max on the blurred distance map,
                     restricted to the Ni mask, with min_distance = 4 VOXELS.
                     min_distance is the dominant over/under-segmentation
                     control, so a sensitivity sweep over
                     min_distance in {2, 3, 4, 6, 8} is run and reported rather
                     than assuming 4 is harmless.
  * watershed      : skimage.segmentation.watershed on the NEGATED blurred
                     distance map, with `mask` = Ni, seeded by those markers.
  * particle size  : equivalent-sphere diameter d = (6V/pi)^(1/3) from each
                     region's voxel volume.  Reported as the plain mean, the
                     median, and the VOLUME-WEIGHTED mean.  The volume-weighted
                     mean is the one comparable to a d50 from laser
                     diffraction, because that technique is volume-weighted.
  * boundary regions: regions touching the ROI boundary are truncated and their
                     size is under-estimated.  They are EXCLUDED from the size
                     statistics (but reported as a count), which is the standard
                     unbiased-stereology choice.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import RAW_POWDER_D50, SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, load_subvolume  # noqa: E402
from cmlib.particles import (  # noqa: E402
    SIGMA_VOX_DEFAULT as SIGMA_VOX, size_stats, watershed_particles,
)
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase4")
os.makedirs(OUT, exist_ok=True)

PRISTINE = ["fine_pre", "medium_pre", "coarse_pre"]

# watershed_particles() and size_stats() moved to cmlib/particles.py, unchanged
# in definition, 2026-08-10 -- see out/next/EXECUTION_SPEC.md Phase 0. Imported
# above so this script keeps producing the same phase4d_particles.csv schema.


def save_overlay(key, mask, labels, edt, dest, min_distance):
    mid = mask.shape[0] // 2
    m2, l2 = mask[mid], labels[mid].astype(float)
    l2[l2 == 0] = np.nan
    rng = np.random.default_rng(1)
    nl = int(np.nanmax(labels)) + 1
    perm = rng.permutation(nl)
    l2s = np.where(np.isnan(l2), np.nan,
                   perm[np.nan_to_num(l2, nan=0).astype(int)])

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2))
    axes[0].imshow(m2, cmap="gray", interpolation="nearest")
    axes[0].set_title(f"Ni phase, slice z={mid}", fontsize=9)
    axes[1].imshow(m2, cmap="gray", interpolation="nearest")
    axes[1].imshow(l2s, cmap="tab20", interpolation="nearest", alpha=0.85)
    axes[1].set_title(f"watershed particles\n(min_distance={min_distance} vox, "
                      f"sigma={SIGMA_VOX})", fontsize=9)
    im = axes[2].imshow(np.where(m2, edt[mid], np.nan), cmap="magma",
                        interpolation="nearest")
    axes[2].set_title("distance transform (nm)", fontsize=9)
    fig.colorbar(im, ax=axes[2], fraction=0.04)
    for a in axes:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"{key} — Ni particle segmentation", fontsize=10)
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roi-um", type=float, default=8.0)
    ap.add_argument("--limit", type=int, default=2)
    ap.add_argument("--samples", nargs="*", default=PRISTINE)
    ap.add_argument("--sweep", nargs="*", type=int, default=[2, 3, 4, 6, 8])
    args = ap.parse_args()

    print("=" * 78)
    print("PHASE 4d — Ni particle sizing by watershed")
    print("=" * 78)
    print(f"sigma={SIGMA_VOX} vox; min_distance sweep {args.sweep} vox; "
          f"border regions excluded\n")

    rows = []
    for key in args.samples:
        _, folder, grain, state, nx_, ny_, nz_, vx, vy, vz = \
            [s for s in SAMPLES if s[0] == key][0]
        spacing = (vz, vy, vx)
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        rois = tile_rois(nz_, ny_, nx_, vz, vy, vx, args.roi_um, args.limit)
        print(f"\n=== {key} ({grain}) — {len(rois)} ROI(s) ===")
        for i, r in enumerate(rois):
            vol = load_subvolume(folder, r["z0"], r["z1"], r["y0"], r["y1"],
                                 r["x0"], r["x1"])
            ni = vol == mapping["Ni"]
            del vol
            for md in args.sweep:
                labels, edt, npk = watershed_particles(ni, spacing,
                                                       min_distance=md)
                st, _ = size_stats(labels, spacing)
                if not st:
                    continue
                rows.append(dict(sample=key, grain=grain, roi=r["roi"],
                                 min_distance=md, n_peaks=npk, **st))
                print(f"  ROI {r['roi']} md={md}: peaks={npk:5d} "
                      f"regions_used={st['n_regions_used']:5d}  "
                      f"d_mean={st['d_mean_nm']:7.0f} nm  "
                      f"d_median={st['d_median_nm']:7.0f}  "
                      f"d_volwt={st['d_volweighted_nm']:7.0f}")
                if i == 0 and md == 4:
                    dest = os.path.join(OUT, f"phase4d_watershed_{key}.png")
                    save_overlay(key, ni, labels, edt, dest, md)
                    print(f"    [overlay] {os.path.basename(dest)}")
                del labels, edt
            del ni

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase4d_particles.csv"), index=False)

    print("\n" + "=" * 78)
    print("Sensitivity to min_distance (volume-weighted mean diameter, nm)")
    print("=" * 78)
    piv = df.pivot_table(index="sample", columns="min_distance",
                         values="d_volweighted_nm", aggfunc="mean")
    piv = piv.reindex([s for s in PRISTINE if s in piv.index])
    print(piv.round(0).to_string())
    print("\nOrdering by column (does the fine<medium<coarse ranking survive?):")
    for c in piv.columns:
        order = piv[c].sort_values().index.tolist()
        print(f"  min_distance={c}: {order}")

    print("\nRaw YSZ powder d50 for context (supplementary Table S1): "
          + ", ".join(f"{k}={v} um" for k, v in RAW_POWDER_D50.items()))

    fig, ax = plt.subplots(figsize=(8, 5))
    for s in piv.index:
        ax.plot(piv.columns, piv.loc[s] / 1000.0, "o-", label=s)
    ax.set_xlabel("peak_local_max min_distance (voxels)")
    ax.set_ylabel("volume-weighted mean Ni particle diameter (um)")
    ax.set_title("Phase 4d: sensitivity of particle size to the marker parameter")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase4d_particle_sensitivity.png"), dpi=145)
    plt.close(fig)
    print(f"\n[saved] {os.path.join(OUT, 'phase4d_particles.csv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
