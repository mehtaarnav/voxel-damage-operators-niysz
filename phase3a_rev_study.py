"""
PHASE 3a — REV (representative elementary volume) study, to CHOOSE the ROI size
from data instead of guessing it.

WHY THIS IS NEEDED
------------------
Memory caps us at roughly 100-150 Mvoxel per analysed sub-volume (the Euclidean
distance transform alone is 8 bytes/voxel).  The published analyses used image
windows of 5660-12246 um^3.  The coarse anode is made from a YSZ powder with
d50 = 10.19 um (supplementary Table S1 of ma8095265), so a small cube may
contain less than one starting particle and would not be representative.

Rather than pick an ROI size by feel, we measure the convergence of two
quantities against sub-volume size L:
    * Ni volume fraction               (the classic REV criterion)
    * Ni percolating (spanning) fraction, 6-connectivity, z-direction
      -- computed with the Phase-0-validated code
at several independent positions per size, so we get both a mean and a spread.

CONVENTIONS
-----------
  * cubes are CUBES IN PHYSICAL SPACE (um), not in voxels, so that samples with
    different voxel sizes are compared like for like.
  * 6-connectivity, face-to-face spanning, free boundaries (as validated in
    Phase 0).
  * positions: a 2x2x2 grid of non-overlapping-as-possible corners within the
    loaded block, plus the centre, giving up to 9 samples per size.

Outputs: out/phase3/phase3a_rev_*.csv/.png
"""

from __future__ import annotations

import os
import sys
import time

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE, ground_truth_frame  # noqa: E402
from cmlib.io import label_histogram, load_subvolume  # noqa: E402
from cmlib.percolation import percolation_report  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase3")
os.makedirs(OUT, exist_ok=True)

# physical sub-volume side lengths to test, um
L_LIST = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]
BLOCK_UM = 13.0          # side of the block we load once per sample
PRISTINE = ["fine_pre", "medium_pre", "coarse_pre"]


def sample_meta(key):
    return [s for s in SAMPLES if s[0] == key][0]


def main():
    only_pristine = "--all" not in sys.argv
    keys = PRISTINE if only_pristine else [s[0] for s in SAMPLES]

    print("=" * 78)
    print("PHASE 3a — REV study (choosing the ROI size from data)")
    print("=" * 78)

    rows = []
    for key in keys:
        _, folder, grain, state, nx, ny, nz, vx, vy, vz = sample_meta(key)
        vx_um, vy_um, vz_um = vx / 1000, vy / 1000, vz / 1000
        ext = (nx * vx_um, ny * vy_um, nz * vz_um)
        print(f"\n--- {key} ({folder}) ---")
        print(f"  voxel {vx:.2f} x {vy:.2f} x {vz:.2f} nm   "
              f"extent {ext[0]:.2f} x {ext[1]:.2f} x {ext[2]:.2f} um")

        # block size in voxels, clipped to the stack
        bx = min(int(round(BLOCK_UM / vx_um)), nx)
        by = min(int(round(BLOCK_UM / vy_um)), ny)
        bz = min(int(round(BLOCK_UM / vz_um)), nz)
        x0 = (nx - bx) // 2
        y0 = (ny - by) // 2
        z0 = (nz - bz) // 2
        print(f"  loading block {bz} x {by} x {bx} voxels "
              f"= {bz*by*bx/1e6:.1f} Mvoxel "
              f"({bz*vz_um:.2f} x {by*vy_um:.2f} x {bx*vx_um:.2f} um)")
        t = time.time()
        vol = load_subvolume(folder, z0, z0 + bz, y0, y0 + by, x0, x0 + bx)
        print(f"  loaded in {time.time()-t:.1f} s, {vol.nbytes/1e6:.0f} MB")

        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        ni = vol == mapping["Ni"]
        del vol
        print(f"  Ni fraction in block: {ni.mean():.4f}")

        for L in L_LIST:
            sx = int(round(L / vx_um))
            sy = int(round(L / vy_um))
            sz = int(round(L / vz_um))
            if sx > bx or sy > by or sz > bz:
                continue
            # positions: corners of a 2x2x2 grid + centre
            xs = sorted(set([0, bx - sx, (bx - sx) // 2]))
            ys = sorted(set([0, by - sy, (by - sy) // 2]))
            zs = sorted(set([0, bz - sz, (bz - sz) // 2]))
            pos = []
            for a in xs[:2] + xs[2:]:
                for b in ys[:2] + ys[2:]:
                    for c in zs[:2] + zs[2:]:
                        pos.append((c, b, a))
            pos = list(dict.fromkeys(pos))[:9]

            phis, pfracs, percs = [], [], []
            for (c, b, a) in pos:
                sub = ni[c:c + sz, b:b + sy, a:a + sx]
                rep = percolation_report(sub, axis=0, connectivity=6)
                phis.append(rep["volume_fraction"])
                pfracs.append(rep["percolating_frac"])
                percs.append(rep["percolates"])
            rows.append(dict(
                sample=key, grain=grain, state=state, L_um=L,
                n_pos=len(pos),
                nvox=sx * sy * sz,
                phi_mean=float(np.mean(phis)), phi_std=float(np.std(phis)),
                pfrac_mean=float(np.mean(pfracs)),
                pfrac_std=float(np.std(pfracs)),
                perc_rate=float(np.mean(percs)),
            ))
            print(f"    L={L:5.1f} um  ({sz}x{sy}x{sx} vox, n={len(pos)})  "
                  f"phi_Ni={np.mean(phis):.4f}+-{np.std(phis):.4f}   "
                  f"perc.frac={np.mean(pfracs):.4f}+-{np.std(pfracs):.4f}   "
                  f"spanning {np.mean(percs)*100:.0f}%")
        del ni

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase3a_rev.csv"), index=False)

    gt = ground_truth_frame().set_index("sample")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    colours = {"fine_pre": "C0", "medium_pre": "C1", "coarse_pre": "C2",
               "fine_post": "C3", "medium_post": "C4", "coarse_post": "C5"}
    for key in df["sample"].unique():
        d = df[df["sample"] == key]
        c = colours.get(key, "k")
        axes[0].errorbar(d.L_um, d.phi_mean, yerr=d.phi_std, marker="o", ms=4,
                         lw=1.2, capsize=3, color=c, label=key)
        axes[0].axhline(gt.loc[key, "Ni_Phi__T-S4"], color=c, ls=":", lw=1.0)
        axes[1].errorbar(d.L_um, d.pfrac_mean, yerr=d.pfrac_std, marker="s",
                         ms=4, lw=1.2, capsize=3, color=c, label=key)
        axes[2].plot(d.L_um, d.phi_std / d.phi_mean, "^-", ms=4, lw=1.2,
                     color=c, label=key)
    axes[0].set_xlabel("sub-volume side L (um)")
    axes[0].set_ylabel(r"Ni volume fraction $\Phi$")
    axes[0].set_title("REV: $\\Phi_{Ni}$ vs L\n(dotted = published Table S4)")
    axes[1].set_xlabel("sub-volume side L (um)")
    axes[1].set_ylabel("Ni percolating fraction (z, 6-conn)")
    axes[1].set_title("REV: spanning fraction vs L")
    axes[2].set_xlabel("sub-volume side L (um)")
    axes[2].set_ylabel(r"relative scatter  $\sigma(\Phi)/\langle\Phi\rangle$")
    axes[2].set_title("REV: scatter between positions")
    axes[2].axhline(0.05, color="crimson", ls="--", lw=1.2,
                    label="5 % scatter")
    for a in axes:
        a.grid(alpha=0.25)
        a.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase3a_rev.png"), dpi=150)
    plt.close(fig)

    print("\n[saved]", os.path.join(OUT, "phase3a_rev.csv"))
    print("[saved]", os.path.join(OUT, "phase3a_rev.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
