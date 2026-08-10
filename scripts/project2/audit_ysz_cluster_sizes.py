"""
AUDIT: is the real YSZ cluster count speckle-dominated or fragment-dominated?

The Trap-2 argument holds that real `n_clusters` is dominated by FIB-SEM
segmentation speckle, because there is not enough isolated volume to support
that many WHOLE particles. The arithmetic is right, but it constrains only the
MEAN size; a heavy-tailed distribution could put the volume in a few large
fragments and the count in many single-voxel specks -- or the fragments could
all be genuinely large. Those two possibilities imply opposite fixes, so this
measures the distribution instead of assuming it.

Reports, for the YSZ phase of each pristine stack: the full cluster-size
distribution (excluding the largest/spanning cluster), the fraction of CLUSTERS
that are voxel-scale speckle, and the fraction of isolated VOLUME they carry.
Also computes the size of a whole particle at that anode's real diameter, so
"how many whole particles could the isolated volume support" is answered
directly.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, slice_paths  # noqa: E402
from cmlib.percolation import structure_for  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402

OUT = os.path.join(ROOT, "out", "project2")
SCRATCH = os.path.join(OUT, "_tmp")
CONN = 6
VOX_NM = 20.0          # nominal; per-sample spacing used for volumes below

# real mean Ni/YSZ particle diameters (nm) -- YSZ grain size is not separately
# measured in this study, so the Ni particle diameter is used as the scale
# proxy and this is stated as an assumption, not a measurement.
REAL_D_NM = {"fine": 1148.0, "medium": 1445.0, "coarse": 1715.0}
SPECKLE_MAX_VOX = 8    # "voxel-scale": <= 8 voxels (a 2x2x2 block)


def load_mask(folder, val):
    ps = slice_paths(folder)
    a0 = tifffile.imread(ps[0])
    out = np.empty((len(ps), a0.shape[0], a0.shape[1]), dtype=bool)
    out[0] = a0 == val
    for i, p in enumerate(ps[1:], start=1):
        out[i] = tifffile.imread(p) == val
    return out


def cluster_sizes(mask, slab=32):
    os.makedirs(SCRATCH, exist_ok=True)
    import tempfile
    fh, path = tempfile.mkstemp(suffix=".i4", dir=SCRATCH)
    os.close(fh)
    try:
        lab = np.memmap(path, dtype=np.int32, mode="w+", shape=mask.shape)
        n = int(ndi.label(mask, structure=structure_for(CONN), output=lab))
        counts = np.zeros(n + 1, dtype=np.int64)
        for z0 in range(0, mask.shape[0], slab):
            counts += np.bincount(np.asarray(lab[z0:z0 + slab]).ravel(),
                                  minlength=n + 1)
        counts[0] = 0
        del lab
    finally:
        try:
            os.remove(path)
        except OSError:
            pass
    return counts[1:]


def main():
    rows = []
    for key, folder, grain, state, nx_, ny_, nz_, vx, vy, vz in SAMPLES:
        if state != "pristine":
            continue
        counts_hist = label_histogram(folder)["counts"]
        mapping = assign_labels(counts_hist, ZENODO_LABEL_NOTE[key])
        mask = load_mask(folder, mapping["YSZ"])
        dom = mask.size
        sizes = cluster_sizes(mask)
        del mask

        vox_nm3 = float(vx) * float(vy) * float(vz)
        sizes = np.sort(sizes)[::-1]
        largest, rest = sizes[0], sizes[1:]
        tot = sizes.sum()
        iso_vox = int(rest.sum())

        d = REAL_D_NM[grain]
        particle_vox = (np.pi / 6.0) * d ** 3 / vox_nm3

        speckle = rest[rest <= SPECKLE_MAX_VOX]
        big = rest[rest > particle_vox]
        halfp = rest[rest > 0.1 * particle_vox]

        r = dict(
            sample=key, grain=grain,
            domain_Mvox=dom / 1e6, ysz_vox=int(tot),
            n_clusters=int(sizes.size),
            n_clusters_per_Mvox=sizes.size / (dom / 1e6),
            P_largest=float(largest / tot),
            iso_frac_of_phase=float(iso_vox / tot),
            iso_vox=iso_vox,
            particle_vox=float(particle_vox),
            iso_in_whole_particles=float(iso_vox / particle_vox),
            whole_particles_per_Mvox=float(
                (iso_vox / particle_vox) / (dom / 1e6)),
            mean_iso_cluster_vox=float(rest.mean()) if rest.size else 0.0,
            median_iso_cluster_vox=float(np.median(rest)) if rest.size else 0.0,
            p90_iso_cluster_vox=float(np.percentile(rest, 90)) if rest.size else 0.0,
            max_iso_cluster_vox=int(rest.max()) if rest.size else 0,
            n_speckle_le8vox=int(speckle.size),
            frac_clusters_speckle=float(speckle.size / rest.size) if rest.size else 0.0,
            frac_isovol_in_speckle=float(speckle.sum() / iso_vox) if iso_vox else 0.0,
            n_gt_particle=int(big.size),
            frac_isovol_gt_particle=float(big.sum() / iso_vox) if iso_vox else 0.0,
            n_gt_tenth_particle=int(halfp.size),
            n_gt_tenth_particle_per_Mvox=float(halfp.size / (dom / 1e6)),
            frac_isovol_gt_tenth=float(halfp.sum() / iso_vox) if iso_vox else 0.0,
        )
        rows.append(r)
        print(f"\n=== {key} ({grain}) ===")
        print(f"  domain {r['domain_Mvox']:.0f} Mvox   YSZ {tot/1e6:.1f} Mvox   "
              f"clusters {r['n_clusters']:,} ({r['n_clusters_per_Mvox']:.2f}/Mvox)")
        print(f"  P_largest={r['P_largest']:.4f}  isolated={r['iso_frac_of_phase']*100:.2f}% "
              f"of phase = {iso_vox:,} vox")
        print(f"  one whole particle = {particle_vox:,.0f} vox "
              f"(d={d:.0f} nm); isolated volume = "
              f"{r['iso_in_whole_particles']:.1f} whole particles "
              f"({r['whole_particles_per_Mvox']:.3f}/Mvox)")
        print(f"  isolated cluster size: median={r['median_iso_cluster_vox']:.0f} "
              f"mean={r['mean_iso_cluster_vox']:.0f} p90={r['p90_iso_cluster_vox']:.0f} "
              f"max={r['max_iso_cluster_vox']:,} vox")
        print(f"  speckle (<= {SPECKLE_MAX_VOX} vox): "
              f"{r['frac_clusters_speckle']*100:.1f}% of CLUSTERS but "
              f"{r['frac_isovol_in_speckle']*100:.2f}% of isolated VOLUME")
        print(f"  clusters > 0.1 particle: {r['n_gt_tenth_particle']:,} "
              f"({r['n_gt_tenth_particle_per_Mvox']:.3f}/Mvox), carrying "
              f"{r['frac_isovol_gt_tenth']*100:.1f}% of isolated volume")
        print(f"  clusters > 1 particle:   {r['n_gt_particle']:,}, carrying "
              f"{r['frac_isovol_gt_particle']*100:.1f}%")
        pd.DataFrame(rows).to_csv(
            os.path.join(OUT, "audit_ysz_cluster_sizes.csv"), index=False)

    try:
        os.rmdir(SCRATCH)
    except OSError:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
