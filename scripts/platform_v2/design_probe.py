"""
PLATFORM V2 — Ni generator qualification, design probe.

Answers the open question posed before any topology-strategy decision: does
PLAIN simple-cubic (6-connectivity, nearest-neighbour-only) lattice adjacency
already land SNOW-measured mean_degree inside the newly-grounded 3.5-4.5 band
at the new Phi_Ni~0.250 / 160-192-vox domain target, or does it require
topology modification (face-diagonal bonds etc.)?

Also solves the geometry (R, pitch, nlat_z, nlat_xy, margin) numerically
(measured via direct rasterization, not analytic approximation) to hit
Phi_Ni~0.250 within the 160-192-vox-per-axis domain target and 24-32-voxel
particle diameter, per preregistration.md #0e / the platform-v2 review.

Does NOT build the qualification run (5+ seeds, gating table) -- that is
scripts/platform_v2/qualification_run.py, run only after this probe settles
the topology question.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from cmlib.api import compute_network_metrics, compute_percolation, extract_network  # noqa: E402
from cmlib.percolation import label_phase  # noqa: E402
from cmlib.synth import (  # noqa: E402
    base_distribution_stats, mixture_neck_widths, platform_v2_lattice_geometry,
    rasterize,
)

OUT = os.path.join(ROOT, "out", "platform_v2")
os.makedirs(OUT, exist_ok=True)

VOXEL_NM = 20.0
R_VOX = 12.1          # tuned by direct search: R=12.0->-0.38%, 12.1->+0.10%,
                      # 12.2->+0.59% deviation from Phi_Ni=0.250 (with base
                      # necks included) -- 12.1 is the closest tested
PITCH_VOX = 32
NLAT_Z = 6
NLAT_XY = 4
MARGIN = 8
JITTER_FRAC = 0.15

# candidate base neck-width mixture, scaled down from the Family B pilot's
# (R=14) design in proportion to the new R=12.5 (scale factor 12.5/14=0.893),
# floor kept at the 4-voxel resolution requirement (cannot shrink further)
FRAC_WEAK = 0.20
WEAK_RANGE = (4, 6)
NORMAL_RANGE = (12, 20)


def main():
    print("=" * 78)
    print("PLATFORM V2 DESIGN PROBE")
    print("=" * 78)
    t0 = time.time()

    rng = np.random.default_rng(999)
    centres, pairs, shape = platform_v2_lattice_geometry(
        NLAT_Z, NLAT_XY, PITCH_VOX, R_VOX, MARGIN, JITTER_FRAC, rng)
    n_sites = len(centres)
    n_pairs = len(pairs)
    raw_mean_degree = 2.0 * n_pairs / n_sites
    print(f"geometry: R={R_VOX}vox (D={2*R_VOX}vox={2*R_VOX*VOXEL_NM:.0f}nm) "
          f"pitch={PITCH_VOX}vox nlat_z={NLAT_Z} nlat_xy={NLAT_XY} margin={MARGIN}")
    print(f"domain shape = {shape}  "
          f"(target: every axis in [160,192]: "
          f"{all(160 <= s <= 192 for s in shape)})")
    print(f"n_sites={n_sites}  n_pairs={n_pairs}  "
          f"RAW TOPOLOGICAL mean_degree = 2*{n_pairs}/{n_sites} = "
          f"{raw_mean_degree:.3f}")
    print(f"  (this is the pairs-list mean degree BEFORE any watershed/SNOW "
          f"segmentation effect -- the actual measured value below may "
          f"differ due to node merging or splitting)\n")

    # ---- spheres-only Phi_Ni (no necks) -----------------------------------
    spheres_only = rasterize(centres, pairs, R_VOX, np.zeros(n_pairs), shape)
    phi_spheres_only = float(spheres_only.mean())
    print(f"Phi_Ni from spheres alone (no necks): {phi_spheres_only:.4f}  "
          f"(target 0.250)")
    del spheres_only

    # ---- representative base structure (spheres + mixture necks) ---------
    base_rng = np.random.default_rng(0)
    widths = mixture_neck_widths(n_pairs, base_rng, FRAC_WEAK, WEAK_RANGE,
                                 NORMAL_RANGE)
    stats = base_distribution_stats(widths)
    print(f"\nbase neck-width mixture ({int(FRAC_WEAK*100)}% weak"
          f"{WEAK_RANGE}vox + {int((1-FRAC_WEAK)*100)}% normal{NORMAL_RANGE}vox), "
          f"n_pairs={n_pairs}:")
    print(f"  p10={stats['p10_vox']:.2f}vox  p50={stats['p50_vox']:.2f}vox  "
          f"p50/p10={stats['ratio']:.2f}  min={stats['min_vox']:.0f}vox  "
          f"max={stats['max_vox']:.0f}vox")
    print(f"  target p50/p10 in [3.0,4.3], validity floor >=2.5: "
          f"{'PASS' if stats['ratio'] >= 2.5 else 'FAIL'}")
    print(f"  base p10 resolved by >={stats['p10_vox']:.1f} voxels "
          f"(target >=3, preferably >=4): "
          f"{'PASS' if stats['p10_vox'] >= 3 else 'FAIL'}")

    ni_mask = rasterize(centres, pairs, R_VOX, widths, shape)
    phi_ni = float(ni_mask.mean())
    print(f"\nPhi_Ni with base necks included: {phi_ni:.4f}  (target 0.250, "
          f"deviation {100*(phi_ni-0.25)/0.25:+.1f}%)  [{time.time()-t0:.0f}s]")

    # ---- percolation / connectivity check ---------------------------------
    perc = compute_percolation(
        np.where(ni_mask, 2, 0).astype(np.uint8),
        {"Ni": 2, "YSZ": 1, "pore": 0}, phase="Ni", axis=0)
    _, n_clusters = label_phase(ni_mask, connectivity=6)
    print(f"P_span={perc['P_span']:.3f}  P_reach={perc['P_reach']:.3f}  "
          f"n_clusters={n_clusters}  "
          f"single-connected-cluster: {'YES' if n_clusters == 1 else 'NO'}")

    # ---- THE OPEN QUESTION: plain-lattice SNOW mean_degree -----------------
    print("\n" + "=" * 78)
    print("OPEN QUESTION: does PLAIN lattice adjacency land in [3.5, 4.5] "
          "without any topology modification?")
    print("=" * 78)
    G, diag = extract_network(ni_mask, (VOXEL_NM,) * 3, axis=0, r_max=4)
    if G is None or G.number_of_edges() == 0:
        print("  NO NETWORK EXTRACTED -- cannot answer; geometry needs revision")
        return 1
    nm = compute_network_metrics(G, G.graph.get("face_lo"), G.graph.get("face_hi"))
    degrees = np.array([d for _, d in G.degree()])
    n_deg0 = int((degrees == 0).sum())    # isolated (shouldn't exist in G by construction)
    n_deg1 = int((degrees == 1).sum())
    print(f"  n_nodes={nm['n_nodes']}  n_edges={nm['n_edges']}  "
          f"SNOW-measured mean_degree={nm['mean_degree']:.3f}  "
          f"[{time.time()-t0:.0f}s]")
    print(f"  degree distribution (min={degrees.min()}, max={degrees.max()}, "
          f"median={np.median(degrees):.0f}):")
    for d in range(int(degrees.max()) + 1):
        c = int((degrees == d).sum())
        if c:
            print(f"      degree {d}: {c:4d} nodes ({100*c/len(degrees):5.1f}%)")
    print(f"  degree-0 (isolated) chambers: {n_deg0}")
    print(f"  degree-1 (dead-end/pendant) chambers: {n_deg1}")

    in_band = 3.5 <= nm["mean_degree"] <= 4.5
    print(f"\n  VERDICT: measured mean_degree {nm['mean_degree']:.3f} is "
          f"{'INSIDE' if in_band else 'OUTSIDE'} the target band [3.5, 4.5].")
    if in_band:
        print("  -> PLAIN lattice adjacency (6-connectivity, nearest-neighbour "
              "only) ALREADY satisfies the coordination target at this "
              "Phi_Ni/domain size. NO topology modification is needed.")
        print("  -> Per instruction: report this plainly, do not build "
              "face-diagonal / proximity-edge topology machinery.")
    else:
        side = "above" if nm["mean_degree"] > 4.5 else "below"
        print(f"  -> Plain lattice adjacency falls {side} the target band. "
              f"Topology modification IS required (see instructions for "
              f"next step -- NOT built by this probe).")

    print(f"\n[done] total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
