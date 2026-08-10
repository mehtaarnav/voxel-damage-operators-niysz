"""
PHASE 3 — extract the Ni-phase graph from the pristine stacks.

Per the agreed hybrid strategy: graph metrics are computed on several
NON-NESTED ROIs of a common physical size per anode, at native resolution, so
that the between-ROI spread becomes a genuine error bar and the coarse anode's
sub-REV status is visible rather than hidden.

Steps (see cmlib/graph.py for every convention and its justification):
  1. extract Ni mask for the ROI
  2. keep only the percolating (x-spanning) cluster
  3. distance transform with true physical voxel spacing
  4. 3D skeletonize
  5. skan Skeleton -> summarize -> skeleton_to_nx
  6. edge weight = neck width w = 2*min(EDT) along the branch; conductance
     (w/2)^2 / L
  7. save a skeleton-on-Ni overlay per anode for visual inspection (GATE)

Usage:
    python phase3_extract_graph.py                 # all pristine ROIs
    python phase3_extract_graph.py --limit 1       # one ROI per anode (smoke test)
    python phase3_extract_graph.py --roi-um 8.0
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
import time

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.graph import build_ni_graph, largest_component, weight_edges  # noqa: E402
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, load_subvolume  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase3")
GRAPHS = os.path.join(HERE, "out", "graphs")
os.makedirs(OUT, exist_ok=True)
os.makedirs(GRAPHS, exist_ok=True)

TRANSPORT_AXIS = 2      # x, through-thickness (see cmlib/graph.py convention 2)
PRISTINE = ["fine_pre", "medium_pre", "coarse_pre"]


def save_overlay(key, roi_name, mask, skel, edt, spacing_nm, dest):
    """Skeleton drawn on the Ni phase, for the visual gate."""
    mid = mask.shape[0] // 2
    m2 = mask[mid]
    s2 = skel[mid]
    e2 = edt[mid]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.4))

    axes[0].imshow(m2, cmap="gray", interpolation="nearest")
    ys, xs = np.nonzero(s2)
    axes[0].plot(xs, ys, ".", ms=1.6, color="#ff2d55")
    axes[0].set_title(f"Ni (white) + skeleton (red)\nslice z={mid}", fontsize=9)

    # thick-slab projection: skeleton within +/-4 voxels, so 3D continuity shows
    lo, hi = max(0, mid - 4), min(mask.shape[0], mid + 5)
    slab = skel[lo:hi].any(axis=0)
    axes[1].imshow(mask[lo:hi].any(axis=0), cmap="gray", interpolation="nearest")
    ys, xs = np.nonzero(slab)
    axes[1].plot(xs, ys, ".", ms=1.6, color="#00e5ff")
    axes[1].set_title(f"9-slice slab projection\n(skeleton cyan)", fontsize=9)

    im = axes[2].imshow(np.where(m2, e2, np.nan), cmap="magma",
                        interpolation="nearest")
    axes[2].set_title("distance transform inside Ni (nm)\n= local inscribed radius",
                      fontsize=9)
    fig.colorbar(im, ax=axes[2], fraction=0.04)

    for a in axes:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"{key}  ROI {roi_name}   voxel "
                 f"{spacing_nm[0]:.2f} x {spacing_nm[1]:.2f} x {spacing_nm[2]:.2f} nm",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def save_neck_hist(key, necks_by_roi, dest):
    fig, ax = plt.subplots(figsize=(7.5, 5))
    allv = np.concatenate([v for v in necks_by_roi.values() if len(v)])
    bins = np.linspace(0, np.percentile(allv, 99.5), 60)
    for name, v in necks_by_roi.items():
        if len(v):
            ax.hist(v, bins=bins, histtype="step", lw=1.1, density=True,
                    label=f"{name} (n={len(v)})")
    ax.axvline(np.percentile(allv, 10), color="crimson", ls="--", lw=1.6,
               label=f"pooled 10th pct = {np.percentile(allv,10):.0f} nm")
    ax.set_xlabel("branch neck width  w = 2·min(EDT)  [nm]")
    ax.set_ylabel("density")
    ax.set_title(f"{key}: distribution of branch neck widths, per ROI")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roi-um", type=float, default=8.0)
    ap.add_argument("--limit", type=int, default=None,
                    help="max ROIs per sample")
    ap.add_argument("--samples", nargs="*", default=PRISTINE)
    args = ap.parse_args()

    print("=" * 78)
    print(f"PHASE 3 — Ni graph extraction   ROI = {args.roi_um} um cubes")
    print("=" * 78)
    print(f"transport axis = {TRANSPORT_AXIS} (x, through-thickness)")
    print("graph built on the x-spanning Ni cluster only; 6-connectivity\n")

    rows = []
    for key in args.samples:
        _, folder, grain, state, nx_, ny_, nz_, vx, vy, vz = \
            [s for s in SAMPLES if s[0] == key][0]
        spacing = (vz, vy, vx)          # array axes are (z, y, x)
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])

        rois = tile_rois(nz_, ny_, nx_, vz, vy, vx, args.roi_um, args.limit)
        print(f"\n=== {key} ({folder}) ===")
        print(f"  voxel {vx:.2f} x {vy:.2f} x {vz:.2f} nm; Ni label {mapping['Ni']}")
        if not rois:
            print(f"  !! no ROI of {args.roi_um} um fits — SKIPPED")
            continue
        print(f"  {len(rois)} ROIs of {rois[0]['nz']}x{rois[0]['ny']}x{rois[0]['nx']}"
              f" voxels = {rois[0]['nvox']/1e6:.1f} Mvoxel "
              f"({rois[0]['size_um'][0]:.2f} x {rois[0]['size_um'][1]:.2f} x "
              f"{rois[0]['size_um'][2]:.2f} um)")

        necks_by_roi = {}
        for i, r in enumerate(rois):
            t0 = time.time()
            vol = load_subvolume(folder, r["z0"], r["z1"], r["y0"], r["y1"],
                                 r["x0"], r["x1"])
            ni = vol == mapping["Ni"]
            del vol
            t_load = time.time() - t0

            t0 = time.time()
            G, diag, arrs = build_ni_graph(ni, spacing, axis=TRANSPORT_AXIS)
            t_graph = time.time() - t0

            row = dict(sample=key, grain=grain, state=state, roi=r["roi"],
                       roi_um=args.roi_um, nvox=r["nvox"],
                       t_load_s=round(t_load, 1), t_graph_s=round(t_graph, 1))
            row.update({k: v for k, v in diag.items()
                        if not isinstance(v, (list, dict))})

            if G is not None and G.number_of_edges() > 0:
                G = weight_edges(G)
                Gc = largest_component(G)
                necks = np.array([d["neck_nm"] for _, _, d in Gc.edges(data=True)])
                necks_by_roi[r["roi"]] = necks
                row.update(
                    lcc_nodes=Gc.number_of_nodes(),
                    lcc_edges=Gc.number_of_edges(),
                    neck_p10_nm=float(np.percentile(necks, 10)),
                    neck_p50_nm=float(np.percentile(necks, 50)),
                    neck_mean_nm=float(necks.mean()),
                )
                with open(os.path.join(
                        GRAPHS, f"{key}__{r['roi']}__{args.roi_um}um.pkl"), "wb") as f:
                    pickle.dump(Gc, f)
                print(f"  ROI {r['roi']}: Ni={diag['perc_volume_fraction']:.4f} "
                      f"span={diag['perc_percolates']} "
                      f"pfrac={diag['perc_percolating_frac']:.4f}  "
                      f"skel={diag.get('skeleton_voxels',0):,} vox  "
                      f"G=({Gc.number_of_nodes():,}n,{Gc.number_of_edges():,}e)  "
                      f"neck_p10={row['neck_p10_nm']:.0f}nm  "
                      f"[{t_load:.0f}s load, {t_graph:.0f}s graph]")

                if i == 0 and arrs:
                    dest = os.path.join(OUT, f"phase3_overlay_{key}.png")
                    save_overlay(key, r["roi"], arrs["mask"], arrs["skel"],
                                 arrs["edt"], spacing, dest)
                    print(f"    [overlay] {os.path.basename(dest)}")
            else:
                print(f"  ROI {r['roi']}: NO GRAPH — {diag.get('reason')}")
            rows.append(row)
            del ni, arrs

        if necks_by_roi:
            dest = os.path.join(OUT, f"phase3_neckhist_{key}.png")
            save_neck_hist(key, necks_by_roi, dest)
            print(f"  [hist] {os.path.basename(dest)}")

    df = pd.DataFrame(rows)
    dest = os.path.join(OUT, f"phase3_graphs_{args.roi_um}um.csv")
    df.to_csv(dest, index=False)
    print(f"\n[saved] {dest}")
    if not df.empty:
        cols = [c for c in ("sample", "roi", "perc_volume_fraction",
                            "perc_percolating_frac", "skeleton_voxels",
                            "lcc_nodes", "lcc_edges", "neck_p10_nm",
                            "neck_p50_nm", "t_graph_s") if c in df.columns]
        print(df[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
