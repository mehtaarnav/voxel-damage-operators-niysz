"""
PHASE 3 (v2) — extract the Ni-phase network by watershed (SNOW) partitioning.

Replaces phase3_extract_graph.py, which failed its visual/statistical gate
(see cmlib/pnm.py header and out/phase3/phase3_GATE_FAILURE_*.png).

Saves, for the visual gate:
    phase3_snow_overlay_<sample>.png   watershed regions on the Ni phase,
                                       plus the distance transform and the
                                       throat-size distribution
Usage:
    python phase3_extract_network.py --limit 1     # smoke test
    python phase3_extract_network.py               # all pristine ROIs
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
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, load_subvolume  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.pnm import extract_ni_network, largest_component  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase3")
NETS = os.path.join(HERE, "out", "networks")
os.makedirs(OUT, exist_ok=True)
os.makedirs(NETS, exist_ok=True)

TRANSPORT_AXIS = 2
PRISTINE = ["fine_pre", "medium_pre", "coarse_pre"]


def save_snow_overlay(key, roi, mask, regions, spacing_nm, G, dest):
    mid = mask.shape[0] // 2
    m2 = mask[mid]
    r2 = regions[mid].astype(float)
    r2[r2 == 0] = np.nan

    edt = ndi.distance_transform_edt(m2, sampling=spacing_nm[1:])

    # randomise region colours so neighbours differ
    rng = np.random.default_rng(0)
    nlab = int(np.nanmax(regions)) + 1
    perm = rng.permutation(nlab)
    r2s = np.where(np.isnan(r2), np.nan, perm[np.nan_to_num(r2, nan=0).astype(int)])

    fig, axes = plt.subplots(1, 4, figsize=(20, 5.2))

    axes[0].imshow(m2, cmap="gray", interpolation="nearest")
    axes[0].set_title(f"percolating Ni cluster\nslice z={mid}", fontsize=9)

    axes[1].imshow(m2, cmap="gray", interpolation="nearest")
    axes[1].imshow(r2s, cmap="tab20", interpolation="nearest", alpha=0.85)
    axes[1].set_title("watershed regions (SNOW)\n= graph nodes", fontsize=9)

    im = axes[2].imshow(np.where(m2, edt, np.nan), cmap="magma",
                        interpolation="nearest")
    axes[2].set_title("distance transform (nm)", fontsize=9)
    fig.colorbar(im, ax=axes[2], fraction=0.04)

    necks = np.array([d["neck_nm"] for _, _, d in G.edges(data=True)])
    axes[3].hist(necks, bins=50, color="C0")
    axes[3].axvline(np.percentile(necks, 10), color="crimson", ls="--", lw=1.6,
                    label=f"p10 = {np.percentile(necks,10):.0f} nm")
    axes[3].axvline(np.median(necks), color="k", ls=":", lw=1.4,
                    label=f"median = {np.median(necks):.0f} nm")
    axes[3].set_xlabel("throat inscribed diameter (neck width) [nm]")
    axes[3].set_ylabel("count")
    axes[3].set_title(f"{G.number_of_nodes():,} nodes, "
                      f"{G.number_of_edges():,} throats", fontsize=9)
    axes[3].legend(fontsize=8, frameon=False)
    axes[3].grid(alpha=0.25)

    for a in axes[:3]:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"{key}  ROI {roi}   SNOW watershed network "
                 f"(sigma=0.4, r_max=4)", fontsize=10)
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roi-um", type=float, default=8.0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--samples", nargs="*", default=PRISTINE)
    ap.add_argument("--r-max", type=int, default=4)
    args = ap.parse_args()

    print("=" * 78)
    print(f"PHASE 3 (v2) — SNOW watershed network   ROI = {args.roi_um} um")
    print("=" * 78)
    print(f"transport axis = {TRANSPORT_AXIS} (x); 6-connectivity; "
          f"sigma=0.4, r_max={args.r_max}\n")

    rows = []
    for key in args.samples:
        _, folder, grain, state, nx_, ny_, nz_, vx, vy, vz = \
            [s for s in SAMPLES if s[0] == key][0]
        spacing = (vz, vy, vx)
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        rois = tile_rois(nz_, ny_, nx_, vz, vy, vx, args.roi_um, args.limit)

        print(f"\n=== {key} ({folder}) ===  {len(rois)} ROIs, "
              f"{rois[0]['nvox']/1e6:.1f} Mvoxel each")
        for i, r in enumerate(rois):
            t0 = time.time()
            vol = load_subvolume(folder, r["z0"], r["z1"], r["y0"], r["y1"],
                                 r["x0"], r["x1"])
            ni = vol == mapping["Ni"]
            del vol
            t1 = time.time()
            G, diag, extras = extract_ni_network(ni, spacing,
                                                 axis=TRANSPORT_AXIS,
                                                 r_max=args.r_max)
            t2 = time.time()

            row = dict(sample=key, grain=grain, state=state, roi=r["roi"],
                       roi_um=args.roi_um, nvox=r["nvox"],
                       t_load_s=round(t1 - t0, 1), t_snow_s=round(t2 - t1, 1))
            row.update({k: v for k, v in diag.items()
                        if not isinstance(v, (list, dict))})

            if G is not None and G.number_of_edges() > 0:
                Gc = largest_component(G)
                # keep face terminal sets, restricted to the LCC
                lo = [n for n in extras["face_lo"] if n in Gc]
                hi = [n for n in extras["face_hi"] if n in Gc]
                Gc.graph["face_lo"] = lo
                Gc.graph["face_hi"] = hi
                Gc.graph["sample"] = key
                Gc.graph["roi"] = r["roi"]
                Gc.graph["roi_um"] = args.roi_um
                Gc.graph["spacing_nm"] = spacing

                necks = np.array([d["neck_nm"] for _, _, d in Gc.edges(data=True)])
                vols = np.array([d["volume_nm3"] for _, d in Gc.nodes(data=True)])
                row.update(
                    lcc_nodes=Gc.number_of_nodes(), lcc_edges=Gc.number_of_edges(),
                    lcc_face_lo=len(lo), lcc_face_hi=len(hi),
                    neck_p10_nm=float(np.percentile(necks, 10)),
                    neck_p50_nm=float(np.percentile(necks, 50)),
                    mean_degree=float(2 * Gc.number_of_edges() / Gc.number_of_nodes()),
                    pore_vol_median_nm3=float(np.median(vols)),
                )
                with open(os.path.join(
                        NETS, f"{key}__{r['roi']}__{args.roi_um}um.pkl"), "wb") as f:
                    pickle.dump(Gc, f)
                print(f"  ROI {r['roi']}: pfrac={diag['perc_percolating_frac']:.4f}  "
                      f"pores={diag['n_pores']:,} throats={diag['n_throats']:,}  "
                      f"LCC=({Gc.number_of_nodes():,},{Gc.number_of_edges():,}) "
                      f"deg={row['mean_degree']:.2f}  "
                      f"neck p10={row['neck_p10_nm']:.0f} p50={row['neck_p50_nm']:.0f} nm  "
                      f"faces=({len(lo)},{len(hi)})  [{t2-t1:.0f}s]")
                if i == 0:
                    dest = os.path.join(OUT, f"phase3_snow_overlay_{key}.png")
                    save_snow_overlay(key, r["roi"], extras["mask"],
                                      extras["regions"], spacing, Gc, dest)
                    print(f"    [overlay] {os.path.basename(dest)}")
            else:
                print(f"  ROI {r['roi']}: NO NETWORK — {diag.get('reason')}")
            rows.append(row)
            del ni, extras

    df = pd.DataFrame(rows)
    dest = os.path.join(OUT, f"phase3_snow_{args.roi_um}um_rmax{args.r_max}.csv")
    df.to_csv(dest, index=False)
    print(f"\n[saved] {dest}")
    cols = [c for c in ("sample", "roi", "perc_percolating_frac", "lcc_nodes",
                        "lcc_edges", "mean_degree", "neck_p10_nm", "neck_p50_nm",
                        "t_snow_s") if c in df.columns]
    if not df.empty:
        print(df[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
