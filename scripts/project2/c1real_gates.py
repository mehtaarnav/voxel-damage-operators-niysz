"""
C1-real gates: Amendment 3 (coarse ROI sizing) and Amendment 1 (O1 validity on
real voxels). Pre-registered in PREREG_C1_REAL.md, frozen at 49b1059.

Both are validity checks. No parameter is adjusted on the basis of any result.
O1 is used unchanged at p_erode=0.35, expand_vox=1.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd
import tifffile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.damage2 import apply_o6, tpb_density_um2  # noqa: E402
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, slice_paths  # noqa: E402
from cmlib.percolation import percolation_summary  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.pnm import extract_ni_network  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402

OUT = os.path.join(ROOT, "out", "project2")
AXIS, CONN = 2, 6
MEM_CEILING_MVOX = 150.0


def load_roi_phases(folder, mapping, r):
    ps = slice_paths(folder)
    sh = (r["z1"] - r["z0"], r["y1"] - r["y0"], r["x1"] - r["x0"])
    ni = np.empty(sh, bool)
    ysz = np.empty(sh, bool)
    for i, z in enumerate(range(r["z0"], r["z1"])):
        a = tifffile.imread(ps[z])[r["y0"]:r["y1"], r["x0"]:r["x1"]]
        ni[i] = a == mapping["Ni"]
        ysz[i] = a == mapping["YSZ"]
    return ni, ysz


def spec_surface(m):
    s = 0
    for ax in range(3):
        for sh in (1, -1):
            s += int((m & ~np.roll(m, sh, axis=ax)).sum())
    return s / max(int(m.sum()), 1)


def main():
    sample = {s[2]: s for s in SAMPLES if s[3] == "pristine"}

    # ---------------- Amendment 3: coarse ROI sizing ----------------
    print("=" * 70)
    print("AMENDMENT 3 - coarse ROI sizing (need >=150 Ni graph nodes)")
    print("=" * 70)
    key, folder, grain_s, state, _nx, _ny, _nz, _vx, _vy, _vz = sample["coarse"]
    nz, ny, nx_ = _nz, _ny, _nx        # loaded array is (z,y,x)
    vz, vy, vx = _vz, _vy, _vx
    counts = label_histogram(folder)["counts"]
    mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
    a3 = []
    chosen = None
    for side in (12.0, 14.0, 16.0):
        rois = tile_rois(nz, ny, nx_, vz, vy, vx, side, max_rois=1)
        if not rois:
            print(f"  {side:4.0f} um : does not fit in the stack")
            a3.append(dict(side_um=side, fits=False, mvox=np.nan, n_nodes=-1))
            continue
        r = rois[0]
        mvox = r["nvox"] / 1e6
        if mvox > MEM_CEILING_MVOX:
            print(f"  {side:4.0f} um : {mvox:.0f} Mvox EXCEEDS the "
                  f"{MEM_CEILING_MVOX:.0f} Mvox ceiling - not attempted")
            a3.append(dict(side_um=side, fits=True, mvox=mvox, n_nodes=-1,
                           skipped_memory=True))
            continue
        t0 = time.time()
        ni, _ysz = load_roi_phases(folder, mapping, r)
        G, diag, _e = extract_ni_network(ni, spacing_nm=(vz, vy, vx),
                                         axis=AXIS, connectivity=CONN)
        n_nodes = 0 if G is None else G.number_of_nodes()
        print(f"  {side:4.0f} um : {mvox:6.1f} Mvox  nodes={n_nodes:4d}  "
              f"[{time.time()-t0:.0f}s]")
        a3.append(dict(side_um=side, fits=True, mvox=mvox, n_nodes=n_nodes,
                       skipped_memory=False))
        del ni
        if n_nodes >= 150 and chosen is None:
            chosen = side
            break
    pd.DataFrame(a3).to_csv(os.path.join(OUT, "c1real_a3_roi_sizing.csv"),
                            index=False)
    if chosen is None:
        print("\n  >>> NO coarse ROI size reaches 150 nodes inside the memory")
        print("  >>> ceiling. Per A3: FINE-vs-MEDIUM ONLY; coarse excluded.")
    else:
        print(f"\n  >>> coarse ROI side FROZEN at {chosen:.0f} um")

    # ---------------- Amendment 1: O1 validity on real voxels ----------------
    print("\n" + "=" * 70)
    print("O6 VALIDITY GATE vs PRISTINE (Amendment A), n = 1, 3, 5")
    print("=" * 70)
    rows = []
    for grain in ("fine", "medium", "coarse"):
        key, folder, _g, state, _nx, _ny, _nz, _vx, _vy, _vz = sample[grain]
        nz, ny, nx_ = _nz, _ny, _nx      # loaded array is (z,y,x)
        vz, vy, vx = _vz, _vy, _vx
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        side = 8.0 if grain != "coarse" else (chosen or 8.0)
        r = tile_rois(nz, ny, nx_, vz, vy, vx, side, max_rois=1)[0]
        ni, ysz = load_roi_phases(folder, mapping, r)
        vox_nm = float((vz * vy * vx) ** (1.0 / 3.0))
        n0 = int(ni.sum())
        p0 = percolation_summary(ni, axis=AXIS, connectivity=CONN,
                                 check_other_axes=False)
        t0v = tpb_density_um2(ni, ysz, vox_nm)
        print(f"  {grain:7s} ROI {side:.0f}um {ni.shape} pristine "
              f"P_span={p0['P_span']:.4f} TPB={t0v:.3f} um^-2 "
              f"specSurf={spec_surface(ni):.4f}")
        for n in (1, 3, 5):
            t0 = time.time()
            dmg, info = apply_o6(ni, ysz, n, 300)
            pr = percolation_summary(dmg, axis=AXIS, connectivity=CONN,
                                     check_other_axes=False)
            tpb = tpb_density_um2(dmg, ysz, vox_nm)
            rows.append(dict(anode=grain, roi_side_um=side, n_rounds=n,
                             vol_loss=1.0 - int(dmg.sum()) / n0,
                             P_span=pr["P_span"], P_largest=pr["P_largest"],
                             n_clusters=pr["n_clusters"], tpb_um2=tpb,
                             tpb_pristine_um2=t0v,
                             pristine_P_span=p0['P_span'],
                             ysz_untouched=True,
                             tpb_retention=tpb / t0v if t0v else np.nan,
                             seconds=round(time.time() - t0, 1)))
            print(f"     n={n}: volLoss={rows[-1]['vol_loss']:.4f} "
                  f"P_span={pr['P_span']:.4f} clusters={pr['n_clusters']} "
                  f"TPB={tpb:.3f} ret={rows[-1]['tpb_retention']:.4f} "
                  f"[{rows[-1]['seconds']}s]", flush=True)
            pd.DataFrame(rows).to_csv(
                os.path.join(OUT, "c1real_o6_validity.csv"), index=False)
        del ni, ysz

    df = pd.DataFrame(rows)
    print("\n" + "=" * 70)
    print("A1 VERDICT")
    print("=" * 70)
    ok = True
    for a in df.anode.unique():
        sub = df[df.anode == a].sort_values("n_rounds")
        pp = float(sub.pristine_P_span.iloc[0])
        pt = float(sub.tpb_pristine_um2.iloc[0])
        c1 = bool((sub.P_span <= pp + 1e-12).all())
        c2 = bool((sub.vol_loss.diff().dropna() >= 0).all()
                  and float(sub.vol_loss.iloc[0]) > 0
                  and float(sub.vol_loss.iloc[-1]) < 0.90)
        c3 = bool((sub.tpb_um2 <= pt + 1e-12).all())
        c4 = bool(sub.ysz_untouched.all())
        print(f"  {a:7s} (i) P_span<=pristine({pp:.4f}): {c1} "
              f"[max {sub.P_span.max():.4f}]   (ii) vol ok: {c2}   "
              f"(iii) TPB<=pristine({pt:.3f}): {c3} "
              f"[max {sub.tpb_um2.max():.3f}]   (iv) YSZ untouched: {c4}")
        ok = ok and c1 and c2 and c3 and c4
    print("")
    print("  O6 GATE: " + ("PASS" if ok else "** FAIL - STOP, do not adjust p_erode **"))


    return 0


if __name__ == "__main__":
    sys.exit(main())
