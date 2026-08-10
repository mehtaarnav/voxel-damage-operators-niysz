"""
Q1: measure REAL Holzer/Pecho c-PSD variability, to test whether the +-5%
synthetic c-PSD ceiling is grounded.

c-PSD was NEVER computed on the real data -- the real-data study used
watershed particle sizing (phase4d) and SNOW chamber diameters (phase4c).
`cmlib.particles.cpsd_r50max` was added later for the synthetic study. So
this script computes it on the real ROIs for the first time, using the SAME
8 um ROI tiling as phase3/phase4 so the numbers are comparable to the rest
of the real-data study.

SCOPE: medium and coarse anodes only -- platform v2 is explicitly scoped to
medium/coarse-like geometry, and fine's ROIs are 69 Mvoxel each (3.5x the
coarse ROIs), which would dominate runtime for an out-of-scope comparison.

sizes=100: chosen after finding sizes=25 (the previous default) is too coarse
-- it quantizes c-PSD into ~6% bins, comparable to the +-5% gate itself, and
gave unreliable per-structure values on the synthetic set (a seed reported
-3.13% at sizes=25 vs -9.18% at sizes=200). 100 is where the synthetic
deviations had substantially converged.
"""
from __future__ import annotations

import os, sys, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, load_subvolume  # noqa: E402
from cmlib.particles import cpsd_r50max  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402

OUT = os.path.join(ROOT, "out", "platform_v2")
SIZES = 100
SAMPLES_IN_SCOPE = ["medium_pre", "coarse_pre"]


def main():
    t0 = time.time()
    rows = []
    for key in SAMPLES_IN_SCOPE:
        _, folder, grain, state, nx_, ny_, nz_, vx, vy, vz = \
            [s for s in SAMPLES if s[0] == key][0]
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        rois = tile_rois(nz_, ny_, nx_, vz, vy, vx, 8.0)
        vox_geo = float((vx * vy * vz) ** (1 / 3))
        print(f"\n=== {key} ({grain}) — {len(rois)} ROIs, voxel_geo={vox_geo:.2f}nm ===")
        for r in rois:
            vol = load_subvolume(folder, r["z0"], r["z1"], r["y0"], r["y1"],
                                 r["x0"], r["x1"])
            ni = vol == mapping["Ni"]
            del vol
            st = cpsd_r50max(ni, vox_geo, sizes=SIZES)
            rows.append(dict(sample=key, grain=grain, roi=r["roi"],
                             nvox=r["nvox"], **st))
            print(f"  {r['roi']}: d_cPSD_r50max={st['d_cPSD_r50max_nm']:8.1f} nm "
                  f"(r50={st['cpsd_r50_nm']:.1f})  [{time.time()-t0:.0f}s]")
            del ni

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "real_cpsd_variability.csv"), index=False)

    print("\n" + "=" * 72)
    print("REAL c-PSD VARIABILITY (the question: is +-5% grounded?)")
    print("=" * 72)
    for key, g in df.groupby("sample"):
        v = g.d_cPSD_r50max_nm
        print(f"\n{key}: n_roi={len(v)}")
        print(f"  mean={v.mean():.1f} nm  sd={v.std(ddof=1):.1f} nm  "
              f"CV={100*v.std(ddof=1)/v.mean():.1f}%")
        print(f"  min={v.min():.1f}  max={v.max():.1f}  "
              f"full spread={100*(v.max()-v.min())/v.mean():.1f}% of mean")
        dev = 100 * (v - v.mean()) / v.mean()
        print(f"  per-ROI deviation from own-anode mean: "
              f"[{dev.min():+.1f}%, {dev.max():+.1f}%]")
        print(f"  ROIs beyond +-5% of own-anode mean: "
              f"{int((dev.abs() > 5).sum())}/{len(v)}")

    if df.sample.nunique() > 1:
        m = df.groupby("sample").d_cPSD_r50max_nm.mean()
        print(f"\nanode-to-anode: {dict(m.round(1))}  "
              f"medium->coarse change = "
              f"{100*(m.get('coarse_pre',np.nan)-m.get('medium_pre',np.nan))/m.get('medium_pre',np.nan):+.1f}%")
    print(f"\n[saved] {os.path.join(OUT,'real_cpsd_variability.csv')}  "
          f"total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
