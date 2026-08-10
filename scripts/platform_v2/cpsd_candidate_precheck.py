"""
Step 2 — CORRELATION PRE-CHECK before building any candidate metric.

Three candidate size measures, cheapest-first:
  (a) local_thickness r50 (`cpsd_r50max`) at sizes=300 -- the incumbent, now
      known unstable/non-monotone. Expensive; already computed for 2 seeds
      (base+high) by the ladder, reused here rather than recomputed.
  (b) raw-EDT volume-weighted mean, NO binning: 2*mean(EDT[Ni])*voxel_nm.
      Each voxel is unit volume, so a plain mean over Ni voxels IS the
      volume-weighted mean. Has no `sizes` parameter at all, so it cannot
      have the binned-median instability -- but that is only useful if it
      actually tracks true local thickness.
  (c) generator_radius_deviation: (r_final - r_base)/r_base, exact and free
      from the qualification CSV. No image processing, no estimator noise --
      potentially the most trustworthy anchor of the three, though it only
      exists for synthetic structures (there is no "generator radius" for a
      real anode, so it can never be the cross-comparison metric).

(b) and (c) are computed for ALL 15 synthetic structures so the correlation
has 10 widened points, not 2. (a) is the 4-point anchor.
"""
from __future__ import annotations

import os, sys, time
import numpy as np, pandas as pd
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from cmlib.synth import (  # noqa: E402
    build_mass_conservative_structure, draw_valid_base_widths,
    max_clip_widths, platform_v2_lattice_geometry, rasterize,
)

OUT = os.path.join(ROOT, "out", "platform_v2")
R_VOX, PITCH, NLAT_Z, NLAT_XY, MARGIN, JITTER = 12.1, 32, 6, 4, 8, 0.15
VOXEL_NM, GEOM_SEED = 20.0, 999
FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO = 0.20, (4, 6), (12, 20), 2.5

# from out/platform_v2/cpsd_convergence.csv (sizes=300), the ladder's top rung
LT300 = {(0, "base"): 458.49, (0, "high"): 414.65,
         (1, "base"): 458.08, (1, "high"): 423.52}


def main():
    t0 = time.time()
    centres, pairs, shape = platform_v2_lattice_geometry(
        NLAT_Z, NLAT_XY, PITCH, R_VOX, MARGIN, JITTER,
        np.random.default_rng(GEOM_SEED))
    qual = pd.read_csv(os.path.join(OUT, "qualification_run.csv"))

    rows = []
    for _, r in qual.iterrows():
        seed = int(r.seed)
        T = 0.0 if r["mode"] == "base" else float(r["intended_T_vox"])
        w, _, _, _ = draw_valid_base_widths(len(pairs), seed, FRAC_WEAK,
                                            WEAK_RANGE, NORMAL_RANGE, MIN_RATIO)
        if T <= 0:
            mask = rasterize(centres, pairs, R_VOX, w, shape)
        else:
            mask = build_mass_conservative_structure(
                centres, pairs, shape, R_VOX, w, max_clip_widths(w, T))["ni_mask"]
        edt = ndi.distance_transform_edt(mask)          # voxel units, no binning
        d_edt = 2.0 * float(edt[mask].mean()) * VOXEL_NM
        rows.append(dict(mode=r["mode"], target_ratio=r.target_ratio, seed=seed,
                         achieved_ratio=r.achieved_ratio,
                         r_final_vox=r.r_final_vox,
                         d_edt_vwmean_nm=d_edt,
                         cpsd_sizes25_nm=r.d_cPSD_r50max_nm))
        print(f"  {r['mode']:10s} tr={r.target_ratio:<5} seed={seed}: "
              f"d_edt_vwmean={d_edt:7.2f} nm  [{time.time()-t0:.0f}s]")
        del mask, edt

    df = pd.DataFrame(rows)
    base = df[df["mode"] == "base"].set_index("seed")
    for col, new in (("d_edt_vwmean_nm", "edt_dev_pct"),
                     ("r_final_vox", "genrad_dev_pct"),
                     ("cpsd_sizes25_nm", "cpsd25_dev_pct")):
        df[new] = df.apply(lambda x: 100 * (x[col] - base.loc[x.seed, col])
                           / base.loc[x.seed, col], axis=1)
    df["lt300_nm"] = df.apply(
        lambda x: LT300.get((x.seed, "base" if x["mode"] == "base"
                             else ("high" if np.isclose(x.target_ratio, 2.0)
                                   else None)), np.nan), axis=1)
    df["lt300_dev_pct"] = df.apply(
        lambda x: (100 * (x.lt300_nm - LT300[(x.seed, "base")])
                   / LT300[(x.seed, "base")])
        if (x.seed in (0, 1) and np.isfinite(x.lt300_nm)) else np.nan, axis=1)
    df.to_csv(os.path.join(OUT, "cpsd_candidate_precheck.csv"), index=False)

    wid = df[df["mode"] != "base"]
    print("\n" + "=" * 74)
    print("DEVIATIONS per widened structure (n=10)")
    print("=" * 74)
    print(wid[["target_ratio", "seed", "achieved_ratio", "edt_dev_pct",
               "genrad_dev_pct", "cpsd25_dev_pct", "lt300_dev_pct"]]
          .round(3).to_string(index=False))

    print("\n" + "=" * 74)
    print("CORRELATION (n=10 widened structures)")
    print("=" * 74)
    a, b = wid.edt_dev_pct.values, wid.genrad_dev_pct.values
    print(f"  raw-EDT dev  vs  generator-radius dev : "
          f"r = {np.corrcoef(a, b)[0,1]:+.4f}")
    print(f"    typical offset (EDT - genrad): "
          f"{np.mean(a-b):+.3f} pp  (sd {np.std(a-b, ddof=1):.3f})")
    c = wid.cpsd25_dev_pct.values
    print(f"  raw-EDT dev  vs  cpsd(sizes=25) dev   : "
          f"r = {np.corrcoef(a, c)[0,1]:+.4f}   (sizes=25 is the BROKEN "
          f"resolution -- shown only for contrast)")

    anchor = wid[(wid.seed.isin([0, 1])) & np.isclose(wid.target_ratio, 2.0)]
    print("\n  4-point anchor vs true local thickness at sizes=300:")
    for _, x in anchor.iterrows():
        print(f"    seed {int(x.seed)}: lt300={x.lt300_dev_pct:+.2f}%  "
              f"raw-EDT={x.edt_dev_pct:+.2f}%  "
              f"genrad={x.genrad_dev_pct:+.2f}%  "
              f"| EDT-lt300 offset={x.edt_dev_pct-x.lt300_dev_pct:+.2f} pp")
    print(f"\n[saved] cpsd_candidate_precheck.csv   total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
