"""
Two jobs, in order:

(A) CONVERGENCE: extend the sizes ladder to 300 on the same two synthetic
    seeds already tested at 25/50/100/200, and check whether the metric has
    stopped moving between SUCCESSIVE resolutions (target: <0.5 pp change).
    Nothing is locked until this passes.

(B) RESOLUTION-MATCHED RECOMPUTE: the platform-v2 qualification run computed
    synthetic c-PSD at sizes=25 (the default at the time). The real-ROI
    variability study used sizes=100. Those are NOT comparable -- sizes=25
    quantizes c-PSD into ~6%-wide bins. This recomputes all 15 synthetic
    structures at the converged sizes so both sides of the comparison use
    the same resolution.

Structures are reconstructed deterministically from qualification_run.csv
(same seeds + recorded max-clip thresholds), the same reconstruction pattern
already verified bit-identical in the E0 spike.
"""
from __future__ import annotations

import os, sys, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from cmlib.particles import cpsd_r50max  # noqa: E402
from cmlib.synth import (  # noqa: E402
    build_mass_conservative_structure, draw_valid_base_widths,
    max_clip_widths, platform_v2_lattice_geometry, rasterize,
)

OUT = os.path.join(ROOT, "out", "platform_v2")
R_VOX, PITCH, NLAT_Z, NLAT_XY, MARGIN, JITTER = 12.1, 32, 6, 4, 8, 0.15
VOXEL_NM, GEOM_SEED = 20.0, 999
FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO = 0.20, (4, 6), (12, 20), 2.5
CONVERGED_SIZES = 100      # provisional; (A) must confirm before this is locked
LADDER = [25, 50, 100, 200, 300]


def geom():
    return platform_v2_lattice_geometry(NLAT_Z, NLAT_XY, PITCH, R_VOX, MARGIN,
                                        JITTER, np.random.default_rng(GEOM_SEED))


def build(centres, pairs, shape, seed, threshold):
    w, _, _, _ = draw_valid_base_widths(len(pairs), seed, FRAC_WEAK,
                                        WEAK_RANGE, NORMAL_RANGE, MIN_RATIO)
    if threshold <= 0:
        return rasterize(centres, pairs, R_VOX, w, shape)
    return build_mass_conservative_structure(
        centres, pairs, shape, R_VOX, w, max_clip_widths(w, threshold))["ni_mask"]


def main():
    t0 = time.time()
    centres, pairs, shape = geom()
    qual = pd.read_csv(os.path.join(OUT, "qualification_run.csv"))

    # ---------------------------------------------------------------- (A)
    print("=" * 76)
    print("(A) CONVERGENCE — successive-resolution change, target <0.5 pp")
    print("=" * 76)
    rows = []
    for seed in (0, 1):
        hi_T = float(qual[(qual["mode"] == "lower_tail")
                          & (np.isclose(qual.target_ratio, 2.0))
                          & (qual.seed == seed)].iloc[0]["intended_T_vox"])
        b = build(centres, pairs, shape, seed, 0.0)
        h = build(centres, pairs, shape, seed, hi_T)
        print(f"\n--- seed {seed} (high-ratio T={hi_T:.2f}vox) ---")
        prev = None
        for s in LADDER:
            cb = cpsd_r50max(b, VOXEL_NM, sizes=s)["d_cPSD_r50max_nm"]
            ch = cpsd_r50max(h, VOXEL_NM, sizes=s)["d_cPSD_r50max_nm"]
            dev = 100 * (ch - cb) / cb
            step = "" if prev is None else f"  step_change={abs(dev-prev):.2f} pp"
            print(f"  sizes={s:4d}: base={cb:7.2f} high={ch:7.2f} "
                  f"dev={dev:+6.2f}%{step}")
            rows.append(dict(seed=seed, sizes=s, base=cb, high=ch, dev_pct=dev,
                             step_change_pp=(np.nan if prev is None
                                             else abs(dev - prev))))
            prev = dev
        del b, h
    conv = pd.DataFrame(rows)
    conv.to_csv(os.path.join(OUT, "cpsd_convergence.csv"), index=False)

    print("\n  successive-step changes (pp):")
    ok_at = {}
    for seed, g in conv.groupby("seed"):
        g = g.sort_values("sizes")
        for _, r in g.iterrows():
            if np.isfinite(r.step_change_pp):
                mark = "OK" if r.step_change_pp < 0.5 else "still moving"
                print(f"    seed {seed}: ->{int(r.sizes):4d}  "
                      f"{r.step_change_pp:5.2f} pp  {mark}")
                ok_at.setdefault(int(r.sizes), []).append(r.step_change_pp < 0.5)
    converged_at = None
    for s in LADDER[1:]:
        if ok_at.get(s) and all(ok_at[s]):
            converged_at = s
            break
    print(f"\n  first ladder step where BOTH seeds change <0.5 pp: "
          f"{converged_at if converged_at else 'NONE — metric still moving'}")

    # ---------------------------------------------------------------- (B)
    use = CONVERGED_SIZES
    print("\n" + "=" * 76)
    print(f"(B) RESOLUTION-MATCHED RECOMPUTE of all 15 synthetic structures "
          f"at sizes={use}")
    print(f"    (real-ROI study used sizes=100; qualification run used 25)")
    print("=" * 76)
    recs = []
    for _, r in qual.iterrows():
        T = 0.0 if r["mode"] == "base" else float(r["intended_T_vox"])
        m = build(centres, pairs, shape, int(r.seed), T)
        st = cpsd_r50max(m, VOXEL_NM, sizes=use)
        recs.append(dict(mode=r["mode"], target_ratio=r.target_ratio,
                         seed=int(r.seed), achieved_ratio=r.achieved_ratio,
                         cpsd_sizes25=r.d_cPSD_r50max_nm,
                         cpsd_converged=st["d_cPSD_r50max_nm"]))
        print(f"  {r['mode']:10s} tr={r.target_ratio:<5} seed={int(r.seed)}: "
              f"sizes25={r.d_cPSD_r50max_nm:7.1f} -> "
              f"sizes{use}={st['d_cPSD_r50max_nm']:7.1f} nm  "
              f"[{time.time()-t0:.0f}s]")
        del m
    df = pd.DataFrame(recs)
    base = df[df["mode"] == "base"].set_index("seed")
    df["cpsd_dev_pct_converged"] = df.apply(
        lambda r: 100 * (r.cpsd_converged - base.loc[r.seed, "cpsd_converged"])
        / base.loc[r.seed, "cpsd_converged"], axis=1)
    df["cpsd_dev_pct_sizes25"] = df.apply(
        lambda r: 100 * (r.cpsd_sizes25 - base.loc[r.seed, "cpsd_sizes25"])
        / base.loc[r.seed, "cpsd_sizes25"], axis=1)
    df.to_csv(os.path.join(OUT, "cpsd_recomputed_converged.csv"), index=False)

    print("\n  c-PSD deviation, sizes=25 (as gated) vs converged:")
    for tr, g in df[df["mode"] != "base"].groupby("target_ratio"):
        print(f"    target {tr}: sizes25 "
              f"[{g.cpsd_dev_pct_sizes25.min():+.2f}%, "
              f"{g.cpsd_dev_pct_sizes25.max():+.2f}%]   converged "
              f"[{g.cpsd_dev_pct_converged.min():+.2f}%, "
              f"{g.cpsd_dev_pct_converged.max():+.2f}%]")
    print(f"\n[saved] cpsd_convergence.csv, cpsd_recomputed_converged.csv"
          f"   total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
