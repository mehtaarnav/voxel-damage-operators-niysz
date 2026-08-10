"""
Bounded null-autopsy for the p10-group null. SECONDARY AND EXPLANATORY --
does not change the frozen primary outcome, introduces no new primary
outcome, and changes no D4 parameter (p_erode=0.35, expand_vox=1 as frozen).

Primary diagnostic question: at n_rounds=8 (the last spanning state), do the
widened structures still retain a measurable lower-tail geometric advantage
over base?

D1 remaining Ni volume at n_rounds 0/8/9
D2 lower-tail thickness proxy (EDT p10/p25 over remaining Ni) at 0/8
D3 spanning-cluster size at n_rounds 8
D4loc failure-step localization: voxels removed 8->9, and their EDT values
    measured in the n=8 mask (one base seed + one high seed)
"""
from __future__ import annotations

import os, sys, time
import numpy as np, pandas as pd
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from cmlib.damage import add_ysz_pore, apply_d4, rebuild_ternary  # noqa: E402
from cmlib.percolation import percolation_summary  # noqa: E402
from cmlib.synth import (  # noqa: E402
    build_mass_conservative_structure, draw_valid_base_widths,
    max_clip_widths, platform_v2_lattice_geometry, rasterize,
)

OUT = os.path.join(ROOT, "out", "platform_v2")
R_VOX, PITCH, NLAT_Z, NLAT_XY, MARGIN, JITTER = 12.1, 32, 6, 4, 8, 0.15
VOXEL_NM, GEOM_SEED = 20.0, 999
FW, WR, NR, MINR = 0.20, (4, 6), (12, 20), 2.5
STRUCT_SEEDS = [0, 1, 2, 3, 4]
DAMAGE_SEEDS = [200, 201, 202, 203, 204]
P_ERODE, EXPAND_VOX = 0.35, 1            # FROZEN
YSZ_REST = 0.388 / (0.388 + 0.362)


def edt_tail(mask):
    """Lower-tail thickness proxy: EDT percentiles over remaining Ni, nm.
    Relative diagnostic only -- NOT a gate metric."""
    if not mask.any():
        return np.nan, np.nan
    d = ndi.distance_transform_edt(mask)[mask] * VOXEL_NM
    return float(np.percentile(d, 10)), float(np.percentile(d, 25))


def main():
    t0 = time.time()
    c, p, shape = platform_v2_lattice_geometry(
        NLAT_Z, NLAT_XY, PITCH, R_VOX, MARGIN, JITTER,
        np.random.default_rng(GEOM_SEED))
    qual = pd.read_csv(os.path.join(OUT, "qualification_run.csv"))
    rows, loc_rows = [], []

    for _, q in qual.iterrows():
        seed, mode, tr = int(q.seed), q["mode"], q.target_ratio
        T = 0.0 if mode == "base" else float(q.intended_T_vox)
        w, _, _, _ = draw_valid_base_widths(len(p), seed, FW, WR, NR, MINR)
        ni0 = (rasterize(c, p, R_VOX, w, shape) if T <= 0 else
               build_mass_conservative_structure(
                   c, p, shape, R_VOX, w, max_clip_widths(w, T))["ni_mask"])
        _, ysz = add_ysz_pore(ni0, seed=seed * 1000 + 7,
                              ysz_frac_of_rest=YSZ_REST)
        v0 = int(ni0.sum())
        t10_0, t25_0 = edt_tail(ni0)

        for ds in DAMAGE_SEEDS:
            rec = dict(mode=mode, nominal_target_ratio=tr, struct_seed=seed,
                       damage_seed=ds, achieved_p10=q.achieved_ratio,
                       ni_vox_n0=v0, edt_p10_n0=t10_0, edt_p25_n0=t25_0)
            masks = {}
            for n in (8, 9):
                dm, _ = apply_d4(ni0, ysz, n, P_ERODE, EXPAND_VOX, ds)
                masks[n] = dm
                v = int(dm.sum())
                rec[f"ni_vox_n{n}"] = v
                rec[f"ni_retained_n{n}"] = v / v0 if v0 else np.nan
                if n == 8:
                    a, b = edt_tail(dm)
                    rec["edt_p10_n8"], rec["edt_p25_n8"] = a, b
                    ps = percolation_summary(dm, axis=0, connectivity=6,
                                             check_other_axes=False)
                    rec["P_span_n8"] = ps["P_span"]
                    rec["span_vox_n8"] = ps["P_span"] * v
                    rec["span_frac_of_orig_n8"] = ps["P_span"] * v / v0
                else:
                    rec["P_span_n9"] = percolation_summary(
                        dm, axis=0, connectivity=6,
                        check_other_axes=False)["P_span"]
            rows.append(rec)

            # D4loc: one base + one high seed, first damage seed only
            if ds == DAMAGE_SEEDS[0] and seed == 0 and (
                    mode == "base" or np.isclose(tr, 2.0)):
                removed = masks[8] & ~masks[9]
                dt8 = ndi.distance_transform_edt(masks[8]) * VOXEL_NM
                vals = dt8[removed]
                allv = dt8[masks[8]]
                loc_rows.append(dict(
                    group=("base" if mode == "base" else "high_2.00x"),
                    n_removed=int(removed.sum()),
                    frac_of_n8=float(removed.sum() / masks[8].sum()),
                    rem_p10=float(np.percentile(vals, 10)),
                    rem_p50=float(np.percentile(vals, 50)),
                    rem_p90=float(np.percentile(vals, 90)),
                    rem_mean=float(vals.mean()),
                    all_p10=float(np.percentile(allv, 10)),
                    all_p50=float(np.percentile(allv, 50)),
                    all_p90=float(np.percentile(allv, 90)),
                    all_mean=float(allv.mean())))
        print(f"  {mode:10s} tr={tr:<5} seed={seed}: v0={v0} "
              f"edt_p10_n0={t10_0:.1f}nm  [{time.time()-t0:.0f}s]", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "null_autopsy.csv"), index=False)
    loc = pd.DataFrame(loc_rows)
    if len(loc):
        loc.to_csv(os.path.join(OUT, "null_autopsy_localization.csv"),
                   index=False)

    g = df.groupby(["mode", "nominal_target_ratio"])
    print("\n" + "=" * 76)
    print("D1  REMAINING Ni VOLUME")
    print("=" * 76)
    t1 = g.agg(achieved_p10=("achieved_p10", "mean"),
               vox_n0=("ni_vox_n0", "mean"),
               vox_n8=("ni_vox_n8", "mean"), vox_n8_sd=("ni_vox_n8", "std"),
               ret_n8=("ni_retained_n8", "mean"),
               ret_n8_sd=("ni_retained_n8", "std"),
               vox_n9=("ni_vox_n9", "mean"),
               ret_n9=("ni_retained_n9", "mean")).round(4)
    print(t1.to_string())

    print("\n" + "=" * 76)
    print("D2  LOWER-TAIL THICKNESS PROXY (EDT percentiles over remaining Ni, nm)")
    print("=" * 76)
    t2 = g.agg(p10_n0=("edt_p10_n0", "mean"), p10_n8=("edt_p10_n8", "mean"),
               p10_n8_sd=("edt_p10_n8", "std"),
               p25_n0=("edt_p25_n0", "mean"),
               p25_n8=("edt_p25_n8", "mean")).round(3)
    print(t2.to_string())

    print("\n" + "=" * 76)
    print("D3  SPANNING-CLUSTER SIZE AT n_rounds=8 (last spanning state)")
    print("=" * 76)
    t3 = g.agg(P_span_n8=("P_span_n8", "mean"),
               span_vox_n8=("span_vox_n8", "mean"),
               span_vox_sd=("span_vox_n8", "std"),
               span_frac_orig=("span_frac_of_orig_n8", "mean"),
               span_frac_sd=("span_frac_of_orig_n8", "std"),
               P_span_n9=("P_span_n9", "mean")).round(4)
    print(t3.to_string())

    if len(loc):
        print("\n" + "=" * 76)
        print("D4loc  FAILURE-STEP LOCALIZATION (voxels removed 8->9;")
        print("       their EDT measured in the n=8 mask, nm)")
        print("=" * 76)
        print(loc.round(2).to_string(index=False))

    print(f"\n[saved] null_autopsy.csv"
          + (", null_autopsy_localization.csv" if len(loc) else "")
          + f"   total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
