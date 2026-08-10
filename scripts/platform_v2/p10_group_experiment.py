"""
Platform v2 — p10-group damage-bisection experiment.

Runs ONLY under the design frozen in preregistration.md §0c–§0g. Key
constraints honoured here:

- §0g/1: 5 damage seeds per structure (compute is cheap: the base-only run
  did 15 brackets in 66 s). Damage seeds {200..204} are INDEPENDENT of
  structure seeds {0..4}. Per-damage-seed midpoints, per-structure means,
  within-structure variance, across-structure variance, and group means with
  seed-level spread are all reported. No single-seed comparison anywhere.
- §0g/2: p_erode=0.35, expand_vox=1 FROZEN. Only n_rounds varies.
- §0c/C: achieved measured p10 ratio is the analysis variable; nominal labels
  are reported alongside but never used for the contrast.
- §0f/4: the causal-interpretation obligation (radius-shrink confound) is
  applied to the OUTPUT, not to the run itself -- generator radius deviation
  is carried through per structure so the confound can be inspected directly.

Outcome variable: the damage intensity (n_rounds) at which the structure
loses percolation, bracketed per (structure, damage seed) to width 1 by the
frozen bisection. HIGHER midpoint = survives more damage = better retention.
"""
from __future__ import annotations

import os, sys, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from cmlib.api import compute_percolation  # noqa: E402
from cmlib.damage import add_ysz_pore, apply_d4, rebuild_ternary  # noqa: E402
from cmlib.synth import (  # noqa: E402
    build_mass_conservative_structure, draw_valid_base_widths,
    max_clip_widths, platform_v2_lattice_geometry, rasterize,
)
from cmlib.synthvol import LABELS  # noqa: E402

OUT = os.path.join(ROOT, "out", "platform_v2")
R_VOX, PITCH, NLAT_Z, NLAT_XY, MARGIN, JITTER = 12.1, 32, 6, 4, 8, 0.15
VOXEL_NM, GEOM_SEED = 20.0, 999
FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO = 0.20, (4, 6), (12, 20), 2.5
STRUCT_SEEDS = [0, 1, 2, 3, 4]
DAMAGE_SEEDS = [200, 201, 202, 203, 204]      # independent of structure seeds
P_ERODE, EXPAND_VOX = 0.35, 1                 # FROZEN (§0g/2)
BRACKET_LO, BRACKET_HI, BRACKET_MAX = 1, 20, 64
YSZ_FRAC_OF_REST = 0.388 / (0.388 + 0.362)


def main():
    t0 = time.time()
    centres, pairs, shape = platform_v2_lattice_geometry(
        NLAT_Z, NLAT_XY, PITCH, R_VOX, MARGIN, JITTER,
        np.random.default_rng(GEOM_SEED))
    qual = pd.read_csv(os.path.join(OUT, "qualification_run.csv"))

    print("=" * 78)
    print("PLATFORM V2 — p10-GROUP DAMAGE-BISECTION EXPERIMENT")
    print("=" * 78)
    print(f"structure seeds {STRUCT_SEEDS}  x  damage seeds {DAMAGE_SEEDS} "
          f"(5 per structure, §0g/1)")
    print(f"p_erode={P_ERODE} expand_vox={EXPAND_VOX} FROZEN (§0g/2); "
          f"only n_rounds varies\n")

    rows = []
    for _, q in qual.iterrows():
        seed = int(q.seed)
        mode, tr = q["mode"], q.target_ratio
        T = 0.0 if mode == "base" else float(q.intended_T_vox)
        w, _, _, _ = draw_valid_base_widths(len(pairs), seed, FRAC_WEAK,
                                            WEAK_RANGE, NORMAL_RANGE, MIN_RATIO)
        if T <= 0:
            ni = rasterize(centres, pairs, R_VOX, w, shape)
        else:
            ni = build_mass_conservative_structure(
                centres, pairs, shape, R_VOX, w, max_clip_widths(w, T))["ni_mask"]
        _, ysz = add_ysz_pore(ni, seed=seed * 1000 + 7,
                              ysz_frac_of_rest=YSZ_FRAC_OF_REST)

        for dseed in DAMAGE_SEEDS:
            def spans(n):
                dm, _ = apply_d4(ni, ysz, n, P_ERODE, EXPAND_VOX, dseed)
                return compute_percolation(rebuild_ternary(dm, ysz), LABELS,
                                           phase="Ni", axis=0)["P_span"] > 0.0
            lo, hi, note = BRACKET_LO, BRACKET_HI, "ok"
            if not spans(lo):
                mid, note = np.nan, "lost at bracket floor"
            else:
                while spans(hi):
                    lo, hi = hi, hi * 2
                    if hi > BRACKET_MAX:
                        break
                if hi > BRACKET_MAX:
                    mid, note = np.nan, f"still spanning at {BRACKET_MAX}"
                else:
                    while hi - lo > 1:
                        m = (lo + hi) // 2
                        if spans(m):
                            lo = m
                        else:
                            hi = m
                    mid = (lo + hi) / 2.0
            rows.append(dict(mode=mode, nominal_target_ratio=tr,
                             struct_seed=seed, damage_seed=dseed,
                             achieved_p10_ratio=q.achieved_ratio,
                             r_final_vox=q.r_final_vox,
                             n_lo=lo, n_hi=hi, midpoint=mid, note=note))
        got = [r["midpoint"] for r in rows[-len(DAMAGE_SEEDS):]]
        print(f"  {mode:10s} tr={tr:<5} seed={seed} "
              f"(achieved p10 {q.achieved_ratio:.2f}x): "
              f"midpoints={[f'{g:.1f}' for g in got]}  "
              f"mean={np.nanmean(got):.2f}  [{time.time()-t0:.0f}s]")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "p10_group_experiment.csv"), index=False)
    analyse(df, qual)
    print(f"\n[saved] p10_group_experiment.csv   total {time.time()-t0:.0f}s")
    return 0


def analyse(df, qual):
    ok = df[df.note == "ok"].copy()
    print("\n" + "=" * 78)
    print("VARIANCE DECOMPOSITION (§0g/1)")
    print("=" * 78)
    per = ok.groupby(["mode", "nominal_target_ratio", "struct_seed"]).agg(
        achieved_p10=("achieved_p10_ratio", "first"),
        mid_mean=("midpoint", "mean"), mid_sd=("midpoint", "std"),
        n=("midpoint", "count")).reset_index()
    print(per.round(3).to_string(index=False))

    within = ok.groupby(["mode", "nominal_target_ratio", "struct_seed"]
                        ).midpoint.var(ddof=1).mean()
    across = per.groupby(["mode", "nominal_target_ratio"]).mid_mean.var(ddof=1).mean()
    print(f"\n  mean WITHIN-structure damage-seed variance : {within:.4f}"
          f"  (sd {np.sqrt(within):.3f} rounds)")
    print(f"  mean ACROSS-structure variance of the means: {across:.4f}"
          f"  (sd {np.sqrt(across):.3f} rounds)")
    print(f"  ratio within/across = {within/across:.2f}"
          if across > 0 else "  across-variance is zero")

    print("\n" + "=" * 78)
    print("GROUP MEANS (analysis variable = ACHIEVED p10 ratio, §0c/C)")
    print("=" * 78)
    grp = per.groupby(["mode", "nominal_target_ratio"]).agg(
        achieved_p10=("achieved_p10", "mean"),
        group_mean_midpoint=("mid_mean", "mean"),
        spread_of_structure_means=("mid_mean", "std"),
        n_structures=("mid_mean", "count")).reset_index()
    print(grp.round(3).to_string(index=False))

    base = grp[grp["mode"] == "base"].group_mean_midpoint.iloc[0]
    print(f"\n  base group mean transition = {base:.3f} rounds")
    print("  differences vs base (positive = survives MORE damage):")
    decisive = True
    for _, g in grp[grp["mode"] != "base"].iterrows():
        d = g.group_mean_midpoint - base
        flag = "" if abs(d) >= 1.0 else "   <-- BELOW the 1.0-round threshold (§0g/1)"
        if abs(d) < 1.0:
            decisive = False
        print(f"    achieved p10 {g.achieved_p10:.2f}x : "
              f"{d:+.3f} rounds{flag}")

    print("\n" + "=" * 78)
    print("PRE-REGISTERED INTERPRETATION")
    print("=" * 78)
    if decisive:
        print("  At least one group difference >= 1.0 damage round.")
    else:
        print("  ALL group differences are SMALLER than 1.0 damage round.")
        print("  Per §0g/1 this comparison is NOT interpretable as a branch")
        print("  decision at 3 seeds; it was already run at the maximum 5")
        print("  damage seeds, so no further seed-completion is available.")
        print("  -> The honest reading is NO RESOLVABLE EFFECT at this")
        print("     structure count and damage-seed count. This is the")
        print("     'no effect' branch, NOT Path B, NOT a negative result on")
        print("     the scientific hypothesis (see §0f/4 and the E0 rules).")
    print("\n  §0f/4 CONFOUND NOTE (mandatory): the high-ratio structures also")
    print("  carry a 5.17-6.79% sphere-radius shrink. Any effect seen here")
    print("  could not be attributed to neck widening alone without the")
    print("  matched-shrink control. Radius deviation per structure:")
    for _, g in per[per["mode"] != "base"].iterrows():
        rb = qual[(qual["mode"] == "base")
                  & (qual.seed == g.struct_seed)].r_final_vox.iloc[0]
        rf = qual[(qual["mode"] != "base")
                  & (qual.seed == g.struct_seed)
                  & np.isclose(qual.target_ratio, g.nominal_target_ratio)
                  ].r_final_vox.iloc[0]
        print(f"    seed {g.struct_seed} p10 {g.achieved_p10:.2f}x: "
              f"radius dev {100*(rf-rb)/rb:+.2f}%  "
              f"transition {g.mid_mean:.2f}+-{g.mid_sd:.2f}")


if __name__ == "__main__":
    sys.exit(main())
