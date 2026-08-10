"""
Platform v2, milestones A + B + C.

A. YSZ/pore placement on the qualified Platform-v2 Ni structures (Ni phase
   untouched; verify Phi_Ni holds at ~0.2502, P_span unchanged, TPB nonzero
   and plausible).
B. D4 re-validation on Platform-v2 ternary structures (reconstruction
   integrity, Ni untouched by placement, YSZ untouched by damage, Ni->pore
   consistency, metric coherence, no NaN, sensible largest-component
   behaviour). E0's intensity grid is explicitly NOT assumed to carry over.
C. Damage-calibration bisection on BASE structures ONLY -- widened
   structures are not inspected. Frozen per-structure bisection: integer
   n_rounds in [1,20], narrow until bracket width <=1, expand bracket only
   if the transition lies outside it (never change D4 parameters or
   criteria).

Stops after C. Does NOT run the p10-group damage experiment.
"""
from __future__ import annotations

import os, sys, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from cmlib.api import compute_percolation, compute_tpb  # noqa: E402
from cmlib.damage import add_ysz_pore, apply_d4, rebuild_ternary  # noqa: E402
from cmlib.synth import (  # noqa: E402
    draw_valid_base_widths, platform_v2_lattice_geometry, rasterize,
)
from cmlib.synthvol import LABELS, volume_fractions_from_volume  # noqa: E402

OUT = os.path.join(ROOT, "out", "platform_v2")
R_VOX, PITCH, NLAT_Z, NLAT_XY, MARGIN, JITTER = 12.1, 32, 6, 4, 8, 0.15
VOXEL_NM, GEOM_SEED = 20.0, 999
FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO = 0.20, (4, 6), (12, 20), 2.5
SEEDS = [0, 1, 2, 3, 4]
SPACING = (VOXEL_NM,) * 3

# medium anode: phi_Ni .250 / phi_YSZ .388 / phi_pore .362
YSZ_FRAC_OF_REST = 0.388 / (0.388 + 0.362)

# D4 params -- UNCHANGED from E0 (only the intensity grid is re-derived)
P_ERODE, EXPAND_VOX = 0.35, 1
DAMAGE_SEEDS = [100, 101, 102]
BRACKET_LO, BRACKET_HI, BRACKET_MAX = 1, 20, 64


def build_base(centres, pairs, shape, seed):
    w, _, _, _ = draw_valid_base_widths(len(pairs), seed, FRAC_WEAK,
                                        WEAK_RANGE, NORMAL_RANGE, MIN_RATIO)
    return rasterize(centres, pairs, R_VOX, w, shape)


def pspan(vol):
    return compute_percolation(vol, LABELS, phase="Ni", axis=0)["P_span"]


def main():
    t0 = time.time()
    centres, pairs, shape = platform_v2_lattice_geometry(
        NLAT_Z, NLAT_XY, PITCH, R_VOX, MARGIN, JITTER,
        np.random.default_rng(GEOM_SEED))
    qual = pd.read_csv(os.path.join(OUT, "qualification_run.csv"))
    qbase = qual[qual["mode"] == "base"].set_index("seed")

    # ================================================== A. YSZ/pore placement
    print("=" * 78)
    print("A. YSZ/PORE PLACEMENT on qualified Platform-v2 Ni structures")
    print("=" * 78)
    print(f"target YSZ share of non-Ni remainder = {YSZ_FRAC_OF_REST:.4f} "
          f"(medium anode 0.388/(0.388+0.362))\n")
    store, rowsA = {}, []
    for seed in SEEDS:
        ni = build_base(centres, pairs, shape, seed)
        phi_ni_before = float(ni.mean())
        recorded = float(qbase.loc[seed, "phi_Ni_final"])
        ps_before = pspan(np.where(ni, LABELS["Ni"], LABELS["pore"]).astype(np.uint8))

        vol, ysz = add_ysz_pore(ni, seed=seed * 1000 + 7,
                                ysz_frac_of_rest=YSZ_FRAC_OF_REST)
        phi = volume_fractions_from_volume(vol)
        ps_after = pspan(vol)
        tpb = compute_tpb(vol, LABELS, SPACING)["tpb_density_um-2"]

        ni_preserved = bool((vol == LABELS["Ni"]).sum() == ni.sum()
                            and np.array_equal(vol == LABELS["Ni"], ni))
        assert ni_preserved, f"seed {seed}: placement modified the Ni phase"
        assert abs(ps_after - ps_before) < 1e-12, f"seed {seed}: P_span changed"
        assert tpb > 0, f"seed {seed}: TPB is zero -- pathological"

        store[seed] = (ni, ysz)
        rowsA.append(dict(seed=seed, phi_Ni=phi["phi_Ni"],
                          phi_Ni_recorded_qual=recorded,
                          phi_Ni_delta_vs_qual=phi["phi_Ni"] - recorded,
                          phi_YSZ=phi["phi_YSZ"], phi_pore=phi["phi_pore"],
                          ysz_share_of_rest=phi["phi_YSZ"] /
                          (phi["phi_YSZ"] + phi["phi_pore"]),
                          P_span_before=ps_before, P_span_after=ps_after,
                          tpb_density=tpb, ni_phase_preserved=ni_preserved))
        print(f"  seed {seed}: Phi_Ni={phi['phi_Ni']:.4f} "
              f"(qual {recorded:.4f}, delta {phi['phi_Ni']-recorded:+.2e})  "
              f"Phi_YSZ={phi['phi_YSZ']:.4f}  Phi_pore={phi['phi_pore']:.4f}  "
              f"YSZ/rest={phi['phi_YSZ']/(phi['phi_YSZ']+phi['phi_pore']):.4f}  "
              f"P_span {ps_before:.3f}->{ps_after:.3f}  "
              f"TPB={tpb:.4f} um^-2  [{time.time()-t0:.0f}s]")
    dfA = pd.DataFrame(rowsA)
    dfA.to_csv(os.path.join(OUT, "ternary_placement.csv"), index=False)
    print(f"\n  Phi_Ni across seeds: mean={dfA.phi_Ni.mean():.4f} "
          f"max|delta vs qualification|={dfA.phi_Ni_delta_vs_qual.abs().max():.2e}")
    print(f"  TPB density: {dfA.tpb_density.min():.4f}-{dfA.tpb_density.max():.4f} um^-2 "
          f"(real anodes 1.07-2.65; see REPORT.md)")
    print("  Ni preserved / P_span unchanged / TPB nonzero: ALL PASS")

    # ================================================== B. D4 re-validation
    print("\n" + "=" * 78)
    print("B. D4 RE-VALIDATION on Platform-v2 ternary geometry")
    print("=" * 78)
    rowsB = []
    for seed in SEEDS:
        ni, ysz = store[seed]
        ni2 = build_base(centres, pairs, shape, seed)
        recon_ok = bool(np.array_equal(ni, ni2))          # reconstruction integrity
        dmg, d = apply_d4(ni, ysz, n_rounds=5, p_erode=P_ERODE,
                          expand_vox=EXPAND_VOX, seed=100)
        volD = rebuild_ternary(dmg, ysz)
        ysz_untouched = bool(np.array_equal(volD == LABELS["YSZ"], ysz))
        lost = ni & ~dmg
        lost_to_pore = bool(np.all(volD[lost] == LABELS["pore"])) if lost.any() else True
        perc = compute_percolation(volD, LABELS, phase="Ni", axis=0)
        tpbD = compute_tpb(volD, LABELS, SPACING)["tpb_density_um-2"]
        finite = bool(np.isfinite([perc["P_span"], perc["P_reach"], tpbD]).all())
        rowsB.append(dict(seed=seed, reconstruction_identical=recon_ok,
                          ysz_untouched_by_damage=ysz_untouched,
                          ni_loss_becomes_pore=lost_to_pore,
                          all_metrics_finite=finite,
                          P_span_post=perc["P_span"], P_reach_post=perc["P_reach"],
                          n_clusters_post=perc["n_clusters"], tpb_post=tpbD,
                          **d))
        print(f"  seed {seed}: recon={recon_ok} ysz_untouched={ysz_untouched} "
              f"Ni_loss->pore={lost_to_pore} finite={finite} | "
              f"n_comp_before_island_removal={d['n_components_before_island_removal']} "
              f"P_span={perc['P_span']:.3f} TPB={tpbD:.4f}  [{time.time()-t0:.0f}s]")
    dfB = pd.DataFrame(rowsB)
    dfB.to_csv(os.path.join(OUT, "d4_revalidation.csv"), index=False)
    allB = bool(dfB[["reconstruction_identical", "ysz_untouched_by_damage",
                     "ni_loss_becomes_pore", "all_metrics_finite"]].all().all())
    print(f"\n  D4 RE-VALIDATION: {'ALL CHECKS PASS' if allB else 'FAILURE'}")

    # ============================== C. bisection on BASE structures only
    print("\n" + "=" * 78)
    print("C. DAMAGE-CALIBRATION BISECTION — BASE structures only")
    print("=" * 78)
    print(f"  frozen procedure: integer n_rounds, bracket [{BRACKET_LO},"
          f"{BRACKET_HI}], expand-only, narrow to width<=1; "
          f"p_erode={P_ERODE} expand_vox={EXPAND_VOX} UNCHANGED\n")
    rowsC = []
    for seed in SEEDS:
        ni, ysz = store[seed]
        for dseed in DAMAGE_SEEDS:
            def spans(n):
                dm, _ = apply_d4(ni, ysz, n, P_ERODE, EXPAND_VOX, dseed)
                return pspan(rebuild_ternary(dm, ysz)) > 0.0

            lo, hi, expands = BRACKET_LO, BRACKET_HI, 0
            if not spans(lo):
                rowsC.append(dict(seed=seed, damage_seed=dseed, n_lo=np.nan,
                                  n_hi=lo, midpoint=np.nan, n_evals=1,
                                  note="already lost at bracket floor"))
                print(f"  seed {seed} dseed {dseed}: lost at n={lo} (floor)")
                continue
            while spans(hi):                      # expand-only, never re-tune
                lo, hi, expands = hi, hi * 2, expands + 1
                if hi > BRACKET_MAX:
                    break
            n_ev = 2 + expands
            if hi > BRACKET_MAX:
                rowsC.append(dict(seed=seed, damage_seed=dseed, n_lo=lo,
                                  n_hi=np.nan, midpoint=np.nan, n_evals=n_ev,
                                  note=f"still spanning at n={BRACKET_MAX}"))
                print(f"  seed {seed} dseed {dseed}: still spanning at "
                      f"n={BRACKET_MAX}")
                continue
            while hi - lo > 1:
                mid = (lo + hi) // 2
                n_ev += 1
                if spans(mid):
                    lo = mid
                else:
                    hi = mid
            rowsC.append(dict(seed=seed, damage_seed=dseed, n_lo=lo, n_hi=hi,
                              midpoint=(lo + hi) / 2.0, n_evals=n_ev, note="ok"))
            print(f"  seed {seed} dseed {dseed}: transition bracketed "
                  f"[{lo},{hi}] midpoint={(lo+hi)/2.0:.1f}  "
                  f"({n_ev} evals)  [{time.time()-t0:.0f}s]")
    dfC = pd.DataFrame(rowsC)
    dfC.to_csv(os.path.join(OUT, "d4_bisection_base.csv"), index=False)

    ok = dfC[dfC.note == "ok"]
    print("\n  TRANSITION-INTENSITY DISTRIBUTION (base structures):")
    if len(ok):
        print(f"    midpoint: mean={ok.midpoint.mean():.2f}  "
              f"sd={ok.midpoint.std(ddof=1):.2f}  "
              f"range=[{ok.midpoint.min():.1f}, {ok.midpoint.max():.1f}]  "
              f"n={len(ok)}")
        per = ok.groupby("seed").midpoint.agg(["mean", "std", "count"])
        print("\n    per seed:")
        print(per.round(2).to_string())
    print(f"\n[saved] ternary_placement.csv, d4_revalidation.csv, "
          f"d4_bisection_base.csv   total {time.time()-t0:.0f}s")
    print("\nSTOPPING per instruction: p10-group damage experiment NOT run.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
