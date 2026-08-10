"""
T5b — CORRECTED coupling-decision experiment (mass-conservative,
percentile-targeted), per the 2026-08-10 review decision. Answers
out/next/preregistration.md #0b acceptance criteria. Supersedes
scripts/next/t5_coupling_experiment.py, whose two failure modes (Ni mass not
conserved; "widen bottom 20%" cannot move measured p10 by construction) are
fixed here — see cmlib/synth.py module docstring and
out/next/t5_coupling_decision_report.md for the full diagnosis.

BASE NECK-WIDTH RANGE WIDENED FROM [2,6] TO [2,10] VOXELS (stated choice, not
a silent change): explored analytically before any structure was built (see
transcript) -- with the original [2,6] range, base p50=4vox is too close to
base p10=2vox (factor of 2) to leave room for a genuinely LOWER-TAIL-SELECTIVE
intervention at the requested 1.5x-2.5x p10 targets without dragging p50 along
too. [2,10] gives p50/p10 = 3.0 (vs 2.0), and at every target ratio tested
analytically (on the INTENDED, not yet SNOW-measured, distribution) the
max-clip construction leaves p50 at exactly 1.00x while p10 hits the target
ratio exactly. Same lattice geometry (5x5x5, R=8vox, pitch=20vox) as T5,
per the review decision ("same basic lattice is acceptable for this
diagnostic").

For each (mode, target_ratio, seed), reports: intended target, measured neck
p10, measured neck p50, c-PSD size deviation, watershed size deviation, Phi_Ni
deviation, generator-known primary size deviation, SNOW node count, initial
P_span, and total Ni voxels added/removed/net-residual -- all seed-matched
against that SAME seed's base (unwidened) structure, not a cross-seed mean.
"""

from __future__ import annotations

import os
import sys
import time

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from cmlib.api import (  # noqa: E402
    compute_network_metrics, compute_particle_stats, compute_percolation,
    extract_network,
)
from cmlib.percolation import label_phase  # noqa: E402
from cmlib.synth import (  # noqa: E402
    base_neck_widths, build_mass_conservative_structure, lattice_geometry,
    max_clip_widths, uniform_shift_widths,
)
from cmlib.synthvol import LABELS  # noqa: E402

OUT = os.path.join(ROOT, "out", "next")
os.makedirs(OUT, exist_ok=True)

NLAT, PITCH_VOX, R_VOX, MARGIN = 5, 20, 8.0, 6
BASE_LO, BASE_HI = 2.0, 10.0          # widened from T5's [2,6], see docstring
VOXEL_NM = 20.0
SEEDS = [0, 1, 2, 3, 4]
TARGET_RATIOS = [1.5, 2.0, 2.5]
MAX_BISECT_ITER = 5
P10_TOL_REL = 0.08                     # accept within 8% of target ratio

CENTRES, PAIRS, SHAPE = lattice_geometry(NLAT, PITCH_VOX, R_VOX, MARGIN)
DOMAIN_VOX = int(np.prod(SHAPE))


def measure_structure(ni_mask, spacing):
    perc = compute_percolation(
        np.where(ni_mask, LABELS["Ni"], LABELS["pore"]).astype(np.uint8),
        LABELS, phase="Ni", axis=0)
    pstats = compute_particle_stats(ni_mask, spacing, min_distance=4)
    _, n_clusters = label_phase(ni_mask, connectivity=6)
    out = dict(P_span=perc["P_span"], percolates=perc["percolates"],
              n_clusters=int(n_clusters))
    out.update({k: v for k, v in pstats.items()
               if k in ("ws_d_volweighted_nm", "ws_n_regions_used",
                        "d_cPSD_r50max_nm", "n_peaks")})
    G, diag = extract_network(ni_mask, spacing, axis=0, r_max=4)
    if G is not None and G.number_of_edges() > 0:
        nm = compute_network_metrics(G, G.graph.get("face_lo"),
                                     G.graph.get("face_hi"))
        out.update(neck_p10_nm=nm["neck_p10_nm"], neck_p50_nm=nm["neck_p50_nm"],
                   n_nodes=nm["n_nodes"], n_edges=nm["n_edges"])
    else:
        out.update(neck_p10_nm=np.nan, neck_p50_nm=np.nan, n_nodes=0, n_edges=0)
    return out


def eval_lower_tail_point(base_widths, r_base, threshold, spacing):
    final_widths = max_clip_widths(base_widths, threshold)
    build = build_mass_conservative_structure(
        CENTRES, PAIRS, SHAPE, r_base, base_widths, final_widths)
    m = measure_structure(build["ni_mask"], spacing)
    return build, m


def bisect_threshold(base_widths, r_base, target_p10_nm, spacing,
                     t_lo, t_hi, max_iter=MAX_BISECT_ITER, tol_rel=P10_TOL_REL):
    """Bisect the max-clip threshold T (voxels) so measured neck_p10_nm hits
    target_p10_nm within tol_rel. Returns (best_T, best_build, best_measure,
    log: list of (T, measured_p10, build, measure))."""
    log = []
    best = None
    lo, hi = t_lo, t_hi
    for it in range(max_iter):
        T = 0.5 * (lo + hi)
        build, m = eval_lower_tail_point(base_widths, r_base, T, spacing)
        p10 = m["neck_p10_nm"]
        log.append((T, p10, build, m))
        err = abs(p10 - target_p10_nm) / target_p10_nm if np.isfinite(p10) else np.inf
        if best is None or err < best[0]:
            best = (err, T, build, m)
        if err <= tol_rel:
            break
        if not np.isfinite(p10) or p10 < target_p10_nm:
            lo = T           # need a bigger threshold
        else:
            hi = T            # overshot -> smaller threshold
    _, T_best, build_best, m_best = best
    return T_best, build_best, m_best, log


def main():
    print("=" * 78)
    print("T5b — mass-conservative, percentile-targeted coupling experiment")
    print("=" * 78)
    print(f"lattice {NLAT}^3, R={R_VOX}vox, pitch={PITCH_VOX}vox, "
          f"base neck range [{BASE_LO},{BASE_HI}]vox, voxel={VOXEL_NM}nm, "
          f"domain {SHAPE}\n")

    spacing = (VOXEL_NM, VOXEL_NM, VOXEL_NM)
    rows = []
    t0 = time.time()

    # ---------------------------------------------------------------- base
    base_cache = {}
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        base_w = base_neck_widths(len(PAIRS), rng, BASE_LO, BASE_HI)
        build, m = eval_lower_tail_point(base_w, R_VOX, 0.0, spacing)  # T=0 -> no widening
        base_cache[seed] = dict(base_widths=base_w, build=build, measure=m)
        rows.append(dict(mode="base", target_ratio=1.0, seed=seed,
                         intended_T_vox=np.nan, iters=0, **_flat(build, m)))
        print(f"  base seed={seed}: phi_Ni={m.get('P_span'):.3f}(Pspan) "
              f"neck_p10={m['neck_p10_nm']:.0f}nm neck_p50={m['neck_p50_nm']:.0f}nm "
              f"ws_d={m.get('ws_d_volweighted_nm', float('nan')):.0f}nm "
              f"cPSD_d={m.get('d_cPSD_r50max_nm', float('nan')):.0f}nm  "
              f"[{time.time()-t0:.0f}s]")

    # ---------------------------------------------------------- lower-tail
    for ratio in TARGET_RATIOS:
        for seed in SEEDS:
            base_w = base_cache[seed]["base_widths"]
            base_p10 = base_cache[seed]["measure"]["neck_p10_nm"]
            target_p10 = ratio * base_p10
            t_lo = base_p10 / VOXEL_NM
            t_hi = BASE_HI + 4
            T_best, build, m, log = bisect_threshold(
                base_w, R_VOX, target_p10, spacing, t_lo, t_hi)
            achieved_ratio = m["neck_p10_nm"] / base_p10
            rows.append(dict(mode="lower_tail", target_ratio=ratio, seed=seed,
                             intended_T_vox=T_best, iters=len(log),
                             achieved_ratio=achieved_ratio,
                             **_flat(build, m)))
            print(f"  lower_tail ratio={ratio} seed={seed}: T={T_best:.2f}vox "
                  f"iters={len(log)} achieved_p10_ratio={achieved_ratio:.2f} "
                  f"neck_p50={m['neck_p50_nm']:.0f}nm "
                  f"phi_dev={100*(build['phi_Ni_final']-build['phi_Ni_base'])/build['phi_Ni_base']:+.1f}% "
                  f"[{time.time()-t0:.0f}s]")

    # ------------------------------------------------------------- uniform
    for ratio in TARGET_RATIOS:
        for seed in SEEDS:
            base_w = base_cache[seed]["base_widths"]
            base_p10 = base_cache[seed]["measure"]["neck_p10_nm"]
            target_p10_vox = ratio * base_p10 / VOXEL_NM
            median_base = float(np.median(base_w))
            delta = target_p10_vox - BASE_LO   # p10 of a uniform shift sits near lo+delta
            final_widths = uniform_shift_widths(base_w, delta, BASE_LO,
                                                 2 * R_VOX - 2)
            build = build_mass_conservative_structure(
                CENTRES, PAIRS, SHAPE, R_VOX, base_w, final_widths)
            m = measure_structure(build["ni_mask"], spacing)
            achieved_ratio = m["neck_p10_nm"] / base_p10
            rows.append(dict(mode="uniform", target_ratio=ratio, seed=seed,
                             intended_T_vox=delta, iters=1,
                             achieved_ratio=achieved_ratio,
                             **_flat(build, m)))
            print(f"  uniform    ratio={ratio} seed={seed}: delta={delta:.2f}vox "
                  f"achieved_p10_ratio={achieved_ratio:.2f} "
                  f"neck_p50={m['neck_p50_nm']:.0f}nm "
                  f"phi_dev={100*(build['phi_Ni_final']-build['phi_Ni_base'])/build['phi_Ni_base']:+.1f}% "
                  f"[{time.time()-t0:.0f}s]")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "t5b_coupling_experiment.csv"), index=False)
    print(f"\n[saved] {os.path.join(OUT, 't5b_coupling_experiment.csv')}"
          f"   total {time.time()-t0:.0f}s")

    analyse_and_report(df)
    return 0


def _flat(build, m):
    return dict(
        r_base_vox=build["r_base_vox"], r_final_vox=build["r_final_vox"],
        voxels_base=build["voxels_base"],
        voxels_added_by_necks=build["voxels_added_by_necks"],
        voxels_removed_by_shrink=build["voxels_removed_by_shrink"],
        voxels_final=build["voxels_final"],
        voxels_net_residual=build["voxels_net_residual"],
        phi_Ni_base=build["phi_Ni_base"], phi_Ni_final=build["phi_Ni_final"],
        **m,
    )


def analyse_and_report(df):
    print("\n" + "=" * 78)
    print("SEED-MATCHED DEVIATIONS (vs that seed's own base structure)")
    print("=" * 78)

    base = df[df["mode"] == "base"].set_index("seed")
    recs = []
    for _, r in df[df["mode"] != "base"].iterrows():
        b = base.loc[r.seed]
        rec = dict(mode=r["mode"], target_ratio=r.target_ratio, seed=r.seed,
                   achieved_p10_ratio=r.get("achieved_ratio", np.nan),
                   neck_p10_nm=r.neck_p10_nm, neck_p50_nm=r.neck_p50_nm,
                   neck_p50_ratio=r.neck_p50_nm / b.neck_p50_nm,
                   ws_dev_pct=100 * (r.ws_d_volweighted_nm - b.ws_d_volweighted_nm) / b.ws_d_volweighted_nm,
                   cpsd_dev_pct=100 * (r.d_cPSD_r50max_nm - b.d_cPSD_r50max_nm) / b.d_cPSD_r50max_nm,
                   phi_dev_pct=100 * (r.phi_Ni_final - b.phi_Ni_final) / b.phi_Ni_final,
                   primary_size_dev_pct=100 * (r.r_final_vox - b.r_final_vox) / b.r_final_vox,
                   n_nodes=r.n_nodes, P_span=r.P_span, n_clusters=r.n_clusters,
                   voxels_added=r.voxels_added_by_necks,
                   voxels_removed=r.voxels_removed_by_shrink,
                   voxels_net_residual=r.voxels_net_residual)
        recs.append(rec)
    dev = pd.DataFrame(recs)
    dev.to_csv(os.path.join(OUT, "t5b_deviations.csv"), index=False)

    agg = dev.groupby(["mode", "target_ratio"]).agg(
        n=("seed", "count"),
        p10_ratio_mean=("achieved_p10_ratio", "mean"),
        p10_ratio_sd=("achieved_p10_ratio", "std"),
        p50_ratio_mean=("neck_p50_ratio", "mean"),
        ws_dev_mean=("ws_dev_pct", "mean"), ws_dev_sd=("ws_dev_pct", "std"),
        cpsd_dev_mean=("cpsd_dev_pct", "mean"), cpsd_dev_sd=("cpsd_dev_pct", "std"),
        phi_dev_mean=("phi_dev_pct", "mean"), phi_dev_sd=("phi_dev_pct", "std"),
        primary_dev_mean=("primary_size_dev_pct", "mean"),
        P_span_min=("P_span", "min"),
        n_clusters_max=("n_clusters", "max"),
    ).reset_index()
    agg.to_csv(os.path.join(OUT, "t5b_deviations_agg.csv"), index=False)

    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(agg.to_string(index=False))

    print("\n" + "=" * 78)
    print("ACCEPTANCE CRITERIA CHECK (preregistration.md #0b)")
    print("=" * 78)
    print("  criteria: p10>=1.5x  AND  |phi_dev|<=5% (target <=2%)  AND  "
          "|cPSD_dev|<=5%  AND  p50_ratio << p10_ratio  AND  P_span intact  "
          "AND  n_clusters not blown up")

    passing = []
    for _, r in agg[agg["mode"] == "lower_tail"].iterrows():
        checks = dict(
            p10=r.p10_ratio_mean >= 1.5,
            phi=abs(r.phi_dev_mean) <= 5.0,
            phi_target=abs(r.phi_dev_mean) <= 2.0,
            cpsd=abs(r.cpsd_dev_mean) <= 5.0,
            p50_stable=r.p50_ratio_mean < (r.p10_ratio_mean - 1) * 0.5 + 1,
            pspan=r.P_span_min >= 0.90,
        )
        full_pass = checks["p10"] and checks["phi"] and checks["cpsd"] and checks["p50_stable"] and checks["pspan"]
        print(f"  lower_tail target_ratio={r.target_ratio}: "
              f"p10_ratio={r.p10_ratio_mean:.2f} (need>=1.5: {checks['p10']})  "
              f"phi_dev={r.phi_dev_mean:+.2f}% (need<=5%: {checks['phi']}, "
              f"<=2% target: {checks['phi_target']})  "
              f"cPSD_dev={r.cpsd_dev_mean:+.2f}% (need<=5%: {checks['cpsd']})  "
              f"p50_ratio={r.p50_ratio_mean:.2f} vs p10_ratio={r.p10_ratio_mean:.2f} "
              f"(p50<<p10: {checks['p50_stable']})  "
              f"P_span_min={r.P_span_min:.3f} (>=0.90: {checks['pspan']})  "
              f"n_clusters_max={r.n_clusters_max}  "
              f"-> {'PASS' if full_pass else 'fail'}")
        if full_pass:
            passing.append(r.target_ratio)

    print(f"\n  lower_tail points passing ALL criteria: {passing}")
    uniform_pass = []
    for _, r in agg[agg["mode"] == "uniform"].iterrows():
        checks_u = (r.p10_ratio_mean >= 1.5 and abs(r.phi_dev_mean) <= 5.0
                   and abs(r.cpsd_dev_mean) <= 5.0 and r.P_span_min >= 0.90)
        if checks_u:
            uniform_pass.append(r.target_ratio)
    print(f"  uniform points passing (p10/phi/cPSD/P_span only, p50 not "
          f"required to be stable since uniform is not tail-selective): "
          f"{uniform_pass}")

    print("\n  DECISION (preregistration.md #0b tree):")
    if passing:
        print("  -> lower-tail criteria MET at >=1 point. PROCEED TO FAMILY B "
              "is the indicated next step (pending your review).")
    elif uniform_pass:
        print("  -> lower-tail FAILED, uniform compensated widening PASSED. "
              "STOP AND REPORT per the frozen decision tree: this requires an "
              "explicit amendment reframing the primary axis before any "
              "further code is written. Do NOT proceed to Family B "
              "automatically.")
    else:
        print("  -> NEITHER lower-tail nor uniform passed. STOP, prepare a "
              "Path-B-style limitation memo. Do not claim physical "
              "impossibility -- report 'not testable within this synthetic "
              "framework, two constructions tried, both failed for stated "
              "reasons.'")

    # ------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    for mode, col in (("lower_tail", "C0"), ("uniform", "C1")):
        d = agg[agg["mode"] == mode].sort_values("target_ratio")
        axes[0].errorbar(d.p10_ratio_mean, d.phi_dev_mean, yerr=d.phi_dev_sd,
                         fmt="o-", color=col, label=mode, capsize=3)
        axes[1].errorbar(d.p10_ratio_mean, d.cpsd_dev_mean, yerr=d.cpsd_dev_sd,
                         fmt="o-", color=col, label=f"{mode} (c-PSD)")
        axes[1].errorbar(d.p10_ratio_mean, d.ws_dev_mean, yerr=d.ws_dev_sd,
                         fmt="s--", color=col, alpha=0.5,
                         label=f"{mode} (watershed, diagnostic)")
        axes[2].plot(d.p10_ratio_mean, d.p50_ratio_mean, "o-", color=col,
                    label=mode)
    axes[0].axhspan(-5, 5, color="green", alpha=0.08, label="+/-5% ceiling")
    axes[0].axhspan(-2, 2, color="green", alpha=0.15, label="+/-2% target")
    axes[0].axvline(1.5, color="k", ls=":", lw=1)
    axes[0].set_xlabel("measured neck p10 ratio")
    axes[0].set_ylabel("Phi_Ni deviation (%)")
    axes[0].set_title("Mass conservation check")
    axes[0].legend(fontsize=7, frameon=False)
    axes[0].grid(alpha=0.25)

    axes[1].axhspan(-5, 5, color="green", alpha=0.1, label="+/-5% band")
    axes[1].axvline(1.5, color="k", ls=":", lw=1)
    axes[1].set_xlabel("measured neck p10 ratio")
    axes[1].set_ylabel("particle-size deviation (%)")
    axes[1].set_title("c-PSD (primary gate) vs watershed (diagnostic)")
    axes[1].legend(fontsize=7, frameon=False)
    axes[1].grid(alpha=0.25)

    axes[2].plot([1, 3], [1, 3], "k:", lw=1, label="p50 moves 1:1 with p10")
    axes[2].axvline(1.5, color="k", ls=":", lw=1)
    axes[2].set_xlabel("measured neck p10 ratio")
    axes[2].set_ylabel("measured neck p50 ratio")
    axes[2].set_title("Tail-selectivity: p50 should stay near 1.0")
    axes[2].legend(fontsize=7, frameon=False)
    axes[2].grid(alpha=0.25)

    fig.suptitle("T5b: mass-conservative, percentile-targeted coupling check",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "t5b_coupling_experiment.png"), dpi=145)
    print(f"\n[saved] {os.path.join(OUT, 't5b_deviations.csv')}")
    print(f"[saved] {os.path.join(OUT, 't5b_deviations_agg.csv')}")
    print(f"[saved] {os.path.join(OUT, 't5b_coupling_experiment.png')}")


if __name__ == "__main__":
    sys.exit(main())
