"""
FAMILY B DISORDERED PILOT — required before Family C / damage models, per
out/next/preregistration.md #0d.

Scope (frozen): 5 seeds, {base, 1.5x, 2.0x} (the primary envelope per #0c
amendment A; 2.5x is a stress case only and is NOT run here), disordered
(jittered-lattice) Ni structures, full per-seed diagnostic table against the
STRICT criteria in #0c amendment B (no aggregate-mean pass criterion --
every verdict is per-seed), base-distribution validity enforced BEFORE
widening with every rejected attempt logged.

STOP after this pilot and report. Do not build ternary YSZ/pore placement or
any damage model until this is reviewed (explicit instruction).

Geometry: `cmlib.synth.jittered_lattice_geometry`, NLAT=4 (64 spheres, 144
candidate neck pairs), R_VOX=14 (particle diameter 28 voxels = 560 nm at
20 nm/voxel, inside the requested 24-32 voxel range), pitch=28 voxels,
domain 85x160x160 (~2.18 Mvoxel) -- a 160-class pilot, per #0d.

Base neck-width distribution: `cmlib.synth.mixture_neck_widths`, a mixture of
a "normal" population (13-22 voxels) and a genuinely narrow-tail
subpopulation (4-6 voxels, weight 20%), giving population p50/p10 ~= 3.2,
comfortably above the pre-registered validity floor of 2.5, with a hard
minimum of 4 voxels satisfying the resolution requirement. Every seed's draw
is checked against the validity floor BEFORE widening
(`cmlib.synth.draw_valid_base_widths`); rejected attempts are logged.
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
    base_distribution_stats, build_mass_conservative_structure,
    draw_valid_base_widths, jittered_lattice_geometry, max_clip_widths,
)
from cmlib.synthvol import LABELS  # noqa: E402

OUT = os.path.join(ROOT, "out", "next")
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------- geometry ---
NLAT, PITCH_VOX, R_VOX, MARGIN, JITTER_FRAC = 4, 28, 14.0, 10, 0.15
VOXEL_NM = 20.0
SEEDS = [0, 1, 2, 3, 4]
TARGET_RATIOS = [1.5, 2.0]           # primary envelope ONLY (2.5x excluded)

# ------------------------------------------------------- base distribution
FRAC_WEAK = 0.20
WEAK_RANGE = (4, 6)
NORMAL_RANGE = (13, 22)
MIN_RATIO = 2.5

# ------------------------------------------------------------------ gates
MAX_BISECT_ITER = 5
P10_TOL_REL = 0.08
BASE_N_NODES_NOMINAL = NLAT ** 3


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
                   n_nodes=nm["n_nodes"], n_edges=nm["n_edges"],
                   mean_degree=nm["mean_degree"])
    else:
        out.update(neck_p10_nm=np.nan, neck_p50_nm=np.nan, n_nodes=0,
                   n_edges=0, mean_degree=np.nan)
    return out


def eval_point(centres, pairs, shape, r_base, base_widths, threshold, spacing):
    final_widths = max_clip_widths(base_widths, threshold)
    build = build_mass_conservative_structure(centres, pairs, shape, r_base,
                                              base_widths, final_widths)
    m = measure_structure(build["ni_mask"], spacing)
    return build, m


def bisect_threshold(centres, pairs, shape, base_widths, r_base,
                     target_p10_nm, spacing, t_lo, t_hi,
                     max_iter=MAX_BISECT_ITER, tol_rel=P10_TOL_REL):
    log = []
    best = None
    lo, hi = t_lo, t_hi
    for _ in range(max_iter):
        T = 0.5 * (lo + hi)
        build, m = eval_point(centres, pairs, shape, r_base, base_widths, T,
                              spacing)
        p10 = m["neck_p10_nm"]
        log.append((T, p10))
        err = abs(p10 - target_p10_nm) / target_p10_nm if np.isfinite(p10) else np.inf
        if best is None or err < best[0]:
            best = (err, T, build, m)
        if err <= tol_rel:
            break
        if not np.isfinite(p10) or p10 < target_p10_nm:
            lo = T
        else:
            hi = T
    _, T_best, build_best, m_best = best
    return T_best, build_best, m_best, log


def gate(row, target_min, base_n_nodes):
    checks = {
        "p10_ge_target": row["achieved_p10_ratio"] >= target_min,
        "phi_le_5pct": abs(row["phi_dev_pct"]) <= 5.0,
        "phi_le_2pct_target": abs(row["phi_dev_pct"]) <= 2.0,
        "cpsd_le_5pct": abs(row["cpsd_dev_pct"]) <= 5.0,
        "p50_ratio_le_1_15": row["neck_p50_ratio"] <= 1.15,
        "n_nodes_ge_95pct": row["n_nodes"] >= 0.95 * base_n_nodes,
        "P_span_intact": row["P_span"] >= 0.90,
    }
    full = (checks["p10_ge_target"] and checks["phi_le_5pct"]
           and checks["cpsd_le_5pct"] and checks["p50_ratio_le_1_15"]
           and checks["n_nodes_ge_95pct"] and checks["P_span_intact"])
    return checks, full


def main():
    print("=" * 78)
    print("FAMILY B DISORDERED PILOT")
    print("=" * 78)
    print(f"lattice {NLAT}^3 spheres (jittered, frac={JITTER_FRAC}), "
          f"R={R_VOX}vox (D={2*R_VOX}vox={2*R_VOX*VOXEL_NM:.0f}nm), "
          f"pitch={PITCH_VOX}vox, voxel={VOXEL_NM}nm")
    print(f"base neck mixture: {int(FRAC_WEAK*100)}% weak{WEAK_RANGE}vox, "
          f"{int((1-FRAC_WEAK)*100)}% normal{NORMAL_RANGE}vox, "
          f"validity floor p50/p10>={MIN_RATIO}")
    print(f"primary envelope: {TARGET_RATIOS} (2.5x excluded per amendment A)\n")

    spacing = (VOXEL_NM, VOXEL_NM, VOXEL_NM)
    t0 = time.time()
    rows = []
    base_cache = {}
    rejection_log_all = []

    print("--- geometry (shared across seeds; only base widths vary) ---")
    geom_rng = np.random.default_rng(12345)
    centres, pairs, shape = jittered_lattice_geometry(
        NLAT, PITCH_VOX, R_VOX, MARGIN, JITTER_FRAC, geom_rng)
    print(f"  {len(centres)} spheres, {len(pairs)} candidate neck pairs, "
          f"domain {shape} = {np.prod(shape)/1e6:.2f} Mvoxel\n")

    print("--- base distribution validity + base structures ---")
    for seed in SEEDS:
        widths, acc_seed, n_attempts, log = draw_valid_base_widths(
            len(pairs), seed, FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO)
        for entry in log:
            rejection_log_all.append(dict(seed=seed, **entry))
        st = base_distribution_stats(widths)
        build, m = eval_point(centres, pairs, shape, R_VOX, widths, 0.0,
                              spacing)
        base_cache[seed] = dict(widths=widths, build=build, measure=m,
                                stats=st, n_attempts=n_attempts)
        rows.append(dict(mode="base", target_ratio=1.0, seed=seed,
                         base_n_attempts=n_attempts,
                         base_p10_vox=st["p10_vox"], base_p50_vox=st["p50_vox"],
                         base_ratio=st["ratio"], intended_T_vox=np.nan,
                         iters=0, **_flat(build, m)))
        print(f"  seed={seed}: n_attempts={n_attempts} "
              f"base(p10={st['p10_vox']:.1f}vox p50={st['p50_vox']:.1f}vox "
              f"ratio={st['ratio']:.2f})  "
              f"measured(P_span={m['P_span']:.2f} n_nodes={m['n_nodes']} "
              f"mean_degree={m.get('mean_degree', float('nan')):.2f} "
              f"neck_p10={m['neck_p10_nm']:.0f}nm neck_p50={m['neck_p50_nm']:.0f}nm "
              f"ws_d={m.get('ws_d_volweighted_nm', float('nan')):.0f}nm "
              f"cPSD_d={m.get('d_cPSD_r50max_nm', float('nan')):.0f}nm)  "
              f"[{time.time()-t0:.0f}s]")

    print("\n--- lower-tail (percentile-targeted, mass-conservative) widening ---")
    for ratio in TARGET_RATIOS:
        for seed in SEEDS:
            widths = base_cache[seed]["widths"]
            base_p10 = base_cache[seed]["measure"]["neck_p10_nm"]
            target_p10 = ratio * base_p10
            t_lo = base_p10 / VOXEL_NM
            t_hi = max(NORMAL_RANGE) + 4
            T_best, build, m, log = bisect_threshold(
                centres, pairs, shape, widths, R_VOX, target_p10, spacing,
                t_lo, t_hi)
            achieved_ratio = m["neck_p10_nm"] / base_p10
            rows.append(dict(mode="lower_tail", target_ratio=ratio, seed=seed,
                             base_n_attempts=base_cache[seed]["n_attempts"],
                             base_p10_vox=base_cache[seed]["stats"]["p10_vox"],
                             base_p50_vox=base_cache[seed]["stats"]["p50_vox"],
                             base_ratio=base_cache[seed]["stats"]["ratio"],
                             intended_T_vox=T_best, iters=len(log),
                             achieved_ratio=achieved_ratio,
                             **_flat(build, m)))
            phi_dev = 100 * (build["phi_Ni_final"] - build["phi_Ni_base"]) / build["phi_Ni_base"]
            print(f"  ratio={ratio} seed={seed}: T={T_best:.2f}vox "
                  f"iters={len(log)} achieved_p10={achieved_ratio:.2f} "
                  f"neck_p50={m['neck_p50_nm']:.0f}nm phi_dev={phi_dev:+.1f}% "
                  f"n_nodes={m['n_nodes']} P_span={m['P_span']:.2f} "
                  f"[{time.time()-t0:.0f}s]")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "familyB_pilot.csv"), index=False)
    rej_df = pd.DataFrame(rejection_log_all)
    rej_df.to_csv(os.path.join(OUT, "familyB_pilot_base_validity_log.csv"),
                  index=False)
    print(f"\n[saved] {os.path.join(OUT, 'familyB_pilot.csv')}")
    print(f"[saved] {os.path.join(OUT, 'familyB_pilot_base_validity_log.csv')}"
          f"   total {time.time()-t0:.0f}s")

    analyse(df)
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


def analyse(df):
    print("\n" + "=" * 78)
    print("PER-SEED STRICT GATING (preregistration.md #0c amendment B — "
          "NO aggregate-mean criterion)")
    print("=" * 78)

    base = df[df["mode"] == "base"].set_index("seed")
    base_n_nodes = {s: base.loc[s, "n_nodes"] for s in SEEDS}
    print("Base n_nodes per seed (nominal 64):",
          {s: int(base_n_nodes[s]) for s in SEEDS})
    print("Base mean_degree per seed (target ~6-8, reported honestly):",
          {s: round(base.loc[s, "mean_degree"], 2) for s in SEEDS})

    recs = []
    for _, r in df[df["mode"] != "base"].iterrows():
        b = base.loc[r.seed]
        rec = dict(mode=r["mode"], nominal_target_ratio=r.target_ratio,
                  seed=r.seed, achieved_p10_ratio=r.achieved_ratio,
                  neck_p10_nm=r.neck_p10_nm, neck_p50_nm=r.neck_p50_nm,
                  neck_p50_ratio=r.neck_p50_nm / b.neck_p50_nm,
                  ws_dev_pct=100 * (r.ws_d_volweighted_nm - b.ws_d_volweighted_nm) / b.ws_d_volweighted_nm,
                  cpsd_dev_pct=100 * (r.d_cPSD_r50max_nm - b.d_cPSD_r50max_nm) / b.d_cPSD_r50max_nm,
                  phi_dev_pct=100 * (r.phi_Ni_final - b.phi_Ni_final) / b.phi_Ni_final,
                  primary_size_dev_pct=100 * (r.r_final_vox - b.r_final_vox) / b.r_final_vox,
                  n_nodes=r.n_nodes, base_n_nodes=int(b.n_nodes),
                  P_span=r.P_span, n_clusters=r.n_clusters,
                  voxels_added=r.voxels_added_by_necks,
                  voxels_removed=r.voxels_removed_by_shrink,
                  voxels_net_residual=r.voxels_net_residual)
        recs.append(rec)
    dev = pd.DataFrame(recs)

    all_rows = []
    for tr in TARGET_RATIOS:
        sub = dev[np.isclose(dev.nominal_target_ratio, tr)]
        print(f"\n--- nominal target ratio = {tr} "
              f"(achieved values are the analysis variable) ---")
        n_pass = 0
        for _, r in sub.iterrows():
            checks, full = gate(r, tr, r.base_n_nodes)
            n_pass += int(full)
            failed = [k for k, v in checks.items() if not v]
            print(f"  seed={r.seed}  achieved_p10={r.achieved_p10_ratio:.2f} "
                  f"(nominal {tr})  phi_dev={r.phi_dev_pct:+.2f}%  "
                  f"cpsd_dev={r.cpsd_dev_pct:+.2f}%  "
                  f"p50_ratio={r.neck_p50_ratio:.2f}  "
                  f"n_nodes={r.n_nodes}/{r.base_n_nodes}  "
                  f"P_span={r.P_span:.3f}  n_clusters={r.n_clusters}  "
                  f"-> {'PASS' if full else 'fail'}"
                  + (f"  (failed: {', '.join(failed)})" if failed else ""))
            all_rows.append(dict(nominal_target_ratio=tr, seed=r.seed,
                                 achieved_p10_ratio=r.achieved_p10_ratio,
                                 full_pass=full,
                                 failed_checks=",".join(failed)))
        feasible = n_pass >= 4
        print(f"  ratio={tr}: {n_pass}/5 seeds pass -> "
              f"{'FEASIBLE' if feasible else 'NOT FEASIBLE'}")

    dev.to_csv(os.path.join(OUT, "familyB_pilot_deviations.csv"), index=False)
    gate_df = pd.DataFrame(all_rows)
    gate_df.to_csv(os.path.join(OUT, "familyB_pilot_gating.csv"), index=False)
    print(f"\n[saved] {os.path.join(OUT, 'familyB_pilot_deviations.csv')}")
    print(f"[saved] {os.path.join(OUT, 'familyB_pilot_gating.csv')}")

    # -------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    colors = {1.5: "C0", 2.0: "C1"}
    for tr in TARGET_RATIOS:
        sub = dev[np.isclose(dev.nominal_target_ratio, tr)]
        c = colors[tr]
        axes[0].scatter(sub.achieved_p10_ratio, sub.phi_dev_pct, color=c,
                        label=f"nominal {tr}x")
        axes[1].scatter(sub.achieved_p10_ratio, sub.cpsd_dev_pct, color=c,
                        label=f"nominal {tr}x")
        axes[2].scatter(sub.achieved_p10_ratio, sub.neck_p50_ratio, color=c,
                        label=f"nominal {tr}x")
    axes[0].axhspan(-5, 5, color="green", alpha=0.08)
    axes[0].axhspan(-2, 2, color="green", alpha=0.15)
    axes[0].set_xlabel("achieved neck p10 ratio")
    axes[0].set_ylabel("Phi_Ni deviation (%)")
    axes[0].set_title("Mass conservation")
    axes[1].axhspan(-5, 5, color="green", alpha=0.1)
    axes[1].set_xlabel("achieved neck p10 ratio")
    axes[1].set_ylabel("c-PSD size deviation (%)")
    axes[1].set_title("Primary size gate (c-PSD)")
    axes[2].axhline(1.15, color="crimson", ls="--", lw=1, label="1.15 ceiling")
    axes[2].set_xlabel("achieved neck p10 ratio")
    axes[2].set_ylabel("measured neck p50 ratio")
    axes[2].set_title("Tail-selectivity")
    for a in axes:
        a.legend(fontsize=8, frameon=False)
        a.grid(alpha=0.25)
    fig.suptitle("Family B disordered pilot — per-seed results "
                 "(achieved p10 ratio is the analysis variable)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "familyB_pilot.png"), dpi=145)
    print(f"[saved] {os.path.join(OUT, 'familyB_pilot.png')}")


if __name__ == "__main__":
    sys.exit(main())
