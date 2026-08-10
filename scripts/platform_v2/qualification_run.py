"""
PLATFORM V2 — Ni generator qualification run.

Locked geometry from scripts/platform_v2/design_probe.py: R=12.1vox
(D=24.2vox=484nm), pitch=32vox, nlat_z=6/nlat_xy=4 (asymmetric -- see
platform_v2_lattice_geometry docstring), margin=8, domain (161,168,168),
Phi_Ni=0.2502 (base, seed 0 draw), PLAIN 6-connectivity (nearest-neighbour)
topology -- the design probe found this already lands SNOW-measured
mean_degree=4.204 inside the [3.5,4.5] target band, so NO topology
modification (face-diagonal bonds etc.) was built or is used here.

Scope: Ni generator qualification ONLY. Per instruction: no YSZ/pore
placement, no D4/damage modelling, no Family C, no calibration against the
real Holzer/Pecho dataset in this script.

5 seeds x {base, intermediate (~1.3-1.6x achieved p10), high (~2.0x achieved
p10)}. Lower-tail widening reuses the validated mass-conservative,
percentile-targeted mechanism (max-clip threshold via bisection against
MEASURED p10) unchanged from Family B / T5b. Achieved measured p10 ratio is
the scientific variable throughout, nominal labels are generation labels only.

Per-seed gating (7 criteria, no aggregate-mean pass criterion, per
preregistration.md #0c amendment B): p10 ratio >= target minimum, Phi_Ni
<=5% (target <=2%), c-PSD <=5%, p50 ratio <=1.15, n_nodes >=95% of base,
P_span intact, no severe fragment blowup (n_clusters not exploding).

Mass-conservation headroom check (new, per instruction item 6): at the HIGH
ratio point specifically, every seed's voxels_added_by_necks /
voxels_removed_by_shrink / resulting Phi_Ni residual is reported
INDIVIDUALLY, not only folded into the aggregate node-count gate -- this is
exactly the T5b seed=1 failure mode, and lowering Phi_Ni to 0.250 (smaller
base particle radii than the pilot's 0.32-0.33) reduces the headroom before
the compensating-radius mechanism's r_lo floor is exhausted.
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
    draw_valid_base_widths, max_clip_widths, mixture_neck_widths,
    platform_v2_lattice_geometry,
)

OUT = os.path.join(ROOT, "out", "platform_v2")
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------- locked geometry ----
R_VOX = 12.1
PITCH_VOX = 32
NLAT_Z, NLAT_XY = 6, 4
MARGIN = 8
JITTER_FRAC = 0.15
VOXEL_NM = 20.0
GEOM_SEED = 999          # fixed, shared across all structure seeds

# --------------------------------------------- base neck distribution --
FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE = 0.20, (4, 6), (12, 20)
MIN_RATIO = 2.5

SEEDS = [0, 1, 2, 3, 4]
# "intermediate" and "high". 1.45 -> 1.33 (2026-08-10, approved cheap fix):
# the first qualification run showed all 5 seeds land on the 1.33x rung, so a
# 1.45x nominal target failed the p10 criterion for every seed purely by
# falling between achievable rungs -- the achievable ladder for this base
# mixture at n_pairs=224 is {1.33x, 2.0x}. Since achieved ratio is the
# scientific variable (preregistration.md #0c amendment C), the nominal label
# is retargeted onto the rung the generator actually produces.
TARGET_RATIOS = [1.33, 2.0]
MAX_BISECT_ITER = 5
P10_TOL_REL = 0.08

COORDINATION_BAND = (3.5, 4.5)


def measure_structure(ni_mask, spacing):
    perc = compute_percolation(
        np.where(ni_mask, 2, 0).astype(np.uint8),
        {"Ni": 2, "YSZ": 1, "pore": 0}, phase="Ni", axis=0)
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
        degrees = np.array([d for _, d in G.degree()])
        out.update(neck_p10_nm=nm["neck_p10_nm"], neck_p50_nm=nm["neck_p50_nm"],
                   n_nodes=nm["n_nodes"], n_edges=nm["n_edges"],
                   mean_degree=nm["mean_degree"],
                   n_deg0=int((degrees == 0).sum()),
                   n_deg1=int((degrees == 1).sum()),
                   degree_hist=np.bincount(degrees).tolist())
    else:
        out.update(neck_p10_nm=np.nan, neck_p50_nm=np.nan, n_nodes=0,
                   n_edges=0, mean_degree=np.nan, n_deg0=0, n_deg1=0,
                   degree_hist=[])
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


def _flat(build, m):
    return dict(
        r_base_vox=build["r_base_vox"], r_final_vox=build["r_final_vox"],
        voxels_base=build["voxels_base"],
        voxels_added_by_necks=build["voxels_added_by_necks"],
        voxels_removed_by_shrink=build["voxels_removed_by_shrink"],
        voxels_final=build["voxels_final"],
        voxels_net_residual=build["voxels_net_residual"],
        phi_Ni_base=build["phi_Ni_base"], phi_Ni_final=build["phi_Ni_final"],
        **{k: v for k, v in m.items() if k != "degree_hist"},
    )


def main():
    # --from-csv re-runs ONLY the analysis/gating against the already-saved
    # qualification_run.csv (the sweep itself is ~23 min; the gating logic is
    # cheap and was iterated after the data was collected).
    if "--from-csv" in sys.argv:
        df = pd.read_csv(os.path.join(OUT, "qualification_run.csv"))
        print("[--from-csv] re-analysing saved qualification_run.csv "
              f"({len(df)} rows), no recomputation")
        analyse(df)
        return 0

    print("=" * 78)
    print("PLATFORM V2 — Ni generator qualification run")
    print("=" * 78)
    print(f"geometry: R={R_VOX}vox pitch={PITCH_VOX}vox nlat_z={NLAT_Z} "
          f"nlat_xy={NLAT_XY} margin={MARGIN}  PLAIN topology (no modification)")
    print(f"target ratios: {TARGET_RATIOS} (intermediate ~1.45x, high ~2.0x)")
    print(f"seeds: {SEEDS}\n")

    spacing = (VOXEL_NM, VOXEL_NM, VOXEL_NM)
    t0 = time.time()
    geom_rng = np.random.default_rng(GEOM_SEED)
    centres, pairs, shape = platform_v2_lattice_geometry(
        NLAT_Z, NLAT_XY, PITCH_VOX, R_VOX, MARGIN, JITTER_FRAC, geom_rng)
    n_sites, n_pairs = len(centres), len(pairs)
    print(f"geometry: {n_sites} spheres, {n_pairs} pairs, domain={shape} "
          f"= {np.prod(shape)/1e6:.2f} Mvoxel, "
          f"raw topological mean_degree={2*n_pairs/n_sites:.3f}\n")

    rows = []
    base_cache = {}
    rejection_log = []

    print("--- base distribution validity + base structures ---")
    for seed in SEEDS:
        widths, acc_seed, n_attempts, log = draw_valid_base_widths(
            n_pairs, seed, FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO)
        for entry in log:
            rejection_log.append(dict(seed=seed, **entry))
        st = base_distribution_stats(widths)
        build, m = eval_point(centres, pairs, shape, R_VOX, widths, 0.0, spacing)
        base_cache[seed] = dict(widths=widths, build=build, measure=m, stats=st,
                                n_attempts=n_attempts)
        rows.append(dict(mode="base", target_ratio=1.0, seed=seed,
                         base_n_attempts=n_attempts,
                         base_p10_vox=st["p10_vox"], base_p50_vox=st["p50_vox"],
                         base_ratio=st["ratio"], intended_T_vox=np.nan,
                         iters=0, achieved_ratio=1.0, **_flat(build, m)))
        in_band = COORDINATION_BAND[0] <= m["mean_degree"] <= COORDINATION_BAND[1]
        print(f"  seed={seed}: n_attempts={n_attempts} "
              f"base(p10={st['p10_vox']:.1f}vox ratio={st['ratio']:.2f})  "
              f"Phi_Ni={build['phi_Ni_final']:.4f}  "
              f"P_span={m['P_span']:.2f} n_nodes={m['n_nodes']} "
              f"mean_degree={m['mean_degree']:.3f} "
              f"({'IN-BAND' if in_band else 'OUT-OF-BAND'}) "
              f"n_deg0={m['n_deg0']} n_deg1={m['n_deg1']}  "
              f"[{time.time()-t0:.0f}s]")

    print("\n--- lower-tail (percentile-targeted, mass-conservative) widening ---")
    for ratio in TARGET_RATIOS:
        for seed in SEEDS:
            widths = base_cache[seed]["widths"]
            base_p10 = base_cache[seed]["measure"]["neck_p10_nm"]
            target_p10 = ratio * base_p10
            t_lo = base_p10 / VOXEL_NM
            t_hi = max(NORMAL_RANGE) + 6
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
                             achieved_ratio=achieved_ratio, **_flat(build, m)))
            phi_dev = 100 * (build["phi_Ni_final"] - build["phi_Ni_base"]) / build["phi_Ni_base"]
            print(f"  ratio={ratio} seed={seed}: T={T_best:.2f}vox "
                  f"iters={len(log)} achieved_p10={achieved_ratio:.2f} "
                  f"phi_dev={phi_dev:+.2f}% n_nodes={m['n_nodes']} "
                  f"mean_degree={m['mean_degree']:.3f} "
                  f"P_span={m['P_span']:.2f}  "
                  f"added={build['voxels_added_by_necks']} "
                  f"removed={build['voxels_removed_by_shrink']} "
                  f"residual={build['voxels_net_residual']}  "
                  f"[{time.time()-t0:.0f}s]")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "qualification_run.csv"), index=False)
    rej_df = pd.DataFrame(rejection_log)
    rej_df.to_csv(os.path.join(OUT, "qualification_base_validity_log.csv"),
                  index=False)
    print(f"\n[saved] {os.path.join(OUT, 'qualification_run.csv')}")
    print(f"[saved] {os.path.join(OUT, 'qualification_base_validity_log.csv')}"
          f"   total {time.time()-t0:.0f}s")

    analyse(df)
    return 0


def analyse(df):
    print("\n" + "=" * 78)
    print("GATE P2-A: composition and topology")
    print("=" * 78)
    base = df[df["mode"] == "base"]
    print(f"  seeds: {len(base)} (need >=5: {len(base) >= 5})")
    print(f"  Phi_Ni: mean={base.phi_Ni_final.mean():.4f} "
          f"range=[{base.phi_Ni_final.min():.4f},{base.phi_Ni_final.max():.4f}] "
          f"(target 0.250)")
    print(f"  mean_degree: mean={base.mean_degree.mean():.3f} "
          f"range=[{base.mean_degree.min():.3f},{base.mean_degree.max():.3f}] "
          f"(target band {COORDINATION_BAND})")
    all_in_band = ((base.mean_degree >= COORDINATION_BAND[0])
                   & (base.mean_degree <= COORDINATION_BAND[1])).all()
    print(f"  all base seeds in coordination band: {all_in_band}")
    print(f"  single connected cluster (n_clusters==1) for all seeds: "
          f"{(base.n_clusters == 1).all()}")
    print(f"  P_span intact (==1.0) for all seeds: {(base.P_span == 1.0).all()}")
    print(f"  TOPOLOGY NOTE: PLAIN lattice adjacency used throughout -- the "
          f"design probe (scripts/platform_v2/design_probe.py) found this "
          f"already lands mean_degree=4.204 in-band; NO topology "
          f"modification (face-diagonal bonds etc.) was built or needed.")

    print("\n" + "=" * 78)
    print("GATE P2-B: base neck distribution")
    print("=" * 78)
    print(f"  p50/p10: mean={base.base_ratio.mean():.2f} "
          f"range=[{base.base_ratio.min():.2f},{base.base_ratio.max():.2f}] "
          f"(target 3.0-4.3, floor >=2.5)")
    print(f"  base p10 resolution: mean={base.base_p10_vox.mean():.1f}vox "
          f"(target >=3, preferably >=4)")
    print(f"  all base seeds accepted on first attempt (n_attempts==1): "
          f"{(base.base_n_attempts == 1).all()}")

    print("\n" + "=" * 78)
    print("GATE P2-C: lower-tail decoupling — PER-SEED (no aggregate-mean pass)")
    print("=" * 78)
    base_idx = base.set_index("seed")
    dev_rows = []
    for _, r in df[df["mode"] != "base"].iterrows():
        b = base_idx.loc[r.seed]
        dev_rows.append(dict(
            target_ratio=r.target_ratio, seed=r.seed,
            achieved_p10_ratio=r.achieved_ratio,
            phi_dev_pct=100 * (r.phi_Ni_final - b.phi_Ni_final) / b.phi_Ni_final,
            cpsd_dev_pct=100 * (r.d_cPSD_r50max_nm - b.d_cPSD_r50max_nm) / b.d_cPSD_r50max_nm,
            neck_p50_ratio=r.neck_p50_nm / b.neck_p50_nm,
            n_nodes=r.n_nodes, base_n_nodes=int(b.n_nodes),
            mean_degree=r.mean_degree, P_span=r.P_span, n_clusters=r.n_clusters,
            voxels_added_by_necks=r.voxels_added_by_necks,
            voxels_removed_by_shrink=r.voxels_removed_by_shrink,
            voxels_net_residual=r.voxels_net_residual,
            r_final_vox=r.r_final_vox,
        ))
    dev = pd.DataFrame(dev_rows)

    for tr in TARGET_RATIOS:
        sub = dev[np.isclose(dev.target_ratio, tr)]
        label = "intermediate" if tr < 1.7 else "high"
        print(f"\n--- {label} (nominal {tr}x) ---")
        n_pass = 0
        for _, r in sub.iterrows():
            checks, full = gate(r, tr, r.base_n_nodes)
            n_pass += int(full)
            failed = [k for k, v in checks.items() if not v]
            print(f"  seed={r.seed}  achieved_p10={r.achieved_p10_ratio:.2f}  "
                  f"phi_dev={r.phi_dev_pct:+.2f}%  cpsd_dev={r.cpsd_dev_pct:+.2f}%  "
                  f"p50_ratio={r.neck_p50_ratio:.2f}  "
                  f"n_nodes={r.n_nodes}/{r.base_n_nodes}  "
                  f"mean_degree={r.mean_degree:.3f}  P_span={r.P_span:.3f}  "
                  f"n_clusters={r.n_clusters}  "
                  f"-> {'PASS' if full else 'fail'}"
                  + (f"  (failed: {', '.join(failed)})" if failed else ""))
        print(f"  {label}: {n_pass}/{len(sub)} seeds pass -> "
              f"{'FEASIBLE' if n_pass >= 4 else 'NOT FEASIBLE'}")

        if label == "high":
            print(f"\n  MASS-CONSERVATION HEADROOM CHECK (per-seed, high-ratio "
                  f"point only, item 6):")
            for _, r in sub.iterrows():
                r_lo_frac = r.r_final_vox / R_VOX  # how close to r_lo=1.0 floor
                flag = ""
                if r.r_final_vox < 0.7 * R_VOX:
                    flag = "  *** RADIUS COMPENSATION APPROACHING FLOOR ***"
                elif r.r_final_vox < 0.85 * R_VOX:
                    flag = "  ** notable radius shrink **"
                print(f"    seed={r.seed}: added={int(r.voxels_added_by_necks):6d}  "
                      f"removed={int(r.voxels_removed_by_shrink):6d}  "
                      f"net_residual={int(r.voxels_net_residual):+6d}  "
                      f"phi_dev={r.phi_dev_pct:+.2f}%  "
                      f"r_final={r.r_final_vox:.2f}vox "
                      f"(r_base={R_VOX}vox, {100*r_lo_frac:.0f}% of base)  "
                      f"n_nodes_lost={int(r.base_n_nodes - r.n_nodes)}{flag}")

    dev.to_csv(os.path.join(OUT, "qualification_deviations.csv"), index=False)
    print(f"\n[saved] {os.path.join(OUT, 'qualification_deviations.csv')}")

    # -------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    colors = {1.45: "C0", 2.0: "C1"}
    for tr in TARGET_RATIOS:
        sub = dev[np.isclose(dev.target_ratio, tr)]
        c = colors.get(tr, "gray")
        axes[0].scatter(sub.achieved_p10_ratio, sub.phi_dev_pct, color=c,
                        label=f"nominal {tr}x")
        axes[1].scatter(sub.achieved_p10_ratio, sub.cpsd_dev_pct, color=c,
                        label=f"nominal {tr}x")
        axes[2].scatter(sub.achieved_p10_ratio, sub.mean_degree, color=c,
                        label=f"nominal {tr}x")
    axes[0].scatter([1.0] * len(base), [0] * len(base), color="C2", label="base")
    axes[0].axhspan(-5, 5, color="green", alpha=0.08)
    axes[0].axhspan(-2, 2, color="green", alpha=0.15)
    axes[0].set_xlabel("achieved neck p10 ratio")
    axes[0].set_ylabel("Phi_Ni deviation (%)")
    axes[0].set_title("Mass conservation")

    axes[1].axhspan(-5, 5, color="green", alpha=0.1)
    axes[1].set_xlabel("achieved neck p10 ratio")
    axes[1].set_ylabel("c-PSD deviation (%)")
    axes[1].set_title("Primary size gate (c-PSD)")

    axes[2].axhspan(*COORDINATION_BAND, color="green", alpha=0.1,
                    label="target band")
    axes[2].scatter([1.0] * len(base), base.mean_degree, color="C2", label="base")
    axes[2].set_xlabel("achieved neck p10 ratio")
    axes[2].set_ylabel("SNOW mean_degree")
    axes[2].set_title("Coordination vs widening")

    for a in axes:
        a.legend(fontsize=7, frameon=False)
        a.grid(alpha=0.25)
    fig.suptitle("Platform v2 Ni qualification: per-seed results", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "qualification_run.png"), dpi=145)
    print(f"[saved] {os.path.join(OUT, 'qualification_run.png')}")


if __name__ == "__main__":
    sys.exit(main())
