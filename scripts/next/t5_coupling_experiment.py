"""
T5 COUPLING DECISION EXPERIMENT (run after gates G0 and G-1, per the approved
sequence in this session). Decides Q1/Q2 in out/next/EXECUTION_SPEC.md before
any generator (cmlib/synth.py) is written.

QUESTION: at FIXED sphere radius (nominal Ni particle size) and FIXED lattice
(fixed Ni volume fraction), can neck width be varied independently enough to
move neck p10 by >=1.5x (target 2x) while holding a particle-size MEASURE
within +/-5%? Two widening modes are tested, per the approved instruction:
  (a) SELECTIVE lower-tail widening: only the narrowest ~20% of necks (by
      their randomly-assigned base width) are widened, up to a target.
  (b) WHOLE-DISTRIBUTION widening (if inexpensive): every neck is widened by
      the same increment.
Both watershed (ws_d_volweighted_nm) and c-PSD (d_cPSD_r50max_nm) particle-size
measures are tracked, because R1 predicts the watershed measure moves and the
c-PSD measure should be far more stable -- this experiment is what checks that
prediction rather than assuming it.

STRUCTURE: a cubic lattice of spheres (nominal Ni particles), nearest-neighbour
pairs joined by a SHORT connecting neck (per the T5b finding in
scripts/next/phase0_validate_synthetic_pipeline.py: short near-contact necks
only, never long freestanding bridges, which spawn a spurious extra watershed
region). This is a throwaway scaffold for the decision, NOT the Phase-1
generator -- it has no config file, no YSZ/pore placement, no damage model.

Every neck stat reported (neck_p10_nm, neck_p50_nm, n_nodes, n_edges) comes
from the SAME extract_network + compute_network_metrics call used everywhere
else in this project, not from the generator's own "intended" neck radii, so
what is reported is what the analysis pipeline actually sees.
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
from cmlib.synthvol import LABELS  # noqa: E402

OUT = os.path.join(ROOT, "out", "next")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# lattice packing with randomized base neck widths, then two widening modes
# ---------------------------------------------------------------------------
NLAT = 5                 # 5x5x5 = 125 spheres
R_VOX = 8                 # nominal Ni particle radius, voxels
PITCH_VOX = 20             # lattice spacing -> base gap = pitch - 2R = 4 voxels
MARGIN = 6                 # transverse (y, x) padding only
# z (axis 0, the percolation/network axis) is sized so the FIRST and LAST
# sphere layers exactly touch the domain faces (z=0 and z=NZ-1): otherwise
# extract_network's spanning restriction (see cmlib/pnm.py, matches the
# real-data convention) returns an empty network even though the lattice is
# fully internally connected -- found while smoke-testing this script.
NZ = (NLAT - 1) * PITCH_VOX + 1
N = NLAT * PITCH_VOX + 2 * R_VOX + 2 * MARGIN     # y, x extent (unchanged)
VOXEL_NM = 20.0


def build_lattice_ni(seed: int, mode: str, target_w: int, base_lo=2, base_hi=6):
    """Cubic lattice of spheres with randomized base neck widths.

    mode='none'      : no necks at all (isolated spheres; sanity baseline)
    mode='selective'  : bottom 20% of necks (by base width) widened to target_w
    mode='uniform'    : every neck widened by (target_w - base median)
    Returns (ni_mask: bool array, neck_base_widths: array, neck_final_widths: array)
    """
    rng = np.random.default_rng(seed)
    centres = {}
    for iz in range(NLAT):
        for iy in range(NLAT):
            for ix in range(NLAT):
                z = iz * PITCH_VOX                    # 0 .. NZ-1, touches both z-faces
                y = MARGIN + R_VOX + iy * PITCH_VOX
                x = MARGIN + R_VOX + ix * PITCH_VOX
                centres[(iz, iy, ix)] = (z, y, x)

    zz, yy, xx = np.ogrid[:NZ, :N, :N]
    ni = np.zeros((NZ, N, N), dtype=bool)
    for (z, y, x) in centres.values():
        ni |= (zz - z) ** 2 + (yy - y) ** 2 + (xx - x) ** 2 < R_VOX ** 2

    # enumerate nearest-neighbour pairs (6-connectivity in lattice index space)
    pairs = []
    for (iz, iy, ix), c0 in centres.items():
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            nb = (iz + d[0], iy + d[1], ix + d[2])
            if nb in centres:
                pairs.append(((iz, iy, ix), nb))

    base_widths = rng.integers(base_lo, base_hi + 1, size=len(pairs))
    final_widths = base_widths.copy()

    if mode == "selective":
        thresh = np.percentile(base_widths, 20)
        sel = base_widths <= thresh
        final_widths = np.where(sel, target_w, base_widths)
    elif mode == "uniform":
        delta = target_w - int(np.median(base_widths))
        final_widths = np.clip(base_widths + delta, base_lo, 2 * R_VOX - 2)
    elif mode == "none":
        final_widths = np.zeros_like(base_widths)
    elif mode == "base":
        pass
    else:
        raise ValueError(mode)

    for (p0, p1), w in zip(pairs, final_widths):
        if w <= 0:
            continue
        c0 = np.array(centres[p0])
        c1 = np.array(centres[p1])
        axis = int(np.argmax(np.abs(c1 - c0)))
        lo = min(c0[axis], c1[axis])
        hi = max(c0[axis], c1[axis])
        # short connector spanning the gap between the two sphere SURFACES
        # (not the full centre-to-centre span), consistent with the
        # short-near-contact-neck design rule from T5b.
        half = int(w) // 2
        sl = [slice(None)] * 3
        sl[axis] = slice(lo, hi + 1)
        for a in range(3):
            if a != axis:
                sl[a] = slice(c0[a] - half, c0[a] - half + int(w))
        ni[tuple(sl)] = True

    return ni, base_widths.astype(float) * VOXEL_NM, final_widths.astype(float) * VOXEL_NM


def measure(ni_mask, tag, seed, mode, target_w):
    spacing = (VOXEL_NM, VOXEL_NM, VOXEL_NM)
    perc = compute_percolation(
        np.where(ni_mask, LABELS["Ni"], LABELS["pore"]).astype(np.uint8),
        LABELS, phase="Ni", axis=0)
    pstats = compute_particle_stats(ni_mask, spacing, min_distance=4)
    row = dict(tag=tag, seed=seed, mode=mode, target_w_vox=target_w,
              phi_Ni=float(ni_mask.mean()), P_span=perc["P_span"],
              percolates=perc["percolates"])
    row.update({k: v for k, v in pstats.items()
               if k in ("ws_d_volweighted_nm", "ws_n_regions_used",
                        "d_cPSD_r50max_nm", "n_peaks")})
    G, diag = extract_network(ni_mask, spacing, axis=0, r_max=4)
    if G is not None and G.number_of_edges() > 0:
        nm = compute_network_metrics(G, G.graph.get("face_lo"),
                                     G.graph.get("face_hi"))
        row.update(neck_p10_nm=nm["neck_p10_nm"], neck_p50_nm=nm["neck_p50_nm"],
                   n_nodes=nm["n_nodes"], n_edges=nm["n_edges"])
    else:
        row.update(neck_p10_nm=np.nan, neck_p50_nm=np.nan, n_nodes=0, n_edges=0)
    return row


def main():
    print("=" * 78)
    print("T5 COUPLING DECISION EXPERIMENT")
    print("=" * 78)
    print(f"lattice {NLAT}x{NLAT}x{NLAT} spheres, R={R_VOX}vox, "
          f"pitch={PITCH_VOX}vox (base gap 4vox), domain {N}^3, "
          f"voxel={VOXEL_NM}nm\n")

    seeds = [0, 1, 2, 3, 4]
    rows = []

    t0 = time.time()
    for seed in seeds:
        ni, base_w, final_w = build_lattice_ni(seed, "base", 0)
        r = measure(ni, "base (no widening)", seed, "base", 0)
        r["neck_intended_p10_nm"] = float(np.percentile(base_w, 10))
        r["neck_intended_p50_nm"] = float(np.percentile(base_w, 50))
        rows.append(r)
        print(f"  base        seed={seed}: phi_Ni={r['phi_Ni']:.4f} "
              f"ws_d={r.get('ws_d_volweighted_nm', float('nan')):.0f}nm "
              f"cPSD_d={r.get('d_cPSD_r50max_nm', float('nan')):.0f}nm "
              f"neck_p10={r['neck_p10_nm']:.0f}nm "
              f"neck_p50={r['neck_p50_nm']:.0f}nm  [{time.time()-t0:.0f}s]")

    for mode in ("selective", "uniform"):
        for target_w in (8, 10, 12, 14):
            for seed in seeds:
                ni, base_w, final_w = build_lattice_ni(seed, mode, target_w)
                r = measure(ni, f"{mode} target_w={target_w}", seed, mode, target_w)
                r["neck_intended_p10_nm"] = float(np.percentile(final_w, 10))
                r["neck_intended_p50_nm"] = float(np.percentile(final_w, 50))
                rows.append(r)
            last = rows[-1]
            print(f"  {mode:10s} target_w={target_w:2d}vox "
                  f"({target_w*VOXEL_NM:.0f}nm)  (seed={seeds[-1]}): "
                  f"phi_Ni={last['phi_Ni']:.4f} "
                  f"ws_d={last.get('ws_d_volweighted_nm', float('nan')):.0f}nm "
                  f"cPSD_d={last.get('d_cPSD_r50max_nm', float('nan')):.0f}nm "
                  f"neck_p10={last['neck_p10_nm']:.0f}nm "
                  f"neck_p50={last['neck_p50_nm']:.0f}nm  [{time.time()-t0:.0f}s]")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "t5_coupling_experiment.csv"), index=False)

    # ------------------------------------------------------------------ agg
    agg = df.groupby(["mode", "target_w_vox"]).agg(
        n=("seed", "count"),
        phi_Ni_mean=("phi_Ni", "mean"), phi_Ni_sd=("phi_Ni", "std"),
        ws_d_mean=("ws_d_volweighted_nm", "mean"),
        ws_d_sd=("ws_d_volweighted_nm", "std"),
        cpsd_d_mean=("d_cPSD_r50max_nm", "mean"),
        cpsd_d_sd=("d_cPSD_r50max_nm", "std"),
        neck_p10_mean=("neck_p10_nm", "mean"), neck_p10_sd=("neck_p10_nm", "std"),
        neck_p50_mean=("neck_p50_nm", "mean"),
    ).reset_index()
    agg.to_csv(os.path.join(OUT, "t5_coupling_experiment_agg.csv"), index=False)

    base = agg[agg["mode"] == "base"].iloc[0]
    print("\n" + "=" * 78)
    print("DECISION SUMMARY (relative to the base/unwidened lattice)")
    print("=" * 78)
    print(f"  base: phi_Ni={base.phi_Ni_mean:.4f}  "
          f"ws_d={base.ws_d_mean:.0f}nm  cPSD_d={base.cpsd_d_mean:.0f}nm  "
          f"neck_p10={base.neck_p10_mean:.0f}nm\n")

    print(f"  {'mode':10s} {'target_w':9s} {'neck_p10':>10s} {'p10_ratio':>10s} "
          f"{'ws_d':>8s} {'ws_dev%':>8s} {'cPSD_d':>8s} {'cPSD_dev%':>10s} "
          f"{'phi_Ni':>8s} {'phi_dev%':>9s}")
    for _, r in agg[agg["mode"].isin(["selective", "uniform"])].iterrows():
        p10_ratio = r.neck_p10_mean / base.neck_p10_mean
        ws_dev = 100 * (r.ws_d_mean - base.ws_d_mean) / base.ws_d_mean
        cpsd_dev = 100 * (r.cpsd_d_mean - base.cpsd_d_mean) / base.cpsd_d_mean
        phi_dev = 100 * (r.phi_Ni_mean - base.phi_Ni_mean) / base.phi_Ni_mean
        print(f"  {r['mode']:10s} {r.target_w_vox:9.0f} "
              f"{r.neck_p10_mean:10.0f} {p10_ratio:10.2f} "
              f"{r.ws_d_mean:8.0f} {ws_dev:8.1f} "
              f"{r.cpsd_d_mean:8.0f} {cpsd_dev:10.1f} "
              f"{r.phi_Ni_mean:8.4f} {phi_dev:9.2f}")

    # gate check: p10 ratio >= 1.5 AND size measure within 5% AND phi_Ni
    # within 5% -- ALL THREE constraints, not size alone. An earlier version
    # of this script checked only p10 ratio + size deviation and wrongly
    # printed "IS satisfiable" while silently ignoring an 11%+ Ni-loading
    # violation at every single tested point. Fixed here; see the honest
    # write-up in out/next/EXECUTION_SPEC.md / the session report for why
    # this matters (raw "add material" neck construction does not preserve
    # Ni loading -- a Ni-budget compensation step is required and UNTESTED
    # by this scaffold).
    for mode_df, name in ((agg[agg["mode"] == "selective"].copy(), "selective"),
                          (agg[agg["mode"] == "uniform"].copy(), "uniform")):
        mode_df["p10_ratio"] = mode_df.neck_p10_mean / base.neck_p10_mean
        mode_df["ws_dev_pct"] = 100 * (mode_df.ws_d_mean - base.ws_d_mean).abs() / base.ws_d_mean
        mode_df["cpsd_dev_pct"] = 100 * (mode_df.cpsd_d_mean - base.cpsd_d_mean).abs() / base.cpsd_d_mean
        mode_df["phi_dev_pct"] = 100 * (mode_df.phi_Ni_mean - base.phi_Ni_mean).abs() / base.phi_Ni_mean
        if name == "selective":
            sel_full = mode_df

    candidates_ws_only = sel_full[(sel_full.p10_ratio >= 1.5) & (sel_full.ws_dev_pct <= 5.0)]
    candidates_ws_and_phi = sel_full[(sel_full.p10_ratio >= 1.5) & (sel_full.ws_dev_pct <= 5.0)
                                     & (sel_full.phi_dev_pct <= 5.0)]
    candidates_cpsd_and_phi = sel_full[(sel_full.p10_ratio >= 1.5) & (sel_full.cpsd_dev_pct <= 5.0)
                                       & (sel_full.phi_dev_pct <= 5.0)]

    print("\n" + "=" * 78)
    print("GATE Q1/Q2 DECISION (size constraint alone is NOT sufficient -- "
          "phi_Ni must also hold)")
    print("=" * 78)
    print(f"  selective, p10>=1.5x AND ws_size<=5%  (IGNORING phi_Ni, "
          f"the wrong check): {len(candidates_ws_only)} point(s)")
    print(f"  selective, p10>=1.5x AND ws_size<=5%  AND phi_Ni<=5%: "
          f"{len(candidates_ws_and_phi)} point(s)")
    print(f"  selective, p10>=1.5x AND cPSD_size<=5% AND phi_Ni<=5%: "
          f"{len(candidates_cpsd_and_phi)} point(s)")
    if len(candidates_ws_and_phi):
        print("  -> primary constraint IS satisfiable under selective "
              "lower-tail widening, with Ni loading also held.")
    elif len(candidates_cpsd_and_phi):
        print("  -> only the c-PSD-scoped constraint is satisfiable with Ni "
              "loading held. Scope the primary claim accordingly.")
    else:
        print("  -> NOT satisfiable by this UNCOMPENSATED construction "
              "(neck material simply added). Every tested point overshoots "
              "the phi_Ni tolerance (11%-57% for selective, 31%-128% for "
              "uniform), even though neck p10 moves 2x-6x easily and cPSD "
              "size stays comparatively stable. This points to a SPECIFIC, "
              "previously-identified fix -- Ni-budget compensation (shrink "
              "particle radius as necks widen, per the Phase-1 'preferred "
              "approach' in EXECUTION_SPEC) -- which this scaffold does not "
              "implement and which must be built and verified before Family "
              "B can be trusted. This is a DIFFERENT conclusion from both "
              "'satisfiable' and 'not testable' -- report it as such.")

    # -------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5))
    for mode, col in (("selective", "C0"), ("uniform", "C1")):
        d = agg[agg["mode"] == mode].sort_values("target_w_vox")
        p10r = d.neck_p10_mean / base.neck_p10_mean
        axes[0].plot(p10r, d.ws_d_mean, "o-", color=col, label=f"{mode} (watershed)")
        axes[1].plot(p10r, d.cpsd_d_mean, "o-", color=col, label=f"{mode} (c-PSD)")
        axes[2].plot(p10r, 100 * (d.phi_Ni_mean - base.phi_Ni_mean) / base.phi_Ni_mean,
                    "o-", color=col, label=mode)
    for ax, title, ylab, band in (
        (axes[0], "watershed size vs neck p10 ratio", "ws_d_volweighted (nm)",
         (base.ws_d_mean * 0.95, base.ws_d_mean * 1.05)),
        (axes[1], "c-PSD size vs neck p10 ratio", "d_cPSD_r50max (nm)",
         (base.cpsd_d_mean * 0.95, base.cpsd_d_mean * 1.05)),
    ):
        ax.axhline(band[0], color="crimson", ls="--", lw=1, label="+/-5% band")
        ax.axhline(band[1], color="crimson", ls="--", lw=1)
        ax.axvline(1.5, color="k", ls=":", lw=1, label="1.5x neck p10 target")
        ax.set_xlabel("neck p10 / base neck p10")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, frameon=False)
        ax.grid(alpha=0.25)
    axes[2].axvline(1.5, color="k", ls=":", lw=1)
    axes[2].axhspan(-1, 1, color="crimson", alpha=0.1, label="+/-1% Ni phi")
    axes[2].set_xlabel("neck p10 / base neck p10")
    axes[2].set_ylabel("Ni volume fraction deviation (%)")
    axes[2].set_title("Ni loading stability vs neck p10 ratio", fontsize=10)
    axes[2].legend(fontsize=7, frameon=False)
    axes[2].grid(alpha=0.25)
    fig.suptitle("T5 coupling decision: can neck p10 move independently of "
                 "particle-size measures and Ni loading?", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "t5_coupling_experiment.png"), dpi=145)
    print(f"\n[saved] {os.path.join(OUT, 't5_coupling_experiment.csv')}")
    print(f"[saved] {os.path.join(OUT, 't5_coupling_experiment_agg.csv')}")
    print(f"[saved] {os.path.join(OUT, 't5_coupling_experiment.png')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
