"""
PROJECT 2 — STEP 1 RE-QUALIFICATION (particulate YSZ generator).

Runs under PREREGISTRATION_V2_1.md (committed e62f30b), with the kill test
(8e77035) approved and p_sinter frozen at 0.955 / 0.523 / 0.416.

5 structure seeds x 3 analog classes = 15 structures. Gates G1-a .. G1-i.

NO DAMAGE OPERATOR is implemented or applied. O1/O2/O3 do not exist.
cmlib/damage.py and cmlib/synth.py are unmodified.

What is new versus the original (random-field) Step 1:
  * YSZ comes from the particulate generator with explicit contacts and a
    frozen sintering yield, not from a thresholded smoothed random field.
  * G1-a solves `neck_scale` PER SEED (amendment A3) instead of freezing it at
    the seed-0 value, which caused the single G1-a failure last time.
  * G1-h is scored on MEAN YSZ EDT, not p50 (amendment A2) -- p50 failed on an
    exact voxel-quantization tie, not on direction.
  * G1-i is scored on Q_YSZ = 1 - P_span with a ratio requirement and an
    anti-outlier clause (amendment A1); raw and filtered cluster counts are
    recorded but NEVER gate.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.percolation import percolation_summary  # noqa: E402
from cmlib.pnm import extract_ni_network  # noqa: E402
from cmlib.project2 import (  # noqa: E402
    LATTICE_NN, draw_sintered, solve_r_ysz_for_phi_v2, ysz_lattice_geometry,
)
from cmlib.synth import (  # noqa: E402
    draw_valid_base_widths, platform_v2_lattice_geometry, rasterize,
)

OUT = os.path.join(ROOT, "out", "project2")
VOXEL_NM, AXIS, CONN = 20.0, 0, 6
LATTICE, YSZ_JITTER, YSZ_SEED_BASE = "FCC", 0.02, 91000
STRUCT_SEEDS = [0, 1, 2, 3, 4]
FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO = 0.20, (4, 6), (12, 20), 2.5

ANALOGS = {
    "fine":   dict(nlat_z=8, nlat_xy=6, pitch=24, R=10.5, margin=8, jitter=0.15,
                   phi_ni_real=0.322, phi_ysz_real=0.421, a_ysz=40.0,
                   w_ysz=8.0, ns0=0.747),
    "medium": dict(nlat_z=7, nlat_xy=5, pitch=30, R=12.1, margin=8, jitter=0.15,
                   phi_ni_real=0.250, phi_ysz_real=0.388, a_ysz=46.0,
                   w_ysz=10.0, ns0=0.807),
    "coarse": dict(nlat_z=6, nlat_xy=4, pitch=36, R=14.0, margin=8, jitter=0.15,
                   phi_ni_real=0.229, phi_ysz_real=0.384, a_ysz=53.0,
                   w_ysz=11.0, ns0=0.988),
}
ORDER = ["fine", "medium", "coarse"]
P_SINTER = json.load(open(os.path.join(OUT, "p_sinter_calibrated.json")))


def solve_neck_scale(centres, pairs, shape, R, w_raw, phi_target, ns0,
                     max_iter=8):
    """A3: bisect `neck_scale` per seed to hit this analog's Phi_Ni target."""
    dom = float(np.prod(shape))
    lo, hi = max(0.2, ns0 - 0.35), ns0 + 0.35
    best, iters = None, 0
    for it in range(max_iter):
        mid = 0.5 * (lo + hi)
        ni = rasterize(centres, pairs, R,
                       np.maximum(np.round(w_raw * mid), 2.0), shape)
        phi = ni.sum() / dom
        iters = it + 1
        if best is None or abs(phi - phi_target) < abs(best[2] - phi_target):
            best = (mid, ni, phi)
        if phi > phi_target:
            hi = mid
        else:
            lo = mid
    return best[0], best[1], best[2], iters


def one(name, cfg, seed):
    t0 = time.time()
    centres, pairs, shape = platform_v2_lattice_geometry(
        cfg["nlat_z"], cfg["nlat_xy"], cfg["pitch"], cfg["R"], cfg["margin"],
        cfg["jitter"], np.random.default_rng(999 + seed))
    w_raw, sub_seed, n_att, _log = draw_valid_base_widths(
        len(pairs), seed, FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO)
    ns, ni, phi_ni, ni_iters = solve_neck_scale(
        centres, pairs, shape, cfg["R"], w_raw, cfg["phi_ni_real"], cfg["ns0"])
    w = np.maximum(np.round(w_raw * ns), 2.0)

    # --- YSZ: particulate generator at the frozen sintering yield ---
    p_s = P_SINTER[name]
    yc, yp = ysz_lattice_geometry(
        shape, cfg["a_ysz"], YSZ_JITTER,
        np.random.default_rng(YSZ_SEED_BASE + 137 * ORDER.index(name) + seed),
        LATTICE)
    sint = draw_sintered(len(yp), p_s,
                         np.random.default_rng(YSZ_SEED_BASE + 7717 + seed))
    nn = LATTICE_NN[LATTICE] * cfg["a_ysz"]
    r_y, ysz, phi_ysz, y_iters = solve_r_ysz_for_phi_v2(
        yc, yp, sint, cfg["w_ysz"], shape, ni, cfg["phi_ysz_real"],
        r_lo=1.0, r_hi=0.5 * nn * 0.999, max_iter=12)

    nip = percolation_summary(ni, axis=AXIS, connectivity=CONN,
                              check_other_axes=False)
    yp_rep = percolation_summary(ysz, axis=AXIS, connectivity=CONN,
                                 check_other_axes=False)

    edt = ndi.distance_transform_edt(ysz) * VOXEL_NM
    ev = edt[ysz]
    edt_mean, edt_p50 = float(ev.mean()), float(np.percentile(ev, 50))
    del edt, ev

    G, _d, _e = extract_ni_network(ni, spacing_nm=(VOXEL_NM,) * 3, axis=AXIS,
                                   connectivity=CONN)
    if G is None or G.number_of_nodes() == 0:
        node_d, n_nodes = np.nan, 0
    else:
        vols = np.array([d["volume_nm3"] for _, d in G.nodes(data=True)])
        eqd = np.array([d["equiv_diam_nm"] for _, d in G.nodes(data=True)])
        node_d, n_nodes = float((eqd * vols).sum() / vols.sum()), G.number_of_nodes()

    lab, ncl = ndi.label(ysz, structure=ndi.generate_binary_structure(3, 1))
    cnt = np.bincount(lab.ravel())
    cnt[0] = 0
    grain_vox = (4.0 / 3.0) * np.pi * r_y ** 3
    rest = np.sort(cnt[1:])[::-1][1:] if ncl > 1 else np.array([], dtype=int)
    filt = rest[rest > 0.1 * grain_vox]
    del lab

    dom = float(np.prod(shape))
    neck_p10 = float(np.percentile(w, 10))
    r = dict(
        analog=name, struct_seed=seed,
        # ---- A3 mandatory per-structure log ----
        neck_scale_solved=ns, ni_solve_iters=ni_iters, ni_solve_failure="",
        phi_ni=phi_ni, phi_ni_target=cfg["phi_ni_real"],
        phi_ni_dev_pct=100 * (phi_ni - cfg["phi_ni_real"]) / cfg["phi_ni_real"],
        phi_ysz=phi_ysz, phi_ysz_target=cfg["phi_ysz_real"],
        phi_ysz_dev_pct=100 * (phi_ysz - cfg["phi_ysz_real"]) / cfg["phi_ysz_real"],
        ysz_solve_iters=y_iters,
        n_ni_particles=len(centres), n_ni_necks=len(pairs),
        n_ysz_grains=len(yc), n_candidate_contacts=len(yp),
        n_sintered_contacts=int(sint.sum()),
        mean_coord_sintered=2.0 * int(sint.sum()) / max(len(yc), 1),
        snow_node_d_volwt_nm=node_d, snow_n_nodes=n_nodes,
        ysz_grain_d_nm=2.0 * r_y * VOXEL_NM, r_ysz_vox=r_y,
        gap_vox=nn - 2.0 * r_y,
        # ---- generation parameters / seeds ----
        p_sinter=p_s, lattice=LATTICE, a_ysz_vox=cfg["a_ysz"],
        w_ysz_vox=cfg["w_ysz"], ysz_jitter=YSZ_JITTER,
        ni_jitter=cfg["jitter"], pitch_vox=cfg["pitch"], R_ni_vox=cfg["R"],
        geom_rng_seed=999 + seed, width_seed_accepted=sub_seed,
        width_draw_attempts=n_att,
        ysz_geom_seed=YSZ_SEED_BASE + 137 * ORDER.index(name) + seed,
        ysz_sinter_seed=YSZ_SEED_BASE + 7717 + seed,
        shape_z=shape[0], shape_y=shape[1], shape_x=shape[2],
        domain_voxels=int(dom), voxel_nm=VOXEL_NM, axis=AXIS, connectivity=CONN,
        # ---- measured ----
        neck_p10_vox=neck_p10, neck_p50_vox=float(np.percentile(w, 50)),
        ni_P_span=nip["P_span"], ni_n_clusters=nip["n_clusters"],
        ni_percolates=nip["percolates"],
        ysz_P_span=yp_rep["P_span"], ysz_P_reach=yp_rep["P_reach"],
        ysz_P_largest=yp_rep["P_largest"], ysz_n_clusters=yp_rep["n_clusters"],
        ysz_percolates=yp_rep["percolates"],
        Q_ysz=1.0 - yp_rep["P_span"],
        ysz_edt_mean_nm=edt_mean, ysz_edt_p50_nm=edt_p50,
        filt_n=int(filt.size), filt_threshold_vox=float(0.1 * grain_vox),
        seconds=round(time.time() - t0, 1))
    # per-structure gates
    r["G1a_phi_ni"] = bool(abs(r["phi_ni_dev_pct"]) <= 2.0)
    r["G1b_phi_ysz"] = bool(abs(r["phi_ysz_dev_pct"]) <= 2.0)
    r["G1c_ni_perc"] = bool(nip["percolates"] and nip["n_clusters"] == 1
                            and abs(nip["P_span"] - 1.0) < 1e-9)
    r["G1d_ysz_perc"] = bool(yp_rep["P_span"] > 0.0)
    r["G1e_particles"] = bool(len(centres) >= 30)
    r["G1f_neck_p10"] = bool(neck_p10 >= 4.0)
    return r


def main():
    print("=" * 78)
    print("PROJECT 2 — STEP 1 RE-QUALIFICATION (particulate YSZ)")
    print("=" * 78)
    print(f"  p_sinter (frozen): {P_SINTER}")
    rows = []
    for name in ORDER:
        for seed in STRUCT_SEEDS:
            r = one(name, ANALOGS[name], seed)
            rows.append(r)
            print(f"  {name:7s} s{seed} ns={r['neck_scale_solved']:.4f} "
                  f"phiNi={r['phi_ni']:.4f}({r['phi_ni_dev_pct']:+.2f}%) "
                  f"phiY={r['phi_ysz']:.4f}({r['phi_ysz_dev_pct']:+.2f}%) "
                  f"NiP={r['ni_P_span']:.4f} Q={r['Q_ysz']:.5f} "
                  f"EDT={r['ysz_edt_mean_nm']:.1f} node_d={r['snow_node_d_volwt_nm']:.0f} "
                  f"d_ysz={r['ysz_grain_d_nm']:.0f} [{r['seconds']}s]")
            pd.DataFrame(rows).to_csv(
                os.path.join(OUT, "step1_requal.csv"), index=False)
    df = pd.DataFrame(rows)

    med = df.groupby("analog").median(numeric_only=True).loc[ORDER]
    inc = lambda v: bool(v[0] < v[1] < v[2])  # noqa: E731

    g1g = inc(med["snow_node_d_volwt_nm"].values)
    g1h = inc(med["ysz_edt_mean_nm"].values)
    q = med["Q_ysz"].values
    g1i_order = inc(q)
    g1i_ratio = float(q[2] / q[0]) if q[0] > 0 else np.inf
    # anti-outlier: dropping any single seed from any class must preserve both
    robust = True
    for a in ORDER:
        for s in STRUCT_SEEDS:
            sub = df[~((df.analog == a) & (df.struct_seed == s))]
            m = sub.groupby("analog").Q_ysz.median().loc[ORDER].values
            if not (inc(m) and (m[2] / m[0]) >= 10):
                robust = False
    g1i = bool(g1i_order and g1i_ratio >= 10 and robust)

    print("\n" + "=" * 78)
    print("CLASS MEDIANS")
    print("=" * 78)
    print(med[["phi_ni", "phi_ysz", "Q_ysz", "ysz_edt_mean_nm",
               "snow_node_d_volwt_nm", "ysz_grain_d_nm", "neck_scale_solved",
               "ysz_P_span", "ysz_n_clusters", "filt_n"]].to_string())

    print("\n" + "=" * 78)
    print("GATES")
    print("=" * 78)
    per = ["G1a_phi_ni", "G1b_phi_ysz", "G1c_ni_perc", "G1d_ysz_perc",
           "G1e_particles", "G1f_neck_p10"]
    for g in per:
        n_ok = int(df[g].sum())
        print(f"  {g:16s} {n_ok}/15 " + ("PASS" if n_ok == 15 else "** FAIL **"))
    print(f"  G1g_size_order   {'PASS' if g1g else '** FAIL **'}  "
          f"({', '.join(f'{v:.0f}' for v in med['snow_node_d_volwt_nm'])} nm)")
    print(f"  G1h_ysz_lenscale {'PASS' if g1h else '** FAIL **'}  mean EDT "
          f"({', '.join(f'{v:.1f}' for v in med['ysz_edt_mean_nm'])} nm)  "
          f"[p50 was {', '.join(f'{v:.1f}' for v in med['ysz_edt_p50_nm'])}]")
    print(f"  G1i_ysz_fragility {'PASS' if g1i else '** FAIL **'}")
    print(f"       Q medians {', '.join(f'{v:.5f}' for v in q)}  "
          f"(real 0.00109, 0.01197, 0.07542)")
    print(f"       ordering {'PASS' if g1i_order else 'FAIL'};  "
          f"ratio {g1i_ratio:.1f}x {'PASS' if g1i_ratio >= 10 else 'FAIL'} "
          f"(need >=10);  anti-outlier {'PASS' if robust else 'FAIL'}")

    gates = pd.DataFrame(
        [dict(gate=g, scope="per-structure", n_pass=int(df[g].sum()), n=15,
              passed=bool(df[g].all())) for g in per] +
        [dict(gate="G1g_size_order", scope="class median", n_pass=int(g1g),
              n=1, passed=g1g),
         dict(gate="G1h_ysz_lenscale_meanEDT", scope="class median",
              n_pass=int(g1h), n=1, passed=g1h),
         dict(gate="G1i_order", scope="class median", n_pass=int(g1i_order),
              n=1, passed=g1i_order),
         dict(gate="G1i_ratio_ge10", scope="class median",
              n_pass=int(g1i_ratio >= 10), n=1, passed=bool(g1i_ratio >= 10)),
         dict(gate="G1i_anti_outlier", scope="leave-one-seed-out",
              n_pass=int(robust), n=1, passed=robust),
         dict(gate="G1i_combined", scope="class median", n_pass=int(g1i), n=1,
              passed=g1i)])
    gates.to_csv(os.path.join(OUT, "step1_requal_gates.csv"), index=False)
    med.to_csv(os.path.join(OUT, "step1_requal_class_medians.csv"))

    allp = bool(df[per].all().all() and g1g and g1h and g1i)
    print("\n" + "=" * 78)
    print("VERDICT: " + ("ALL GATES PASS" if allp else "** GATE FAILURE **"))
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
