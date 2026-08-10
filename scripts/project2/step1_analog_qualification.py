"""
PROJECT 2 — STEP 1: analog qualification.

Builds 5 structure seeds x 3 analog classes (fine / medium / coarse) at the
EXACT parameters frozen in out/project2/DESIGN_MEMO.md sec 1.2, applies the
A2.3 YSZ sigma-scaling rule (cmlib/project2.py), and runs gates G1-a .. G1-i.

NOT DONE HERE, per instruction: no damage operator (O1/O2/O3) is implemented or
applied; `cmlib/damage.py` and `cmlib/synth.py` are unmodified (add_ysz_pore
already exposes `smooth_sigma_vox`, so A2.3 needed no edit to either).

NO TUNING RULE. neck_scale is frozen per analog at the seed-0 solved value from
the design memo; sigma is fixed by the A2.3 rule. If a gate fails, it is
reported as a failure. Nothing here re-solves, re-seeds, or re-tunes on the
basis of a gate outcome.

Gates
-----
G1-a Phi_Ni within +-2% of the real anode value
G1-b Phi_YSZ within +-2% of the real anode value
G1-c Ni percolates, single cluster, P_span == 1.000
G1-d YSZ percolates (P_span > 0); exact value recorded, no threshold  [revised]
G1-e Ni particle count >= 30
G1-f Ni neck p10 >= 4 voxels
G1-g Ni particle size orders fine < medium < coarse on SNOW node size (class)
G1-h YSZ length scale orders fine < medium < coarse on YSZ EDT p50    (class)
G1-i pristine YSZ P_span orders fine > medium > coarse AND
     pristine YSZ n_clusters orders fine < medium < coarse            (class)

G1-a..f are per-structure; G1-g..i are cross-class and evaluated on class means.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.damage import add_ysz_pore  # noqa: E402
from cmlib.percolation import percolation_summary  # noqa: E402
from cmlib.pnm import extract_ni_network  # noqa: E402
from cmlib.project2 import ysz_sigma_for_analog  # noqa: E402
from cmlib.synth import (  # noqa: E402
    draw_valid_base_widths, platform_v2_lattice_geometry, rasterize,
)

OUT = os.path.join(ROOT, "out", "project2")
os.makedirs(OUT, exist_ok=True)

VOXEL_NM = 20.0
AXIS, CONN = 0, 6
STRUCT_SEEDS = [0, 1, 2, 3, 4]
GEOM_SEED_BASE = 999
YSZ_SEED_BASE = 5000

# Frozen neck-width mixture (inherited from Platform v2, unchanged)
FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO = 0.20, (4, 6), (12, 20), 2.5

# DESIGN_MEMO sec 1.2, verbatim. neck_scale is the seed-0 solved value.
ANALOGS = {
    "fine":   dict(nlat_z=8, nlat_xy=6, pitch=24, R=10.5, margin=8,
                   jitter=0.15, neck_scale=0.747,
                   phi_ni_real=0.322, phi_ysz_real=0.421,
                   ysz_frac_of_rest=0.421 / (0.421 + 0.254)),
    "medium": dict(nlat_z=7, nlat_xy=5, pitch=30, R=12.1, margin=8,
                   jitter=0.15, neck_scale=0.807,
                   phi_ni_real=0.250, phi_ysz_real=0.388,
                   ysz_frac_of_rest=0.388 / (0.388 + 0.362)),
    "coarse": dict(nlat_z=6, nlat_xy=4, pitch=36, R=14.0, margin=8,
                   jitter=0.15, neck_scale=0.988,
                   phi_ni_real=0.229, phi_ysz_real=0.384,
                   ysz_frac_of_rest=0.384 / (0.384 + 0.387)),
}
ORDER = ["fine", "medium", "coarse"]


def build_structure(cfg, seed):
    centres, pairs, shape = platform_v2_lattice_geometry(
        cfg["nlat_z"], cfg["nlat_xy"], cfg["pitch"], cfg["R"], cfg["margin"],
        cfg["jitter"], np.random.default_rng(GEOM_SEED_BASE + seed))
    w_raw, sub_seed, n_att, _log = draw_valid_base_widths(
        len(pairs), seed, FRAC_WEAK, WEAK_RANGE, NORMAL_RANGE, MIN_RATIO)
    w = np.maximum(np.round(w_raw * cfg["neck_scale"]), 2.0)
    ni = rasterize(centres, pairs, cfg["R"], w, shape)
    return centres, pairs, shape, w, ni, sub_seed, n_att


def qualify_one(name, cfg, seed):
    t0 = time.time()
    centres, pairs, shape, w, ni, sub_seed, n_att = build_structure(cfg, seed)
    dom = int(np.prod(shape))
    d_particle_nm = 2.0 * cfg["R"] * VOXEL_NM
    sigma = ysz_sigma_for_analog(d_particle_nm)
    ysz_seed = YSZ_SEED_BASE + 100 * ORDER.index(name) + seed

    _vol, ysz = add_ysz_pore(ni, seed=ysz_seed,
                             ysz_frac_of_rest=cfg["ysz_frac_of_rest"],
                             smooth_sigma_vox=sigma)

    nip = percolation_summary(ni, axis=AXIS, connectivity=CONN,
                              check_other_axes=False)
    yp = percolation_summary(ysz, axis=AXIS, connectivity=CONN,
                             check_other_axes=False)

    # YSZ length scale (G1-h): EDT over YSZ voxels, in nm
    edt = ndi.distance_transform_edt(ysz) * VOXEL_NM
    ev = edt[ysz]
    ysz_edt_p50 = float(np.percentile(ev, 50))
    ysz_edt_mean = float(ev.mean())
    ysz_edt_p90 = float(np.percentile(ev, 90))
    del edt, ev

    # Ni SNOW node size (G1-g)
    G, diag, _extras = extract_ni_network(ni, spacing_nm=(VOXEL_NM,) * 3,
                                          axis=AXIS, connectivity=CONN)
    if G is None or G.number_of_nodes() == 0:
        node_d_volwt = np.nan
        node_d_mean = np.nan
        n_nodes = 0
    else:
        vols = np.array([d["volume_nm3"] for _, d in G.nodes(data=True)])
        eqd = np.array([d["equiv_diam_nm"] for _, d in G.nodes(data=True)])
        node_d_volwt = float((eqd * vols).sum() / vols.sum())
        node_d_mean = float(eqd.mean())
        n_nodes = int(G.number_of_nodes())

    phi_ni = ni.sum() / dom
    phi_ysz = ysz.sum() / dom
    neck_p10_vox = float(np.percentile(w, 10))

    row = dict(
        analog=name, struct_seed=seed,
        # --- generation parameters (fully recorded) ---
        nlat_z=cfg["nlat_z"], nlat_xy=cfg["nlat_xy"], pitch_vox=cfg["pitch"],
        R_vox=cfg["R"], margin_vox=cfg["margin"], jitter_frac=cfg["jitter"],
        neck_scale=cfg["neck_scale"], frac_weak=FRAC_WEAK,
        weak_lo=WEAK_RANGE[0], weak_hi=WEAK_RANGE[1],
        normal_lo=NORMAL_RANGE[0], normal_hi=NORMAL_RANGE[1],
        min_p50_p10_ratio=MIN_RATIO,
        geom_rng_seed=GEOM_SEED_BASE + seed,
        width_seed_requested=seed, width_seed_accepted=sub_seed,
        width_draw_attempts=n_att,
        ysz_seed=ysz_seed, ysz_sigma_vox=sigma,
        ysz_frac_of_rest=cfg["ysz_frac_of_rest"],
        voxel_nm=VOXEL_NM, axis=AXIS, connectivity=CONN,
        # --- geometry ---
        shape_z=shape[0], shape_y=shape[1], shape_x=shape[2],
        domain_voxels=dom, n_particles=len(centres), n_pairs=len(pairs),
        d_particle_nm=d_particle_nm,
        neck_p10_vox=neck_p10_vox, neck_p50_vox=float(np.percentile(w, 50)),
        neck_p10_nm=neck_p10_vox * VOXEL_NM,
        neck_p50_nm=float(np.percentile(w, 50)) * VOXEL_NM,
        # --- measured ---
        phi_ni=phi_ni, phi_ni_real=cfg["phi_ni_real"],
        phi_ni_dev_pct=100.0 * (phi_ni - cfg["phi_ni_real"]) / cfg["phi_ni_real"],
        phi_ysz=phi_ysz, phi_ysz_real=cfg["phi_ysz_real"],
        phi_ysz_dev_pct=100.0 * (phi_ysz - cfg["phi_ysz_real"]) / cfg["phi_ysz_real"],
        ni_P_span=nip["P_span"], ni_n_clusters=nip["n_clusters"],
        ni_percolates=nip["percolates"], ni_P_largest=nip["P_largest"],
        ysz_P_span=yp["P_span"], ysz_P_reach=yp["P_reach"],
        ysz_P_largest=yp["P_largest"], ysz_n_clusters=yp["n_clusters"],
        ysz_percolates=yp["percolates"],
        ysz_n_clusters_per_Mvox=yp["n_clusters"] / (dom / 1e6),
        ysz_edt_p50_nm=ysz_edt_p50, ysz_edt_mean_nm=ysz_edt_mean,
        ysz_edt_p90_nm=ysz_edt_p90,
        snow_n_nodes=n_nodes, snow_node_d_volwt_nm=node_d_volwt,
        snow_node_d_mean_nm=node_d_mean,
        seconds=round(time.time() - t0, 1),
    )
    # --- per-structure gates ---
    row["G1a_phi_ni"] = bool(abs(row["phi_ni_dev_pct"]) <= 2.0)
    row["G1b_phi_ysz"] = bool(abs(row["phi_ysz_dev_pct"]) <= 2.0)
    row["G1c_ni_perc"] = bool(nip["percolates"] and nip["n_clusters"] == 1
                              and abs(nip["P_span"] - 1.0) < 1e-9)
    row["G1d_ysz_perc"] = bool(yp["P_span"] > 0.0)
    row["G1e_particles"] = bool(len(centres) >= 30)
    row["G1f_neck_p10"] = bool(neck_p10_vox >= 4.0)
    return row


def main():
    print("=" * 78)
    print("PROJECT 2 STEP 1 — analog qualification (15 structures)")
    print("=" * 78)
    print("A2.3 sigma:", {n: round(ysz_sigma_for_analog(2 * c["R"] * VOXEL_NM), 4)
                          for n, c in ANALOGS.items()})

    rows = []
    for name in ORDER:
        cfg = ANALOGS[name]
        for seed in STRUCT_SEEDS:
            r = qualify_one(name, cfg, seed)
            rows.append(r)
            print(f"\n[{name} seed {seed}] {r['seconds']}s  "
                  f"shape=({r['shape_z']},{r['shape_y']},{r['shape_x']}) "
                  f"particles={r['n_particles']} sigma={r['ysz_sigma_vox']:.3f}")
            print(f"   Phi_Ni={r['phi_ni']:.4f} ({r['phi_ni_dev_pct']:+.2f}%)  "
                  f"Phi_YSZ={r['phi_ysz']:.4f} ({r['phi_ysz_dev_pct']:+.2f}%)")
            print(f"   Ni  P_span={r['ni_P_span']:.4f} clusters={r['ni_n_clusters']}"
                  f"   neck_p10={r['neck_p10_vox']:.1f} vox"
                  f"   SNOW nodes={r['snow_n_nodes']} d_volwt="
                  f"{r['snow_node_d_volwt_nm']:.0f} nm")
            print(f"   YSZ P_span={r['ysz_P_span']:.4f} clusters={r['ysz_n_clusters']}"
                  f" ({r['ysz_n_clusters_per_Mvox']:.2f}/Mvox)"
                  f"  EDT p50={r['ysz_edt_p50_nm']:.1f} mean="
                  f"{r['ysz_edt_mean_nm']:.1f} nm")
            print("   gates: " + " ".join(
                f"{k.split('_')[0]}={'P' if r[k] else 'FAIL'}"
                for k in ("G1a_phi_ni", "G1b_phi_ysz", "G1c_ni_perc",
                          "G1d_ysz_perc", "G1e_particles", "G1f_neck_p10")))
            pd.DataFrame(rows).to_csv(
                os.path.join(OUT, "step1_analog_qualification.csv"), index=False)

    df = pd.DataFrame(rows)

    # ---------------- cross-class gates ----------------
    cls = df.groupby("analog").agg(
        snow_node_d_volwt_nm=("snow_node_d_volwt_nm", "mean"),
        ysz_edt_p50_nm=("ysz_edt_p50_nm", "mean"),
        ysz_edt_mean_nm=("ysz_edt_mean_nm", "mean"),
        ysz_P_span=("ysz_P_span", "mean"),
        ysz_n_clusters=("ysz_n_clusters", "mean"),
        ysz_n_clusters_per_Mvox=("ysz_n_clusters_per_Mvox", "mean"),
    ).loc[ORDER]

    def inc(v):
        return bool(v[0] < v[1] < v[2])

    def dec(v):
        return bool(v[0] > v[1] > v[2])

    g1g = inc(cls["snow_node_d_volwt_nm"].values)
    g1h = inc(cls["ysz_edt_p50_nm"].values)
    g1h_mean = inc(cls["ysz_edt_mean_nm"].values)
    g1i_span = dec(cls["ysz_P_span"].values)
    g1i_clust_raw = inc(cls["ysz_n_clusters"].values)
    g1i_clust_den = inc(cls["ysz_n_clusters_per_Mvox"].values)
    g1i = bool(g1i_span and g1i_clust_raw)

    print("\n" + "=" * 78)
    print("CLASS MEANS")
    print("=" * 78)
    print(cls.to_string())

    print("\n" + "=" * 78)
    print("GATES")
    print("=" * 78)
    per = ["G1a_phi_ni", "G1b_phi_ysz", "G1c_ni_perc", "G1d_ysz_perc",
           "G1e_particles", "G1f_neck_p10"]
    for g in per:
        n_ok = int(df[g].sum())
        print(f"  {g:16s} {n_ok}/15 " + ("PASS" if n_ok == 15 else "** FAIL **"))
    print(f"  G1g_size_order   {'PASS' if g1g else '** FAIL **'}"
          f"  (SNOW node d_volwt: "
          f"{', '.join(f'{v:.0f}' for v in cls['snow_node_d_volwt_nm'])} nm)")
    print(f"  G1h_ysz_lenscale {'PASS' if g1h else '** FAIL **'}"
          f"  (YSZ EDT p50: "
          f"{', '.join(f'{v:.1f}' for v in cls['ysz_edt_p50_nm'])} nm; "
          f"mean-based would be {'PASS' if g1h_mean else 'FAIL'}: "
          f"{', '.join(f'{v:.1f}' for v in cls['ysz_edt_mean_nm'])})")
    print(f"  G1i_ysz_pristine {'PASS' if g1i else '** FAIL **'}")
    print(f"       P_span decreasing  {'PASS' if g1i_span else '** FAIL **'}"
          f"  ({', '.join(f'{v:.4f}' for v in cls['ysz_P_span'])})")
    print(f"       n_clusters raw inc {'PASS' if g1i_clust_raw else '** FAIL **'}"
          f"  ({', '.join(f'{v:.0f}' for v in cls['ysz_n_clusters'])})")
    print(f"       n_clusters/Mvox    {'PASS' if g1i_clust_den else '** FAIL **'}"
          f"  ({', '.join(f'{v:.2f}' for v in cls['ysz_n_clusters_per_Mvox'])})")

    gates = pd.DataFrame([
        dict(gate="G1a_phi_ni", scope="per-structure",
             n_pass=int(df.G1a_phi_ni.sum()), n=15,
             passed=bool(df.G1a_phi_ni.all())),
        dict(gate="G1b_phi_ysz", scope="per-structure",
             n_pass=int(df.G1b_phi_ysz.sum()), n=15,
             passed=bool(df.G1b_phi_ysz.all())),
        dict(gate="G1c_ni_perc", scope="per-structure",
             n_pass=int(df.G1c_ni_perc.sum()), n=15,
             passed=bool(df.G1c_ni_perc.all())),
        dict(gate="G1d_ysz_perc", scope="per-structure",
             n_pass=int(df.G1d_ysz_perc.sum()), n=15,
             passed=bool(df.G1d_ysz_perc.all())),
        dict(gate="G1e_particles", scope="per-structure",
             n_pass=int(df.G1e_particles.sum()), n=15,
             passed=bool(df.G1e_particles.all())),
        dict(gate="G1f_neck_p10", scope="per-structure",
             n_pass=int(df.G1f_neck_p10.sum()), n=15,
             passed=bool(df.G1f_neck_p10.all())),
        dict(gate="G1g_size_order", scope="class", n_pass=int(g1g), n=1,
             passed=g1g),
        dict(gate="G1h_ysz_lenscale", scope="class", n_pass=int(g1h), n=1,
             passed=g1h),
        dict(gate="G1h_ysz_lenscale_meanEDT", scope="class (supplementary)",
             n_pass=int(g1h_mean), n=1, passed=g1h_mean),
        dict(gate="G1i_ysz_P_span_order", scope="class",
             n_pass=int(g1i_span), n=1, passed=g1i_span),
        dict(gate="G1i_ysz_nclusters_raw", scope="class",
             n_pass=int(g1i_clust_raw), n=1, passed=g1i_clust_raw),
        dict(gate="G1i_ysz_nclusters_density", scope="class (supplementary)",
             n_pass=int(g1i_clust_den), n=1, passed=g1i_clust_den),
        dict(gate="G1i_combined", scope="class", n_pass=int(g1i), n=1,
             passed=g1i),
    ])
    gates.to_csv(os.path.join(OUT, "step1_gates.csv"), index=False)
    cls.to_csv(os.path.join(OUT, "step1_class_means.csv"))

    all_pass = bool(df[per].all().all() and g1g and g1h and g1i)
    print("\n" + "=" * 78)
    print("STEP 1 VERDICT: " + ("ALL GATES PASS" if all_pass
                                else "** ONE OR MORE GATES FAILED **"))
    print("=" * 78)
    if not (g1h and g1i):
        print("G1-h and/or G1-i failed. Per the Step 1 instruction this is a")
        print("STOP condition: the YSZ placement model is insufficient to")
        print("capture the real-data morphological trend. NOT tuning sigma,")
        print("jitter or any other placement parameter to force a pass.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
