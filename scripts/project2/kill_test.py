"""
PROJECT 2 — KILL TEST for the particulate YSZ generator (Option A).

Runs under PREREGISTRATION_V2_1.md (committed e62f30b) §5. Gates K0–K5.
No damage operator is implemented or applied anywhere in this file.

ARCHITECTURE NOTE, recorded before any result. The first architecture tried was
an SC (simple-cubic, body-centred-on-Ni) YSZ lattice. It is INFEASIBLE for the
fine analog by a closed-form packing bound, not by experiment:

    Phi_YSZ_max ~= f_pack * (1 - Phi_Ni);  SC gives 0.5236 * 0.678 = 0.355,
    against a required 0.421 -- i.e. 119% of its own ceiling.

An empirical probe agreed (radius bisection saturated at the touching limit,
gap 0.01 voxels, Phi stuck at 0.331 vs 0.388 target for medium). SC was
therefore REJECTED BEFORE K0 on the derivation, rather than run and failed.
BCC (cap 0.461) and FCC (cap 0.502) both clear all three analogs; FCC is used
because at equal inter-grain gap it yields the most grains per domain.

Usage:
    python kill_test.py --mode scope     # 1 seed, wide p grid, fast
    python kill_test.py --mode full      # 3 seeds, frozen grid, gates K0-K5
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.percolation import percolation_summary  # noqa: E402
from cmlib.project2 import (  # noqa: E402
    LATTICE_NN, draw_sintered, max_phi_ysz, solve_r_ysz_for_phi_v2,
    ysz_lattice_geometry,
)
from cmlib.synth import (  # noqa: E402
    draw_valid_base_widths, platform_v2_lattice_geometry, rasterize,
)

OUT = os.path.join(ROOT, "out", "project2")
VOXEL_NM, AXIS, CONN = 20.0, 0, 6
LATTICE = "FCC"
YSZ_JITTER = 0.02          # fraction of nearest-neighbour distance
YSZ_SEED_BASE = 91000

# Ni analog parameters: DESIGN_MEMO 1.2 (unchanged). YSZ cube edge a_Y scales
# with the analog's Ni particle diameter, anchored at fine = 40 vox, which is
# the smallest a_Y giving a >= 1.5 voxel inter-grain gap for the fine analog.
ANALOGS = {
    "fine":   dict(nlat_z=8, nlat_xy=6, pitch=24, R=10.5, margin=8, jitter=0.15,
                   neck_scale=0.747, phi_ni_real=0.322, phi_ysz_real=0.421,
                   a_ysz=40.0, w_ysz=8.0),
    "medium": dict(nlat_z=7, nlat_xy=5, pitch=30, R=12.1, margin=8, jitter=0.15,
                   neck_scale=0.807, phi_ni_real=0.250, phi_ysz_real=0.388,
                   a_ysz=46.0, w_ysz=10.0),
    "coarse": dict(nlat_z=6, nlat_xy=4, pitch=36, R=14.0, margin=8, jitter=0.15,
                   neck_scale=0.988, phi_ni_real=0.229, phi_ysz_real=0.384,
                   a_ysz=53.0, w_ysz=11.0),
}
ORDER = ["fine", "medium", "coarse"]
FROZEN_P = {"fine": 0.70, "medium": 0.56, "coarse": 0.42}   # SC-derived priors


def build_ni(cfg, seed):
    centres, pairs, shape = platform_v2_lattice_geometry(
        cfg["nlat_z"], cfg["nlat_xy"], cfg["pitch"], cfg["R"], cfg["margin"],
        cfg["jitter"], np.random.default_rng(999 + seed))
    w, _s, _n, _l = draw_valid_base_widths(len(pairs), seed, 0.20, (4, 6),
                                           (12, 20), 2.5)
    ni = rasterize(centres, pairs, cfg["R"],
                   np.maximum(np.round(w * cfg["neck_scale"]), 2.0), shape)
    return ni, shape


def one(name, cfg, seed, p_sinter, ni, shape):
    t0 = time.time()
    rng_g = np.random.default_rng(YSZ_SEED_BASE + 137 * ORDER.index(name) + seed)
    yc, yp = ysz_lattice_geometry(shape, cfg["a_ysz"], YSZ_JITTER, rng_g,
                                  LATTICE)
    rng_s = np.random.default_rng(YSZ_SEED_BASE + 7717 + seed)
    sint = draw_sintered(len(yp), p_sinter, rng_s)
    nn = LATTICE_NN[LATTICE] * cfg["a_ysz"]
    r_hi = 0.5 * nn * 0.999
    r, mask, phi, nit = solve_r_ysz_for_phi_v2(
        yc, yp, sint, cfg["w_ysz"], shape, ni, cfg["phi_ysz_real"],
        r_lo=1.0, r_hi=r_hi, max_iter=12)
    rep = percolation_summary(mask, axis=AXIS, connectivity=CONN,
                              check_other_axes=False)
    edt = ndi.distance_transform_edt(mask) * VOXEL_NM
    ev = edt[mask]
    edt_mean = float(ev.mean()) if ev.size else 0.0
    edt_p50 = float(np.percentile(ev, 50)) if ev.size else 0.0
    del edt, ev

    # K5 filtered diagnostics
    lab, ncl = ndi.label(mask, structure=ndi.generate_binary_structure(3, 1))
    cnt = np.bincount(lab.ravel())
    cnt[0] = 0
    grain_vox = (4.0 / 3.0) * np.pi * r ** 3
    rest = np.sort(cnt[1:])[::-1][1:] if ncl > 1 else np.array([], dtype=int)
    filt = rest[rest > 0.1 * grain_vox]
    del lab

    dom = float(np.prod(shape))
    return dict(
        analog=name, seed=seed, p_sinter=p_sinter, lattice=LATTICE,
        a_ysz_vox=cfg["a_ysz"], w_ysz_vox=cfg["w_ysz"], ysz_jitter=YSZ_JITTER,
        n_grains=len(yc), n_candidate_contacts=len(yp),
        n_sintered=int(sint.sum()),
        mean_coord_sintered=2.0 * int(sint.sum()) / max(len(yc), 1),
        r_ysz_vox=r, grain_d_nm=2.0 * r * VOXEL_NM,
        gap_vox=nn - 2.0 * r, solve_iters=nit,
        phi_ysz=phi, phi_ysz_target=cfg["phi_ysz_real"],
        phi_ysz_rel_dev=abs(phi - cfg["phi_ysz_real"]) / cfg["phi_ysz_real"],
        phi_ni=float(ni.sum() / dom),
        P_span=rep["P_span"], P_reach=rep["P_reach"],
        P_largest=rep["P_largest"], Q_ysz=1.0 - rep["P_span"],
        n_clusters=rep["n_clusters"], percolates=rep["percolates"],
        edt_mean_nm=edt_mean, edt_p50_nm=edt_p50,
        filt_n=int(filt.size), filt_volfrac=float(filt.sum() / max(cnt[1:].sum(), 1)),
        filt_threshold_vox=float(0.1 * grain_vox),
        domain_vox=int(dom), seconds=round(time.time() - t0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["scope", "full"], default="scope")
    args = ap.parse_args()

    if args.mode == "scope":
        seeds, grid, analogs = [0], [0.0, 0.15, 0.25, 0.35, 0.50, 1.0], ORDER
        outfile = "killtest_scope.csv"
    else:
        seeds, grid, analogs = [0, 1, 2], None, ORDER
        outfile = "killtest_full.csv"

    print("=" * 78)
    print(f"KILL TEST ({args.mode}) — lattice={LATTICE}, jitter={YSZ_JITTER}")
    print("=" * 78)
    for n in ORDER:
        c = ANALOGS[n]
        print(f"  {n:7s} cap={max_phi_ysz(LATTICE, c['phi_ni_real']):.3f} "
              f"need={c['phi_ysz_real']:.3f} "
              f"({100*c['phi_ysz_real']/max_phi_ysz(LATTICE, c['phi_ni_real']):.0f}% of cap)"
              f"  a_ysz={c['a_ysz']}")

    if grid is None:
        gp = pd.read_csv(os.path.join(OUT, "killtest_scope.csv"))
        grid = sorted(gp.p_sinter.unique())

    rows = []
    for name in analogs:
        cfg = ANALOGS[name]
        for seed in seeds:
            ni, shape = build_ni(cfg, seed)
            for p in grid:
                r = one(name, cfg, seed, p, ni, shape)
                rows.append(r)
                print(f"  {name:7s} s{seed} p={p:.2f}  r={r['r_ysz_vox']:.2f} "
                      f"gap={r['gap_vox']:.2f}  phi={r['phi_ysz']:.4f} "
                      f"(dev {100*r['phi_ysz_rel_dev']:.2f}%)  "
                      f"P_span={r['P_span']:.4f} Q={r['Q_ysz']:.5f} "
                      f"ncl={r['n_clusters']} (grains {r['n_grains']})  "
                      f"[{r['seconds']}s]")
                pd.DataFrame(rows).to_csv(os.path.join(OUT, outfile),
                                          index=False)
    print(f"\n[saved] {os.path.join(OUT, outfile)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
