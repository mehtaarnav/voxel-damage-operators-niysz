"""
Step 2 preparation: rebuild the 15 qualified Step-1 structures ONCE, cache
everything the operators need, and compute the R2 pristine TPB baseline.

Caching matters: the main arm is ~180 bisections and a structure rebuild costs
~90 s, so rebuilding per evaluation would dominate everything. Cached per
structure: Ni mask, YSZ mask, the YSZ contact graph (centres/pairs/sintered/
radius), and the Ni SNOW regions + throat table that O2 needs.

R2 also lands here: pristine TPB on all 15 structures, corrected units (um^-2),
using the unit-tested `cmlib.damage2.tpb_density_um2`.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.damage2 import tpb_density_um2  # noqa: E402
from cmlib.pnm import extract_ni_network  # noqa: E402
from cmlib.project2 import (  # noqa: E402
    LATTICE_NN, draw_sintered, solve_r_ysz_for_phi_v2, ysz_lattice_geometry,
)
from cmlib.synth import (  # noqa: E402
    draw_valid_base_widths, platform_v2_lattice_geometry, rasterize,
)

OUT = os.path.join(ROOT, "out", "project2")
CACHE = os.path.join(OUT, "_cache")
os.makedirs(CACHE, exist_ok=True)
VOXEL_NM, AXIS, CONN = 20.0, 0, 6
LATTICE, YSZ_JITTER, YSZ_SEED_BASE = "FCC", 0.02, 91000
ORDER = ["fine", "medium", "coarse"]
SEEDS = [0, 1, 2, 3, 4]

CFG = {
    "fine":   dict(nlat_z=8, nlat_xy=6, pitch=24, R=10.5, margin=8, jitter=0.15,
                   phi_ysz=0.421, a_ysz=40.0, w_ysz=8.0),
    "medium": dict(nlat_z=7, nlat_xy=5, pitch=30, R=12.1, margin=8, jitter=0.15,
                   phi_ysz=0.388, a_ysz=46.0, w_ysz=10.0),
    "coarse": dict(nlat_z=6, nlat_xy=4, pitch=36, R=14.0, margin=8, jitter=0.15,
                   phi_ysz=0.384, a_ysz=53.0, w_ysz=11.0),
}
P_SINTER = json.load(open(os.path.join(OUT, "p_sinter_calibrated.json")))


def path(name, seed, p_tag="main"):
    return os.path.join(CACHE, f"{name}_s{seed}_{p_tag}.npz")


def build(name, seed, p_sinter, p_tag="main"):
    cfg = CFG[name]
    req = pd.read_csv(os.path.join(OUT, "step1_requal.csv"))
    ns = float(req[(req.analog == name) &
                   (req.struct_seed == seed)].neck_scale_solved.iloc[0])
    centres, pairs, shape = platform_v2_lattice_geometry(
        cfg["nlat_z"], cfg["nlat_xy"], cfg["pitch"], cfg["R"], cfg["margin"],
        cfg["jitter"], np.random.default_rng(999 + seed))
    w_raw, _s, _n, _l = draw_valid_base_widths(len(pairs), seed, 0.20, (4, 6),
                                               (12, 20), 2.5)
    ni = rasterize(centres, pairs, cfg["R"],
                   np.maximum(np.round(w_raw * ns), 2.0), shape)

    yc, yp = ysz_lattice_geometry(
        shape, cfg["a_ysz"], YSZ_JITTER,
        np.random.default_rng(YSZ_SEED_BASE + 137 * ORDER.index(name) + seed),
        LATTICE)
    sint = draw_sintered(len(yp), p_sinter,
                         np.random.default_rng(YSZ_SEED_BASE + 7717 + seed))
    nn = LATTICE_NN[LATTICE] * cfg["a_ysz"]
    r_y, ysz, phi_y, _ = solve_r_ysz_for_phi_v2(
        yc, yp, sint, cfg["w_ysz"], shape, ni, cfg["phi_ysz"],
        r_lo=1.0, r_hi=0.5 * nn * 0.999, max_iter=12)

    G, diag, extras = extract_ni_network(ni, spacing_nm=(VOXEL_NM,) * 3,
                                         axis=AXIS, connectivity=CONN)
    regions = extras["regions"].astype(np.int16)
    if G is not None and G.number_of_edges():
        conns, diam = [], []
        for u, v, d in G.edges(data=True):
            conns.append((u, v))
            diam.append(d["neck_nm"])
        conns = np.asarray(conns, dtype=np.int32)
        diam = np.asarray(diam, dtype=np.float64)
    else:
        conns = np.zeros((0, 2), np.int32)
        diam = np.zeros(0)

    keys = np.asarray(sorted(yc.keys()), dtype=np.int64)
    ycoords = np.asarray([yc[k] for k in keys], dtype=np.float64)
    ypairs = np.asarray([(a, b) for a, b in yp], dtype=np.int64)

    np.savez_compressed(
        path(name, seed, p_tag), ni=ni, ysz=ysz, regions=regions,
        throat_conns=conns, throat_diam=diam, ysz_keys=keys,
        ysz_coords=ycoords, ysz_pairs=ypairs, sintered=sint,
        r_ysz=np.array([r_y]), w_ysz=np.array([cfg["w_ysz"]]),
        shape=np.asarray(shape), phi_ysz=np.array([phi_y]),
        neck_scale=np.array([ns]), p_sinter=np.array([p_sinter]))
    return ni, ysz, shape, phi_y, conns.shape[0]


def load(name, seed, p_tag="main"):
    z = np.load(path(name, seed, p_tag))
    centres = {int(k): tuple(c) for k, c in zip(z["ysz_keys"], z["ysz_coords"])}
    pairs = [(int(a), int(b)) for a, b in z["ysz_pairs"]]
    return dict(ni=z["ni"], ysz=z["ysz"], regions=z["regions"],
                throat_conns=z["throat_conns"], throat_diam=z["throat_diam"],
                ysz_centres=centres, ysz_pairs=pairs,
                sintered=z["sintered"], r_ysz=float(z["r_ysz"][0]),
                w_ysz=float(z["w_ysz"][0]), shape=tuple(z["shape"]),
                phi_ysz=float(z["phi_ysz"][0]))


def main():
    rows = []
    for name in ORDER:
        for seed in SEEDS:
            t0 = time.time()
            ni, ysz, shape, phi_y, nthr = build(name, seed, P_SINTER[name])
            tpb = tpb_density_um2(ni, ysz, VOXEL_NM)
            dom = float(np.prod(shape))
            rows.append(dict(analog=name, struct_seed=seed,
                             phi_ni=ni.sum() / dom, phi_ysz=phi_y,
                             tpb_um2=tpb, n_ni_throats=nthr,
                             seconds=round(time.time() - t0, 1)))
            print(f"  {name:7s} s{seed}  TPB={tpb:7.3f} um^-2  "
                  f"phiNi={ni.sum()/dom:.4f} phiYSZ={phi_y:.4f} "
                  f"throats={nthr} [{rows[-1]['seconds']}s]", flush=True)
            pd.DataFrame(rows).to_csv(
                os.path.join(OUT, "step2_r2_tpb_baseline.csv"), index=False)
    df = pd.DataFrame(rows)
    print("\nR2 PRISTINE TPB BASELINE (frozen, um^-2)")
    med = df.groupby("analog").tpb_um2.median().loc[ORDER]
    real = {"fine": 3.624, "medium": 2.109, "coarse": 1.473}
    for a in ORDER:
        v = df[df.analog == a].tpb_um2
        print(f"  {a:7s} median {med[a]:7.3f}  (seeds {v.min():.3f}-{v.max():.3f})"
              f"  real {real[a]:.3f}  ratio {med[a]/real[a]:.1f}x")
    print(f"  ordering fine>medium>coarse: "
          f"{bool(med['fine'] > med['medium'] > med['coarse'])}   "
          f"fine/coarse {med['fine']/med['coarse']:.2f} (real "
          f"{real['fine']/real['coarse']:.2f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
