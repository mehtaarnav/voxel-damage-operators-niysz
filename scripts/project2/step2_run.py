"""
Step 2 runner: R1 scoping, main arm (O1/O2/O3), and the R3 matched-fragility
control arm.

Protocol per the advisor's Step 2 instruction and PREREGISTRATION_V2_1:
  * integer bisection, bracket [1,20], expand-only, final width <= 1
  * 5 structure seeds per analog class, 3 damage seeds per structure
  * damage-seed averaging mandatory; group differences < 1.0 round are
    UNRESOLVED, never passes
  * boundary rules: transition at <=1 round for all seeds, or still percolating
    at 20 rounds for all seeds -> stop and report, do not retune

Phase scored per operator: O1 and O2 damage Ni -> Ni percolation loss.
O3 damages YSZ -> YSZ percolation loss. NOTE (structural, reported not worked
around): no single operator can satisfy C3, because C3 requires C1 (Ni) and C2
(YSZ) simultaneously and no single operator touches both phases. The
pre-registered pair {O1 or O2} x O3 is therefore required for C3 by
construction, not as a rescue.

    python step2_run.py --mode r1scope
    python step2_run.py --mode main --op O1|O2|O3
    python step2_run.py --mode control
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cmlib.damage2 import apply_o1, apply_o2, apply_o3, tpb_density_um2  # noqa: E402
from cmlib.percolation import percolation_summary  # noqa: E402
from step2_build_cache import CFG, ORDER, SEEDS, build, load, path  # noqa: E402

OUT = os.path.join(ROOT, "out", "project2")
VOXEL_NM, AXIS, CONN = 20.0, 0, 6
DAMAGE_SEEDS = [300, 301, 302]
BRACKET_LO, BRACKET_HI = 1, 20
R1_GRID = [1, 2, 3, 5, 8, 13, 20]
OP_PHASE = {"O1": "Ni", "O2": "Ni", "O3": "YSZ"}


def damage(op, st, n, dseed, region_slices=None):
    """Apply one operator at intensity n. Returns (ni, ysz, info)."""
    if op == "O1":
        ni, info = apply_o1(st["ni"], st["ysz"], n, dseed)
        return ni, st["ysz"], info
    if op == "O2":
        ni, info = apply_o2(st["ni"], st["regions"], st["throat_conns"],
                            st["throat_diam"], n, dseed,
                            region_slices=region_slices)
        return ni, st["ysz"], info
    ysz, intact, info = apply_o3(st["ysz_centres"], st["ysz_pairs"],
                                 st["sintered"], st["r_ysz"], st["w_ysz"],
                                 st["shape"], st["ni"], n, dseed)
    return st["ni"], ysz, info


def spans(op, st, n, dseed, rs=None):
    ni, ysz, _ = damage(op, st, n, dseed, rs)
    m = ni if OP_PHASE[op] == "Ni" else ysz
    return percolation_summary(m, axis=AXIS, connectivity=CONN,
                               check_other_axes=False)["P_span"] > 0.0


def bisect(op, st, dseed, rs=None):
    """Integer bisection for the loss-of-percolation intensity."""
    lo, hi = BRACKET_LO, BRACKET_HI
    if not spans(op, st, lo, dseed, rs):
        return float(lo) - 0.5, "floor"
    if spans(op, st, hi, dseed, rs):
        return float(hi) + 0.5, "ceiling"
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if spans(op, st, mid, dseed, rs):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), "ok"


def metrics(op, st, n, dseed, rs=None):
    ni, ysz, info = damage(op, st, n, dseed, rs)
    dom = float(np.prod(st["shape"]))
    nip = percolation_summary(ni, axis=AXIS, connectivity=CONN,
                              check_other_axes=False)
    yp = percolation_summary(ysz, axis=AXIS, connectivity=CONN,
                             check_other_axes=False)
    return dict(
        phi_ni_post=ni.sum() / dom, phi_ysz_post=ysz.sum() / dom,
        ni_P_span=nip["P_span"], ni_P_reach=nip["P_reach"],
        ni_P_largest=nip["P_largest"], ni_n_clusters=nip["n_clusters"],
        ysz_P_span=yp["P_span"], ysz_P_largest=yp["P_largest"],
        ysz_n_clusters=yp["n_clusters"],
        tpb_um2=tpb_density_um2(ni, ysz, VOXEL_NM),
        ni_vol_loss_frac=1.0 - ni.sum() / max(st["ni"].sum(), 1),
        ysz_vol_loss_frac=1.0 - ysz.sum() / max(st["ysz"].sum(), 1))


def r1scope():
    rows = []
    for op in ("O1", "O2", "O3"):
        for name in ORDER:
            st = load(name, 0)
            rs = ndi.find_objects(st["regions"].astype(np.int32)) if op == "O2" else None
            for n in R1_GRID:
                t0 = time.time()
                ni, ysz, _ = damage(op, st, n, DAMAGE_SEEDS[0], rs)
                m = ni if OP_PHASE[op] == "Ni" else ysz
                r = percolation_summary(m, axis=AXIS, connectivity=CONN,
                                        check_other_axes=False)
                rows.append(dict(op=op, analog=name, n_rounds=n,
                                 phase=OP_PHASE[op], P_span=r["P_span"],
                                 seconds=round(time.time() - t0, 1)))
                print(f"  {op} {name:7s} n={n:2d} P_span={r['P_span']:.4f} "
                      f"[{rows[-1]['seconds']}s]", flush=True)
                pd.DataFrame(rows).to_csv(
                    os.path.join(OUT, "step2_r1_scope.csv"), index=False)
    df = pd.DataFrame(rows)
    nsec = {}
    for op in ("O1", "O2", "O3"):
        brackets = []
        for name in ORDER:
            s = df[(df.op == op) & (df.analog == name)].sort_values("n_rounds")
            lost = s[s.P_span <= 0.0]
            if not len(lost):
                continue
            first = int(lost.n_rounds.iloc[0])
            prev = s[s.n_rounds < first].n_rounds
            brackets.append((int(prev.max()) if len(prev) else 0, first))
        if not brackets:
            nsec[op] = 5
        else:
            a, b = min(brackets, key=lambda t: t[1])   # shallowest transition
            nsec[op] = int(np.clip((a + b) // 2, 2, 15))
        print(f"  n_secondary[{op}] = {nsec[op]}  (brackets {brackets})")
    json.dump(nsec, open(os.path.join(OUT, "step2_n_secondary.json"), "w"),
              indent=1)
    return 0


def run_arm(op, p_tag="main", outfile=None, nsec=None):
    outfile = outfile or f"step2_{op}_{p_tag}.csv"
    rows = []
    for name in ORDER:
        for seed in SEEDS:
            st = load(name, seed, p_tag)
            rs = ndi.find_objects(st["regions"].astype(np.int32)) if op == "O2" else None
            pris_q = 1.0 - percolation_summary(
                st["ysz"], axis=AXIS, connectivity=CONN,
                check_other_axes=False)["P_span"]
            for ds in DAMAGE_SEEDS:
                t0 = time.time()
                mid, flag = bisect(op, st, ds, rs)
                n_tr = int(np.ceil(mid))
                m = metrics(op, st, n_tr, ds, rs)
                sec = spans(op, st, nsec[op], ds, rs) if nsec else None
                ni2, ysz2, _ = damage(op, st, nsec[op], ds, rs) if nsec else (None, None, None)
                if nsec:
                    mm = ni2 if OP_PHASE[op] == "Ni" else ysz2
                    sec_pspan = percolation_summary(
                        mm, axis=AXIS, connectivity=CONN,
                        check_other_axes=False)["P_span"]
                else:
                    sec_pspan = np.nan
                rows.append(dict(op=op, arm=p_tag, analog=name,
                                 struct_seed=seed, damage_seed=ds,
                                 transition_midpoint=mid, bracket_flag=flag,
                                 n_transition=n_tr, phase=OP_PHASE[op],
                                 pristine_Q_ysz=pris_q,
                                 n_secondary=nsec[op] if nsec else np.nan,
                                 retained_P_span_at_nsec=sec_pspan,
                                 seconds=round(time.time() - t0, 1), **m))
                print(f"  {op} {p_tag} {name:7s} s{seed} d{ds} "
                      f"mid={mid:.1f}({flag}) TPB={m['tpb_um2']:.2f} "
                      f"NiP={m['ni_P_span']:.3f} YP={m['ysz_P_span']:.3f} "
                      f"sec={sec_pspan:.3f} [{rows[-1]['seconds']}s]",
                      flush=True)
                pd.DataFrame(rows).to_csv(os.path.join(OUT, outfile),
                                          index=False)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True,
                    choices=["r1scope", "main", "control", "buildctrl"])
    ap.add_argument("--op", default="O1")
    a = ap.parse_args()
    if a.mode == "r1scope":
        return r1scope()
    nsec = json.load(open(os.path.join(OUT, "step2_n_secondary.json")))
    if a.mode == "buildctrl":
        pc = json.load(open(os.path.join(OUT, "step2_p_control.json")))
        for name in ORDER:
            for seed in SEEDS:
                if not os.path.exists(path(name, seed, "ctrl")):
                    build(name, seed, pc[name], "ctrl")
                    print(f"  built control {name} s{seed} p={pc[name]}",
                          flush=True)
        return 0
    if a.mode == "control":
        return run_arm("O3", "ctrl", "step2_O3_control.csv", nsec)
    return run_arm(a.op, "main", f"step2_{a.op}_main.csv", nsec)


if __name__ == "__main__":
    sys.exit(main())
