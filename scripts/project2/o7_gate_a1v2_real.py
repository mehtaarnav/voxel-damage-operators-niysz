"""Gate A1v2 on REAL ROIs with the sequential greedy operator.

PREREG_O5V2_OPTIONB.md specifies Gate A1v2 as "1 ROI/anode, n = 1,3,5, seeds
300/301/302" with five conditions:
  (i) |dPhi_Ni| <= 0.005; (ii) S_spec(1) < S_spec(0) strict+monotonic;
  (iii) TPB(n) <= TPB(0); (iv) R_Ni non-increasing; (v) YSZ untouched.

The run that closed the agglomeration route used a synthetic structure on which
(i), (iii), (iv) and (v) are all vacuous (see O7_O5V2B_RERUN_REPORT.md). Real
ROIs make TPB(0) > 0 and R_Ni(0) > 0, so every condition carries information.

Conventions follow scripts/project2/c1real_run.py: AXIS=2, CONN=6,
SIDE = {fine 8, medium 8, coarse 12} um, k = 0.03 * (surface voxels) per round.

Usage:  python o7_gate_a1v2_real.py [fine|medium|coarse ...]
"""
import os
import sys
import time

import numpy as np
import tifffile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.damage2 import tpb_density_um2                     # noqa: E402
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE     # noqa: E402
from cmlib.io import label_histogram, slice_paths             # noqa: E402
from cmlib.percolation import percolating_mask, percolation_summary  # noqa: E402
from cmlib.phases import assign_labels                        # noqa: E402
from cmlib.roi import tile_rois                               # noqa: E402
from cmlib.seqgreedy import SeqGreedy                         # noqa: E402

AXIS, CONN = 2, 6
SIDE = {"fine": 8.0, "medium": 8.0, "coarse": 12.0}
SEEDS = (300, 301, 302)
NS = (1, 3, 5)


def load(folder, mapping, r):
    ps = slice_paths(folder)
    sh = (r["z1"] - r["z0"], r["y1"] - r["y0"], r["x1"] - r["x0"])
    ni = np.empty(sh, bool)
    ysz = np.empty(sh, bool)
    for i, z in enumerate(range(r["z0"], r["z1"])):
        a = tifffile.imread(ps[z])[r["y0"]:r["y1"], r["x0"]:r["x1"]]
        ni[i] = a == mapping["Ni"]
        ysz[i] = a == mapping["YSZ"]
    return ni, ysz


def spec_surf(m):
    s = 0
    for ax in range(3):
        for sh in (1, -1):
            s += int((m & ~np.roll(m, sh, axis=ax)).sum())
    return s / max(int(m.sum()), 1)


def main(which):
    sample = {s[2]: s for s in SAMPLES if s[3] == "pristine"}
    for grain in which:
        k = sample[grain]
        nz, ny, nx_ = k[6], k[5], k[4]
        vz, vy, vx = k[9], k[8], k[7]
        mapping = assign_labels(label_histogram(k[1])["counts"],
                                ZENODO_LABEL_NOTE[k[0]])
        r = tile_rois(nz, ny, nx_, vz, vy, vx, SIDE[grain], max_rois=1)[0]

        t0 = time.time()
        ni0, ysz = load(k[1], mapping, r)
        vox = float((vz * vy * vx) ** (1 / 3))
        n0 = int(ni0.sum())
        s0 = spec_surf(ni0)
        phi0 = n0 / ni0.size
        tpb0 = tpb_density_um2(ni0, ysz, vox)
        p0 = percolation_summary(ni0, axis=AXIS, connectivity=CONN,
                                 check_other_axes=False)["P_span"]
        rni0 = int(percolating_mask(ni0, axis=AXIS,
                                    connectivity=CONN).sum()) / n0

        print(f"\n=== {grain.upper()}  ROI {ni0.shape}  "
              f"{ni0.size/1e6:.1f} Mvox  loaded in {time.time()-t0:.0f}s ===")
        print(f"  Phi_Ni={phi0:.6f}  S_spec={s0:.5f}  TPB(0)={tpb0:.4f} um^-2")
        print(f"  P_span={p0:.4f}  R_Ni(0)={rni0:.4f}")
        # PREREG_RNI_METRIC.md mandatory sanity check
        print(f"  SANITY R_Ni(0) == P_span: {abs(rni0 - p0) < 1e-9}")
        print(f"  informative? TPB(0)>0: {tpb0 > 0}   R_Ni(0)>0: {rni0 > 0}")

        probe = SeqGreedy(ni0, ysz, seed=0)
        k_round = int(round(0.03 * probe.n_surf0))
        print(f"  surface voxels={probe.n_surf0}  k/round={k_round}")
        del probe

        hdr = (f"  {'seed':>5} {'n':>2} {'moves':>9} {'neutral':>9} "
               f"{'S_spec':>9} {'dS':>10} {'dPhi':>9} {'TPB':>8} "
               f"{'R_Ni':>7} {'YSZ':>5} {'sec':>6}")
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        res = {}
        for seed in SEEDS:
            op = SeqGreedy(ni0, ysz, seed=seed)
            done = 0
            for n in NS:
                target = k_round * n
                t1 = time.time()
                op.run(target - done)
                done = target
                ni = op.ni
                s1 = spec_surf(ni)
                tpb = tpb_density_um2(ni, ysz, vox)
                rni = int(percolating_mask(ni, axis=AXIS,
                                           connectivity=CONN).sum()) / n0
                dphi = abs(int(ni.sum()) / ni.size - phi0)
                ysz_ok = not bool((ni & ysz).any())
                res[(seed, n)] = (s1, dphi, tpb, rni, ysz_ok)
                print(f"  {seed:>5} {n:>2} {op.accepted:>9} {op.neutral:>9} "
                      f"{s1:>9.5f} {s1-s0:>+10.5f} {dphi:>9.2e} "
                      f"{tpb:>8.4f} {rni:>7.4f} {str(ysz_ok):>5} "
                      f"{time.time()-t1:>6.0f}")
            del op

        ok_i = all(v[1] <= 0.005 for v in res.values())
        strict = all(res[(s, 1)][0] < s0 for s in SEEDS)
        mono = all(res[(s, 1)][0] >= res[(s, 3)][0] >= res[(s, 5)][0]
                   for s in SEEDS)
        ok_iii = all(v[2] <= tpb0 + 1e-12 for v in res.values())
        ok_iv = all(res[(s, 1)][3] <= rni0 + 1e-12
                    and res[(s, 3)][3] <= res[(s, 1)][3] + 1e-12
                    and res[(s, 5)][3] <= res[(s, 3)][3] + 1e-12
                    for s in SEEDS)
        ok_v = all(v[4] for v in res.values())
        print(f"\n  GATE A1v2 [{grain}]")
        print(f"    (i)   |dPhi_Ni| <= 0.005            : "
              f"{'PASS' if ok_i else 'FAIL'}")
        print(f"    (ii)  S_spec strict + monotonic     : "
              f"{'PASS' if strict and mono else 'FAIL'}"
              f"  (strict={strict}, monotonic={mono})")
        print(f"    (iii) TPB(n) <= TPB(0)              : "
              f"{'PASS' if ok_iii else 'FAIL'}")
        print(f"    (iv)  R_Ni non-increasing           : "
              f"{'PASS' if ok_iv else 'FAIL'}")
        print(f"    (v)   YSZ untouched                 : "
              f"{'PASS' if ok_v else 'FAIL'}")


if __name__ == "__main__":
    args = sys.argv[1:] or ["fine", "medium", "coarse"]
    main(args)
