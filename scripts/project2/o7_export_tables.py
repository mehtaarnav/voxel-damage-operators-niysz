"""Export the sequential-swap results as committed CSVs.

The manuscript's convention is that every plotted value reads from a committed
table rather than from a live run. This script produces the three tables behind
the new figures:

  o7_trajectory.csv      specific surface area, TPB and R_Ni against accepted
                         moves, each indexed to its pristine value
  o7_da_histogram.csv    the distribution of dA over accepted moves
  o7_counterfactual.csv  TPB ratio when only the contact-adjacent moves, or
                         only the remaining moves, are re-applied to pristine

All three are measured on the fine anode ROI at seed 300.
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cmlib.damage2 import tpb_density_um2                      # noqa: E402
from cmlib.percolation import percolating_mask                 # noqa: E402
from cmlib.seqgreedy import SeqGreedy, ST7                     # noqa: E402
import primer_figures as pf                                    # noqa: E402

OUT = os.path.join(ROOT, "out", "project2")
AXIS, CONN = 2, 6


def main():
    ni0, ysz, vox = pf.load_roi("fine")
    s0 = pf.spec_surf(ni0)
    t0 = tpb_density_um2(ni0, ysz, vox)
    n0 = int(ni0.sum())
    r0 = int(percolating_mask(ni0, axis=AXIS, connectivity=CONN).sum()) / n0
    print(f"pristine  S_spec {s0:.5f}  TPB {t0:.4f}  R_Ni {r0:.4f}")

    op = SeqGreedy(ni0, ysz, seed=300)
    op.dA_log = []
    total = int(round(0.03 * op.n_surf0)) * 5
    rows = [dict(accepted=0, s_spec=s0, tpb=t0, r_ni=r0,
                 s_ratio=1.0, tpb_ratio=1.0, r_ratio=1.0)]
    done = 0
    tic = time.time()
    for cp in np.unique(np.linspace(0, total, 13).astype(int))[1:]:
        op.run(int(cp) - done)
        done = int(cp)
        s = pf.spec_surf(op.ni)
        t = tpb_density_um2(op.ni, ysz, vox)
        r = int(percolating_mask(op.ni, axis=AXIS,
                                 connectivity=CONN).sum()) / n0
        rows.append(dict(accepted=op.accepted, s_spec=s, tpb=t, r_ni=r,
                         s_ratio=s / s0, tpb_ratio=t / t0, r_ratio=r / r0))
        print(f"  {op.accepted:>7}  S {s/s0:.4f}  TPB {t/t0:.3f}")
    pd.DataFrame(rows).to_csv(f"{OUT}/o7_trajectory.csv", index=False)
    print(f"trajectory written ({time.time()-tic:.0f}s)")

    dA = np.array(op.dA_log)
    vals, counts = np.unique(dA, return_counts=True)
    pd.DataFrame({"dA": vals, "count": counts,
                  "share": counts / counts.sum()}).to_csv(
        f"{OUT}/o7_da_histogram.csv", index=False)
    print(f"histogram written; neutral share {(dA == 0).mean():.4f}")

    ni1 = op.ni
    added, removed = ni1 & ~ni0, ni0 & ~ni1
    nearY = ndi.binary_dilation(ysz, ST7)
    nc = ni0.copy(); nc[added & nearY] = True; nc[removed & nearY] = False
    nf = ni0.copy(); nf[added & ~nearY] = True; nf[removed & ~nearY] = False
    cf = [
        dict(case="all_moves", tpb=tpb_density_um2(ni1, ysz, vox)),
        dict(case="contact_adjacent_only", tpb=tpb_density_um2(nc, ysz, vox)),
        dict(case="away_from_contact_only", tpb=tpb_density_um2(nf, ysz, vox)),
    ]
    for r in cf:
        r["tpb_ratio"] = r["tpb"] / t0
    cf.append(dict(case="pristine", tpb=t0, tpb_ratio=1.0))
    pd.DataFrame(cf).to_csv(f"{OUT}/o7_counterfactual.csv", index=False)
    print("counterfactual written:",
          {r["case"]: round(r["tpb_ratio"], 3) for r in cf})

    stats = dict(
        accepted=int(op.accepted), neutral=int(op.neutral),
        neutral_share=float((dA == 0).mean()),
        added=int(added.sum()), removed=int(removed.sum()),
        removed_near_ysz=int((removed & nearY).sum()),
        added_near_ysz=int((added & nearY).sum()),
    )
    pd.DataFrame([stats]).to_csv(f"{OUT}/o7_move_stats.csv", index=False)
    print("stats:", stats)


if __name__ == "__main__":
    main()
