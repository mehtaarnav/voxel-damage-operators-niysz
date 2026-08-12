"""Sensitivity of the sequential greedy operator to its tie-breaking policy.

Among moves of identical dA the rule must choose one. That choice is not
neutral: area-neutral moves dominate, so which of many equal-cost moves is
taken determines how much area is actually removed.

Both policies are run here at a MATCHED move budget so the comparison is not
confounded by how long each ran. "random" is the frozen production policy;
"lifo" exists only to keep this measurement reproducible.
"""
import os
import sys

import numpy as np
import scipy.ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.seqgreedy import SeqGreedy, ST6                    # noqa: E402


def original_structure():
    z, y, x = np.ogrid[:60, :60, :60]
    big = ((((z - 15) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100) |
           (((z - 45) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100))
    big = np.array(big)
    big[15:46, 29:32, 29:32] = True
    ysz = np.zeros_like(big)
    ysz[:, :, 55:] = True
    return big, ysz


def s_spec(ni):
    nb = ndi.convolve(ni.astype(np.int16), ST6.astype(np.int16),
                      mode="constant", cval=0)
    return float((6 - nb)[ni].sum()) / max(int(ni.sum()), 1)


def main(budget=1220):
    ni0, ysz = original_structure()
    s0 = s_spec(ni0)
    print(f"recovered O5v2 structure, matched budget {budget} moves")
    print(f"pristine S_spec = {s0:.5f}\n")
    print(f"{'tiebreak':<10} {'accepted':>9} {'neutral':>9} {'S_spec':>9} "
          f"{'dS':>10}")
    print("-" * 52)
    out = {}
    for tb in ("lifo", "random"):
        op = SeqGreedy(ni0, ysz, seed=300, tiebreak=tb)
        op.run(budget)
        s1 = s_spec(op.ni)
        out[tb] = s1 - s0
        assert int(op.ni.sum()) == int(ni0.sum()), "volume not conserved"
        print(f"{tb:<10} {op.accepted:>9} {op.neutral:>9} {s1:>9.5f} "
              f"{s1-s0:>+10.5f}")
    print(f"\nratio of area reductions (random / lifo): "
          f"{out['random']/out['lifo']:.1f}x")


if __name__ == "__main__":
    main()
