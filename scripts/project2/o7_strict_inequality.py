"""Is the area-neutral plateau the sole driver of TPB manufacture?

Under dA <= 0 the sequential swap operator inflates TPB density 3.7-5.0x on
real ROIs while reducing specific surface area monotonically. About 80% of its
accepted moves carry dA = 0 exactly and are unpriced by the area criterion.

Test: forbid them. Accept only dA < 0 (strict), so every move must strictly
lower area, and measure TPB.

CONFOUND THIS CONTROLS FOR. The strict operator exhausts its admissible moves
far sooner, so a low TPB inflation could simply mean "it barely moved". Two
comparisons are therefore reported:

  A. strict (dA < 0) run to exhaustion, M accepted moves
  B. standard (dA <= 0) TRUNCATED to the same M accepted moves

If the plateau is the driver, A shows little inflation while B, at the identical
move count, shows substantially more.
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
from cmlib.percolation import percolating_mask                # noqa: E402
from cmlib.phases import assign_labels                        # noqa: E402
from cmlib.roi import tile_rois                               # noqa: E402
from cmlib.seqgreedy import SeqGreedy                         # noqa: E402

AXIS, CONN = 2, 6
SIDE = {"fine": 8.0, "medium": 8.0, "coarse": 12.0}
SEED = 300


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
    print(f"{'anode':<8} {'rule':<22} {'accepted':>9} {'neutral':>9} "
          f"{'S_spec':>9} {'dS':>10} {'TPB':>8} {'TPB/TPB0':>9} {'R_Ni':>7}")
    print("-" * 100)
    for grain in which:
        k = sample[grain]
        mapping = assign_labels(label_histogram(k[1])["counts"],
                                ZENODO_LABEL_NOTE[k[0]])
        r = tile_rois(k[6], k[5], k[4], k[9], k[8], k[7],
                      SIDE[grain], max_rois=1)[0]
        ni0, ysz = load(k[1], mapping, r)
        vox = float((k[9] * k[8] * k[7]) ** (1 / 3))
        n0 = int(ni0.sum())
        s0 = spec_surf(ni0)
        tpb0 = tpb_density_um2(ni0, ysz, vox)
        rni0 = int(percolating_mask(ni0, axis=AXIS, connectivity=CONN).sum()) / n0
        print(f"{grain:<8} {'pristine':<22} {'':>9} {'':>9} "
              f"{s0:>9.5f} {'':>10} {tpb0:>8.4f} {1.0:>9.3f} {rni0:>7.4f}")

        # A. strict, run to exhaustion
        t0 = time.time()
        op = SeqGreedy(ni0, ysz, seed=SEED, strict=True)
        budget = 40 * int(round(0.03 * op.n_surf0))     # generous upper bound
        op.run(budget)
        M = op.accepted
        niA = op.ni
        sA, tA = spec_surf(niA), tpb_density_um2(niA, ysz, vox)
        rA = int(percolating_mask(niA, axis=AXIS, connectivity=CONN).sum()) / n0
        del op
        print(f"{grain:<8} {'A strict dA<0':<22} {M:>9} {0:>9} "
              f"{sA:>9.5f} {sA-s0:>+10.5f} {tA:>8.4f} {tA/tpb0:>9.3f} "
              f"{rA:>7.4f}")

        # B. standard, truncated to the same accepted-move count
        opB = SeqGreedy(ni0, ysz, seed=SEED, strict=False)
        while opB.accepted < M:
            if not opB.step():
                break
        niB = opB.ni
        sB, tB = spec_surf(niB), tpb_density_um2(niB, ysz, vox)
        rB = int(percolating_mask(niB, axis=AXIS, connectivity=CONN).sum()) / n0
        print(f"{grain:<8} {'B dA<=0, matched M':<22} {opB.accepted:>9} "
              f"{opB.neutral:>9} {sB:>9.5f} {sB-s0:>+10.5f} {tB:>8.4f} "
              f"{tB/tpb0:>9.3f} {rB:>7.4f}   [{time.time()-t0:.0f}s]")
        del opB
        print()


if __name__ == "__main__":
    main(sys.argv[1:] or ["fine", "medium", "coarse"])
