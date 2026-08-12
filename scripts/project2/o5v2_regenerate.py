"""Regenerate out/project2/o5v2_area_barrier.csv.

SUPERSEDES scripts/project2/o5v2_transcribe.py, which copied numbers out of the
run reports. The `greedy_area` rows in that transcription record
`acceptance_rate = 0.0`, produced by a defect: cmlib/damage2.py convolved with
STRUCT6 (centre INCLUDED, sum 7) to count Ni 6-neighbours, biasing Ni sites by
+1 and turning the frozen predicate nN(a) <= nN(b) into nN(a) < nN(b). Every
area-neutral move was rejected. See out/project2/O7_O5V2B_RERUN_REPORT.md.

The `greedy_area` rows are therefore RECOMPUTED here with the corrected kernel
(cmlib.damage2.NB6), on the recovered original structure.

The `curvature_ranked` rows are RETAINED AS TRANSCRIBED. Option A is provably
unaffected by the defect: it ranks Ni sites by (nmax - 2*nb) and pore sites by
nb in two SEPARATE orderings, and a constant +1 on every Ni site cannot reorder
either. Its results also depend on a jitter RNG whose original draw is not
recoverable, so recomputing them would change numbers for no correctness gain.

NECK COLUMN. The original reports state neck = 63 voxels pristine. That value is
not reproducible from the recovered structure under any definition tried: the
3x3 bar contributes 99 voxels outside the spheres, over 11 slices of 9. The
original definition is not recorded anywhere in the repo. Rather than invent a
match, regenerated rows carry a stated definition (`free_span_3x3_column`) and
the transcribed rows are marked `original_unrecoverable`.
"""
import csv
import os
import sys

import numpy as np
import scipy.ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.damage2 import apply_o5v2b, NB6            # noqa: E402

OUT = os.path.join(ROOT, "out", "project2", "o5v2_area_barrier.csv")
A = "o5v2_report.md"
B = "o5v2_optionB_report.md"
THIS = "o5v2_regenerate.py"


def original_structure():
    z, y, x = np.ogrid[:60, :60, :60]
    big = ((((z - 15) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100) |
           (((z - 45) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100))
    big = np.array(big)
    spheres = big.copy()
    big[15:46, 29:32, 29:32] = True
    ysz = np.zeros_like(big)
    ysz[:, :, 55:] = True
    return big, ysz, spheres


def s_spec(ni):
    nb = ndi.convolve(ni.astype(np.int16), NB6.astype(np.int16),
                      mode="constant", cval=0)
    return float((6 - nb)[ni].sum()) / max(int(ni.sum()), 1)


def neck_free_span(ni, spheres):
    """Ni voxels in the 3x3 connecting column that lie outside both spheres."""
    return int((ni[:, 29:32, 29:32] & ~spheres[:, 29:32, 29:32]).sum())


# curvature_ranked: unaffected by the defect, retained verbatim from the reports
TRANSCRIBED = [
    ("curvature_ranked", "6",  0, 0.0, 0.45052, 63, "", "", "", A,
     "original_unrecoverable"),
    ("curvature_ranked", "6",  1, 0.0, 0.45195, 63, "", "", "", A,
     "original_unrecoverable"),
    ("curvature_ranked", "6",  3, 0.0, 0.45219, 63, "", "", "", A,
     "original_unrecoverable"),
    ("curvature_ranked", "6",  5, 0.0, 0.45362, 63, "", "", "", A,
     "original_unrecoverable"),
    ("curvature_ranked", "26", 0, 0.0, 0.45052, 63, "", "", "", A,
     "original_unrecoverable"),
    ("curvature_ranked", "26", 1, 0.0, 0.45195, 57, "", "", "", A,
     "original_unrecoverable"),
    ("curvature_ranked", "26", 3, 0.0, 0.44431, 15, "", "", "", A,
     "original_unrecoverable"),
    ("curvature_ranked", "26", 5, 0.0, 0.44216,  0, "", "", "", A,
     "original_unrecoverable"),
]


def main():
    ni0, ysz, spheres = original_structure()
    rows = list(TRANSCRIBED)

    print("Regenerating greedy_area rows with the corrected NB6 kernel.")
    print(f"pristine S_spec = {s_spec(ni0):.5f} "
          f"(original reports 0.45052)")
    print(f"{'n':>3} {'proposed':>9} {'accepted':>9} {'rate':>7} "
          f"{'S_spec':>9} {'neck':>6} {'volerr':>8}")
    for n in (0, 1, 3, 5, 8):
        if n == 0:
            ni, info = ni0, dict(proposed=0, accepted=0, acceptance_rate="",
                                 volume_error=0.0)
        else:
            ni, info = apply_o5v2b(ni0, ysz, n, seed=300)
        ss = s_spec(ni)
        nk = neck_free_span(ni, spheres)
        rate = info["acceptance_rate"]
        rows.append(("greedy_area", "6", n, info["volume_error"],
                     round(ss, 5), nk, rate, info["proposed"],
                     info["accepted"], THIS, "free_span_3x3_column"))
        rs = f"{rate:.4f}" if rate != "" else "n/a"
        print(f"{n:>3} {info['proposed']:>9} {info['accepted']:>9} {rs:>7} "
              f"{ss:>9.5f} {nk:>6} {info['volume_error']:>8.1e}")

    with open(OUT, "w", newline="", encoding="utf8") as fh:
        w = csv.writer(fh)
        w.writerow(["operator", "stencil", "n_rounds", "volume_error",
                    "S_spec", "neck_voxels", "acceptance_rate",
                    "proposed", "accepted", "source", "neck_metric"])
        w.writerows(rows)
    print("\nwrote", OUT)
    print("NOTE: the greedy stencil column is '6', not the '26' recorded in the "
          "transcription;\n      apply_o5v2b has no conn26 argument and always "
          "used the 6-neighbour kernel.")


if __name__ == "__main__":
    main()
