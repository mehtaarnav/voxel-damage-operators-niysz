"""
PROJECT 2 — STEP 0: real-data percolation for BOTH phases, all six stacks.

Two jobs in one pass over the data:

1. PHASE 5 DEFECT RESOLUTION. `out/phase5/phase5_percolation.csv` was committed
   in 47b08ee with 1 of 6 rows and `phase5_retention.csv` empty -- an
   interrupted write, root-caused to memory exhaustion (the int32 label array
   for the largest stack is 4.46 GB on top of a 1.115 GB mask, on a 16.8 GB
   machine). Both files are regenerated here, complete, with the EXACT original
   column schema, so nothing downstream has to change.

   Two durable fixes, not just a re-run:
     - `cmlib.percolation.percolation_summary_lowmem` puts the label array in a
       disk-backed memmap and reduces it slab-wise. Definitions unchanged;
       equivalence is ASSERTED against the frozen `percolation_summary` on
       randomised volumes below, before any real stack is touched.
     - every CSV is rewritten after EACH sample, so an interruption truncates
       to a valid short file rather than an empty or malformed one.

2. STEP 0 OF THE PROJECT 2 PILOT (out/project2/DESIGN_MEMO.md sec 0.3, 4.1).
   YSZ percolation has never been measured on these stacks -- Phase 5 was
   Ni-only -- so the target ordinal signature for the YSZ half of Project 2 is
   currently a literature paraphrase rather than something we hold. This script
   measures it: P_span, P_reach, P_largest, n_clusters, n_phase_voxels for the
   YSZ phase of all six stacks, and the retained values per anode.

CONVENTIONS (inherited, not re-chosen): transport axis = x = array axis 2;
6-connectivity; full stacks at native resolution, no sub-volume.

Each stack is read twice (once per phase) rather than held once as a uint8
volume: reading is I/O-bound but the boolean mask is then the only large
resident array, which is the binding constraint on this machine.
"""

from __future__ import annotations

import gc
import os
import sys
import time

import numpy as np
import pandas as pd
import tifffile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE, ground_truth_frame  # noqa: E402
from cmlib.io import label_histogram, slice_paths  # noqa: E402
from cmlib.percolation import (  # noqa: E402
    percolation_summary, percolation_summary_lowmem,
)
from cmlib.phases import assign_labels  # noqa: E402

OUT5 = os.path.join(ROOT, "out", "phase5")
OUT2 = os.path.join(ROOT, "out", "project2")
os.makedirs(OUT5, exist_ok=True)
os.makedirs(OUT2, exist_ok=True)

AXIS, CONN = 2, 6
SCRATCH = os.environ.get("P2_SCRATCH") or os.path.join(OUT2, "_tmp")


def test_percolation_lowmem_equivalence(n_cases: int = 12) -> None:
    """Gate: the low-memory path must agree EXACTLY with the frozen function.

    Randomised volumes spanning below / at / above the site-percolation
    threshold, plus degenerate all-empty and all-full cases. Any mismatch
    raises before a single real stack is read.
    """
    rng = np.random.default_rng(20260810)
    cases = []
    for p in (0.05, 0.25, 0.3116, 0.35, 0.5, 0.8):
        cases.append(rng.random((37, 41, 29)) < p)
        cases.append(rng.random((23, 19, 53)) < p)
    cases.append(np.zeros((11, 13, 17), dtype=bool))
    cases.append(np.ones((11, 13, 17), dtype=bool))
    n_bad = 0
    for i, m in enumerate(cases[:max(n_cases, len(cases))]):
        for ax in (0, 1, 2):
            a = percolation_summary(m, axis=ax, connectivity=CONN,
                                    check_other_axes=True)
            b = percolation_summary_lowmem(m, axis=ax, connectivity=CONN,
                                           check_other_axes=True,
                                           workdir=SCRATCH, slab=7)
            if set(a) != set(b):
                raise AssertionError(f"key mismatch case {i} axis {ax}")
            for k in a:
                va, vb = a[k], b[k]
                same = (abs(va - vb) < 1e-12 if isinstance(va, float)
                        else va == vb)
                if not same:
                    n_bad += 1
                    print(f"  MISMATCH case {i} axis {ax} {k}: {va} != {vb}")
    if n_bad:
        raise AssertionError(f"{n_bad} mismatches: low-memory path is NOT "
                             f"equivalent; refusing to proceed")
    print(f"[gate] low-memory percolation == frozen percolation_summary on "
          f"{len(cases)} volumes x 3 axes: PASS")


def load_phase_mask(folder: str, label_value: int) -> np.ndarray:
    """Full-stack boolean mask for one phase, built slice by slice."""
    ps = slice_paths(folder)
    a0 = tifffile.imread(ps[0])
    out = np.empty((len(ps), a0.shape[0], a0.shape[1]), dtype=bool)
    out[0] = a0 == label_value
    for i, p in enumerate(ps[1:], start=1):
        out[i] = tifffile.imread(p) == label_value
    return out


def analyse(mask: np.ndarray) -> dict:
    return percolation_summary_lowmem(mask, axis=AXIS, connectivity=CONN,
                                      check_other_axes=True,
                                      workdir=SCRATCH)


def write_all(ni_rows, ysz_rows):
    """Rewrite every CSV from scratch after each sample (crash-safe)."""
    if ni_rows:
        df = pd.DataFrame(ni_rows)
        df.to_csv(os.path.join(OUT5, "phase5_percolation.csv"), index=False)
        _write_retention(df, os.path.join(OUT5, "phase5_retention.csv"),
                         published=True)
    if ysz_rows:
        dy = pd.DataFrame(ysz_rows)
        dy.to_csv(os.path.join(OUT2, "step0_ysz_percolation.csv"), index=False)
        _write_retention(dy, os.path.join(OUT2, "step0_ysz_retention.csv"),
                         published=False)


def _write_retention(df, path, published: bool):
    ret = []
    for g in ("fine", "medium", "coarse"):
        a = df[(df.grain == g) & (df.state == "pristine")]
        b = df[(df.grain == g) & (df.state == "degraded")]
        if not (len(a) and len(b)):
            continue
        a, b = a.iloc[0], b.iloc[0]
        row = dict(
            grain=g,
            P_span_pre=a.P_span, P_span_post=b.P_span,
            P_span_retained=b.P_span / a.P_span if a.P_span else np.nan,
            P_reach_pre=a.P_reach, P_reach_post=b.P_reach,
            P_reach_retained=b.P_reach / a.P_reach if a.P_reach else np.nan,
            P_largest_pre=a.P_largest, P_largest_post=b.P_largest,
            P_largest_retained=(b.P_largest / a.P_largest
                                if a.P_largest else np.nan),
            n_clusters_pre=a.n_clusters, n_clusters_post=b.n_clusters,
            n_clusters_ratio=(b.n_clusters / a.n_clusters
                              if a.n_clusters else np.nan),
        )
        if published:
            row.update(P_pub_pre=a.P_published, P_pub_post=b.P_published,
                       P_pub_retained=b.P_published / a.P_published)
        ret.append(row)
    pd.DataFrame(ret).to_csv(path, index=False)


def main():
    os.makedirs(SCRATCH, exist_ok=True)
    print("=" * 78)
    print("PROJECT 2 STEP 0 — Ni + YSZ percolation, full stacks, 6-conn, x axis")
    print("=" * 78)
    test_percolation_lowmem_equivalence()

    gt = ground_truth_frame().set_index("sample")
    ni_rows, ysz_rows = [], []

    for key, folder, grain, state, *_rest in SAMPLES:
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        print(f"\n--- {key} ({grain}, {state}) ---")
        print(f"    labels: Ni={mapping['Ni']} YSZ={mapping['YSZ']} "
              f"pore={mapping['pore']}")

        for phase, rows in (("Ni", ni_rows), ("YSZ", ysz_rows)):
            t0 = time.time()
            mask = load_phase_mask(folder, mapping[phase])
            tload = time.time() - t0
            t0 = time.time()
            r = analyse(mask)
            dt = time.time() - t0
            shape, size = mask.shape, mask.size
            del mask
            gc.collect()

            print(f"  [{phase}] {shape} = {size/1e6:.0f} Mvox  "
                  f"load {tload:.0f}s  analyse {dt:.0f}s")
            print(f"       phi={r['volume_fraction']:.4f} "
                  f"clusters={r['n_clusters']:,}  spans_x={r['percolates']}  "
                  f"P_span={r['P_span']:.4f}  P_reach={r['P_reach']:.4f}  "
                  f"P_largest={r['P_largest']:.4f}")

            base = dict(sample=key, grain=grain, state=state,
                        n_phase_voxels=r["n_phase_voxels"],
                        volume_fraction=r["volume_fraction"],
                        n_clusters=r["n_clusters"],
                        percolates_x=r["percolates"],
                        n_spanning_clusters=r["n_spanning_clusters"],
                        P_span=r["P_span"], P_reach=r["P_reach"],
                        P_largest=r["P_largest"],
                        percolates_z=r["percolates_axis0"],
                        percolates_y=r["percolates_axis1"],
                        seconds=round(dt, 1))
            if phase == "Ni":
                pub_P = float(gt.loc[key, "Ni_P__T-S4"])
                pub_phi = float(gt.loc[key, "Ni_Phi__T-S4"])
                # original phase5 schema: n_ni_voxels, phi/P published columns
                rows.append(dict(
                    sample=key, grain=grain, state=state,
                    phi_ni_published=pub_phi, P_published=pub_P,
                    n_ni_voxels=r["n_phase_voxels"],
                    n_clusters=r["n_clusters"],
                    percolates_x=r["percolates"],
                    n_spanning_clusters=r["n_spanning_clusters"],
                    P_span=r["P_span"], P_reach=r["P_reach"],
                    P_largest=r["P_largest"],
                    percolates_z=r["percolates_axis0"],
                    percolates_y=r["percolates_axis1"],
                    P_reach_over_published=r["P_reach"] / pub_P,
                    seconds=round(dt, 1)))
                print(f"       published Ni P={pub_P:.3f}  "
                      f"P_reach/published={r['P_reach']/pub_P:.3f}")
            else:
                rows.append(base)

            write_all(ni_rows, ysz_rows)

    print("\n" + "=" * 78)
    print("NI  (phase5 schema, regenerated complete)")
    print("=" * 78)
    dn = pd.DataFrame(ni_rows)
    print(dn[["sample", "grain", "state", "n_clusters", "percolates_x",
              "P_span", "P_reach", "P_largest", "P_published"]]
          .to_string(index=False))

    print("\n" + "=" * 78)
    print("YSZ  (Step 0 — never measured before)")
    print("=" * 78)
    dy = pd.DataFrame(ysz_rows)
    print(dy[["sample", "grain", "state", "volume_fraction", "n_clusters",
              "percolates_x", "P_span", "P_reach", "P_largest"]]
          .to_string(index=False))

    for nm, path in (("Ni", os.path.join(OUT5, "phase5_retention.csv")),
                     ("YSZ", os.path.join(OUT2, "step0_ysz_retention.csv"))):
        if os.path.exists(path):
            print(f"\n--- {nm} retention (degraded / pristine) ---")
            print(pd.read_csv(path).to_string(index=False))

    try:
        os.rmdir(SCRATCH)
    except OSError:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
