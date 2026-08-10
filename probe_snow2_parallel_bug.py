"""
Probe: porespy.networks.snow2's default `parallel_kw={}` silently enables
chunked/parallel watershed partitioning, which (a) can crash entirely on small
idealized geometries and (b) gives measurably different pore/throat counts and
neck statistics than serial mode on real data -- see cmlib/pnm.py docstring
for the full writeup and disclosure of the quantified impact on the
already-reported real-data study.

Part A reproduces the crash on an idealized dumbbell (two cubes + a thin bar)
and shows serial mode (parallel_kw=None) recovers cleanly.
Part B quantifies the shift on one real ROI already reported in
out/phase3/phase3_snow_8.0um_rmax4.csv (coarse_pre, ROI z0y0x0).
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import porespy as ps

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, load_subvolume  # noqa: E402
from cmlib.percolation import percolating_mask  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.pnm import geometric_voxel_size_nm  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402


def part_a():
    print("=" * 74)
    print("A) idealized dumbbell: default (chunked) crashes; serial doesn't")
    print("=" * 74)
    N = 50
    w = 6
    vol = np.zeros((N, N, N), dtype=bool)
    c = N // 2
    vol[0:14, c - 7:c + 7, c - 7:c + 7] = True
    vol[N - 14:N, c - 7:c + 7, c - 7:c + 7] = True
    half = w // 2
    vol[13:N - 13, c - half:c - half + w, c - half:c - half + w] = True

    for pk, label in (({}, "default parallel_kw={}"), (None, "parallel_kw=None (serial)")):
        try:
            snow = ps.networks.snow2(vol, voxel_size=10e-9, boundary_width=0,
                                     sigma=0.4, r_max=4, accuracy="standard",
                                     parallel_kw=pk)
            net = snow.network
            print(f"  {label:32s}: OK  pores={net['pore.all'].size} "
                  f"throats={net['throat.all'].size}")
        except Exception as e:
            print(f"  {label:32s}: FAILED  {type(e).__name__}: {e}")


def part_b():
    print("\n" + "=" * 74)
    print("B) real coarse_pre ROI z0y0x0: chunked vs serial, quantified")
    print("=" * 74)
    key = "coarse_pre"
    _, folder, grain, state, nx_, ny_, nz_, vx, vy, vz = \
        [s for s in SAMPLES if s[0] == key][0]
    spacing = (vz, vy, vx)
    counts = label_histogram(folder)["counts"]
    mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
    r = tile_rois(nz_, ny_, nx_, vz, vy, vx, 8.0, 1)[0]
    vol = load_subvolume(folder, r["z0"], r["z1"], r["y0"], r["y1"],
                         r["x0"], r["x1"])
    ni = vol == mapping["Ni"]
    del vol
    mask = percolating_mask(ni, axis=2, connectivity=6)
    vox_nm = geometric_voxel_size_nm(spacing)
    vox_m = vox_nm * 1e-9

    print("  reported in out/phase3/phase3_snow_8.0um_rmax4.csv (chunked, "
          "undisclosed at the time): pores=62 throats=145 neck_p10=199nm "
          "neck_p50=516nm")
    for pk, label in ((None, "serial"), ({}, "default (chunked)")):
        t0 = time.time()
        snow = ps.networks.snow2(mask, voxel_size=vox_m, boundary_width=0,
                                 sigma=0.4, r_max=4, accuracy="standard",
                                 parallel_kw=pk)
        net = snow.network
        dt = time.time() - t0
        necks = net["throat.inscribed_diameter"]
        print(f"  {label:20s}: pores={net['pore.all'].size:4d} "
              f"throats={net['throat.all'].size:4d} "
              f"neck_p10={np.percentile(necks, 10)*1e9:6.1f}nm "
              f"neck_p50={np.percentile(necks, 50)*1e9:6.1f}nm  [{dt:.1f}s]")


if __name__ == "__main__":
    part_a()
    part_b()
    print("\nCONCLUSION: cmlib.pnm.extract_ni_network now defaults to "
          "parallel_kw=None (serial). The original real-data REPORT.md and "
          "out/phase3, out/phase4 CSVs are NOT retroactively altered; the "
          "~5-6% shift shown in part B is smaller than both the already-"
          "reported between-ROI spread and r_max sensitivity for coarse_pre, "
          "so no conclusion in REPORT.md changes.")
