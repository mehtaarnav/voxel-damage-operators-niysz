"""Probe: how much RAM can we actually get, and how big a volume can the
Phase-3 pipeline (EDT + skeletonize + label) handle?

This decides the ROI size, which is a load-bearing methodological choice.
"""
from __future__ import annotations

import time

import numpy as np
from scipy import ndimage as ndi

try:
    import psutil
    vm = psutil.virtual_memory()
    print(f"psutil: total {vm.total/1e9:.1f} GB   available {vm.available/1e9:.1f} GB")
except ImportError:
    print("psutil not installed")

# how large a float64 block can we actually allocate and touch?
print("\nAllocation probe (allocate AND write, so it is really committed):")
biggest = 0
for gb in (1, 2, 3, 4, 6, 8):
    try:
        a = np.empty(int(gb * 1e9 / 8), dtype=np.float64)
        a[::1000] = 1.0
        del a
        print(f"  {gb} GB float64 : OK")
        biggest = gb
    except (MemoryError, np.core._exceptions._ArrayMemoryError) as e:
        print(f"  {gb} GB float64 : FAILED ({type(e).__name__})")
        break
print(f"  -> largest successful single allocation: {biggest} GB")

# realistic pipeline probe on a modest cube
print("\nPipeline timing/memory probe on a synthetic 300^3 two-phase volume:")
rng = np.random.default_rng(0)
v = ndi.gaussian_filter(rng.random((300, 300, 300), dtype=np.float32), 4) > 0.5
print(f"  mask                     {v.nbytes/1e6:8.1f} MB  frac={v.mean():.3f}")
t = time.time()
d = ndi.distance_transform_edt(v, sampling=(1.0, 1.0, 1.0))
print(f"  distance_transform_edt   {d.nbytes/1e6:8.1f} MB  dtype={d.dtype}  "
      f"{time.time()-t:.1f} s")
t = time.time()
lab, n = ndi.label(v, structure=ndi.generate_binary_structure(3, 1))
print(f"  label                    {lab.nbytes/1e6:8.1f} MB  dtype={lab.dtype}  "
      f"n={n}  {time.time()-t:.1f} s")

n_vox = 300 ** 3
print(f"\n  per-voxel cost of just EDT(float64)+label(int32): "
      f"{(d.nbytes + lab.nbytes)/n_vox:.1f} bytes/voxel")
print("  (plus the mask, the skeleton, and skan's internal graph)")
