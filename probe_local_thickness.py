"""
Probe: porespy.filters.local_thickness unit handling.

Finding (see cmlib/particles.py module docstring for the consequence): calling
local_thickness with a distance transform computed in PHYSICAL units silently
returns wrong, `sizes`-dependent answers with no error or warning. Calling it
with a distance transform in VOXEL units (spacing=1) recovers the exact answer,
provided `sizes` is an explicit int/array (`sizes=None` is separately broken).

Analytic case: a solid sphere of exact voxel-radius 10 at isotropic 20 nm
voxels, so the true physical radius is 200 nm.
"""
from __future__ import annotations

import numpy as np
import porespy as ps
from scipy import ndimage as ndi

N = 40
zz, yy, xx = np.ogrid[:N, :N, :N]
mask = ((zz - 20) ** 2 + (yy - 20) ** 2 + (xx - 20) ** 2) < 10 ** 2
spacing_nm = (20.0, 20.0, 20.0)

print("=" * 74)
print("A) local_thickness on a PHYSICAL-unit distance transform (WRONG)")
print("=" * 74)
dt_phys = ndi.distance_transform_edt(mask, sampling=spacing_nm)
print(f"  dt_phys.max() = {dt_phys.max():.1f} nm  (true radius 200 nm)")
for sizes in (25, None, 50):
    lt = ps.filters.local_thickness(mask, dt=dt_phys, method="dt", sizes=sizes)
    v = lt[mask]
    print(f"  sizes={str(sizes):5s}: max={v.max():7.1f} nm  "
          f"median={np.median(v):7.1f} nm   <- WRONG, and silently so")

print("\n" + "=" * 74)
print("B) local_thickness on a VOXEL-unit distance transform (CORRECT)")
print("=" * 74)
dt_vox = ndi.distance_transform_edt(mask)
print(f"  dt_vox.max() = {dt_vox.max():.2f} voxels  (true radius 10 voxels)")
for sizes in (25, None, 50):
    lt = ps.filters.local_thickness(mask, dt=dt_vox, method="dt", sizes=sizes)
    v = lt[mask]
    flag = "BUGGY (sizes=None)" if sizes is None else (
        "PASS" if abs(v.max() - 10.0) < 1e-6 else "FAIL")
    print(f"  sizes={str(sizes):5s}: max={v.max():6.2f} vox  "
          f"median={np.median(v):6.2f} vox   {flag}")
    print(f"           -> rescaled to nm: max={v.max()*20:.1f} nm "
          f"(true 200 nm)")

print("\nCONCLUSION: cmlib.particles.cpsd_r50max always computes the distance")
print("transform at voxel spacing, passes an explicit integer `sizes`, and")
print("rescales the OUTPUT to physical units -- never passes a physical-unit")
print("distance transform into local_thickness, and never uses sizes=None.")
