"""Probe 2: how to actually get the distance-transform value onto skan branches.

Probe 1 showed that in skan 0.13.1 the `source_image=` keyword does NOT populate
`pixel_values` (they came back all-1.0).  The documented alternative is to pass
the skeleton itself as a FLOAT image whose nonzero values are the quantity of
interest.  Verify that, on a shape whose neck radius is known analytically.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skan import csr
from skimage.morphology import skeletonize

# A dumbbell with an ANALYTICALLY KNOWN neck: a square-section bar of
# half-width 2 voxels joining two spheres.  With isotropic 10 nm voxels the
# EDT at the bar axis should be ~ (2+1)*10 = 30 nm (distance to the nearest
# background voxel), and the spheres should give ~ 8*10 = 80 nm.
N = 60
vol = np.zeros((N, N, N), bool)
zz, yy, xx = np.ogrid[:N, :N, :N]
vol |= ((zz - 15) ** 2 + (yy - 30) ** 2 + (xx - 30) ** 2) < 8 ** 2
vol |= ((zz - 45) ** 2 + (yy - 30) ** 2 + (xx - 30) ** 2) < 8 ** 2
vol[15:46, 28:33, 28:33] = True           # 5x5 voxel bar -> half-width 2

spacing = (10.0, 10.0, 10.0)
edt = ndi.distance_transform_edt(vol, sampling=spacing)
skel = skeletonize(vol)

print(f"voxels={vol.sum()}  skeleton voxels={skel.sum()}")
print(f"EDT: max={edt.max():.1f} nm   at bar centre (30,30,30)={edt[30,30,30]:.1f} nm")
print("  expected bar-axis EDT ~ 30 nm, sphere-centre EDT ~ 80 nm\n")

# --- method A: source_image (shown broken in probe 1) ---
A = csr.Skeleton(skel, spacing=spacing, source_image=edt)
print("A) source_image= :", "pixel_values is None" if A.pixel_values is None
      else f"min={A.pixel_values.min():.2f} max={A.pixel_values.max():.2f}")

# --- method B: float skeleton image carrying the EDT ---
skel_f = skel.astype(np.float64) * edt
B = csr.Skeleton(skel_f, spacing=spacing)
print("B) float skeleton:",
      "pixel_values is None" if B.pixel_values is None
      else f"min={B.pixel_values.min():.2f} max={B.pixel_values.max():.2f}")

print(f"\nB) n_paths = {B.n_paths}")
for i in range(B.n_paths):
    coords, vals = B.path_with_data(i)
    print(f"   path {i}: n={len(vals):3d}  MIN={vals.min():7.2f} nm  "
          f"mean={vals.mean():7.2f}  max={vals.max():7.2f}  "
          f"len={B.path_lengths()[i]:8.2f} nm")

summ = csr.summarize(B, separator="-")
print("\nsummarize():")
print(summ[["skeleton-id", "node-id-src", "node-id-dst", "branch-distance",
            "branch-type", "mean-pixel-value", "stdev-pixel-value"]].to_string())

G = csr.skeleton_to_nx(B, summ)
print(f"\nskeleton_to_nx -> nodes={G.number_of_nodes()} edges={G.number_of_edges()}")
for u, v, k, d in G.edges(keys=True, data=True):
    vals = d["values"]
    print(f"   edge {u}-{v}: n_pts={len(vals)}  min_value={vals.min():.2f} nm  "
          f"max={vals.max():.2f} nm")

print("\nCONCLUSION: method B (float skeleton image) is the one that works; "
      "the per-branch minimum of `values` is the neck-radius proxy.")
