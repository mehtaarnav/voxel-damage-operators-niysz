"""Probe the installed skan API so Phase 3 uses it correctly rather than from memory."""
from __future__ import annotations

import inspect

import numpy as np
import skan
from skan import csr

print("skan version:", skan.__version__)
print("\ncsr public names:", [n for n in dir(csr) if not n.startswith("_")])

print("\n--- Skeleton.__init__ ---")
print(inspect.signature(csr.Skeleton.__init__))
print("\nSkeleton methods/attrs:",
      [n for n in dir(csr.Skeleton) if not n.startswith("_")])

for fn in ("summarize", "skeleton_to_nx"):
    if hasattr(csr, fn):
        print(f"\n--- csr.{fn} ---")
        print(inspect.signature(getattr(csr, fn)))
        doc = (getattr(csr, fn).__doc__ or "").strip().splitlines()
        print("\n".join("   " + l for l in doc[:25]))

# tiny 3D worked example: a dumbbell (two blobs joined by a thin neck)
vol = np.zeros((40, 40, 40), bool)
zz, yy, xx = np.ogrid[:40, :40, :40]
vol |= ((zz - 12) ** 2 + (yy - 20) ** 2 + (xx - 20) ** 2) < 7 ** 2
vol |= ((zz - 28) ** 2 + (yy - 20) ** 2 + (xx - 20) ** 2) < 7 ** 2
vol[12:29, 19:22, 19:22] = True          # 3x3 voxel neck

from scipy import ndimage as ndi
from skimage.morphology import skeletonize

spacing = (20.0, 19.53, 19.53)           # nm, deliberately anisotropic
edt = ndi.distance_transform_edt(vol, sampling=spacing)
skel = skeletonize(vol)
print(f"\nvolume voxels={vol.sum()}  skeleton voxels={skel.sum()}")
print(f"EDT max = {edt.max():.1f} nm")

S = csr.Skeleton(skel, spacing=spacing, source_image=edt)
print("\nSkeleton built.")
print("  n_paths:", S.n_paths)
print("  path_lengths():", np.round(S.path_lengths(), 1))
for attr in ("coordinates", "pixel_values", "distances", "graph", "degrees"):
    if hasattr(S, attr):
        v = getattr(S, attr)
        print(f"  {attr}: type={type(v).__name__} "
              f"shape={getattr(v, 'shape', None)}")

print("\n  per-path source-image values (min = narrowest point):")
for i in range(S.n_paths):
    coords, vals = S.path_with_data(i)
    print(f"    path {i}: n={len(vals):3d}  min={vals.min():7.2f}  "
          f"mean={vals.mean():7.2f}  max={vals.max():7.2f}  "
          f"length={S.path_lengths()[i]:7.2f} nm")

summ = csr.summarize(S, separator="-")
print("\n  summarize() columns:", list(summ.columns))
print(summ.to_string())

if hasattr(csr, "skeleton_to_nx"):
    G = csr.skeleton_to_nx(S)
    print(f"\n  skeleton_to_nx -> {type(G).__name__}  "
          f"nodes={G.number_of_nodes()} edges={G.number_of_edges()}")
    for u, v, k, d in list(G.edges(keys=True, data=True))[:6]:
        print(f"    edge {u}-{v} key={k} data={d}")
