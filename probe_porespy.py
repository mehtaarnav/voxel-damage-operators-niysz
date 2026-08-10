"""Probe porespy's SNOW network extraction against a shape with a KNOWN throat.

Geometry: two cubes joined by a square-section bar 6 voxels across, at 10 nm
voxels.  The true throat width is 60 nm; the inscribed-diameter measure should
recover ~60-70 nm (the EDT carries a ~+1 voxel bias, see cmlib/graph.py note 4).
We must confirm:
  * which key holds the throat (neck) width
  * that voxel_size is honoured
  * that the network is a proper 2-node/1-throat graph for this shape
"""
from __future__ import annotations

import numpy as np
import porespy as ps

print("porespy version:", ps.__version__)
print("\nnetworks module:", [n for n in dir(ps.networks) if not n.startswith("_")])

import inspect
print("\nsnow2 signature:")
print(inspect.signature(ps.networks.snow2))

# ---- dumbbell with an analytically known throat -------------------------
N = 60
im = np.zeros((N, N, N), bool)
im[10:26, 20:40, 20:40] = True          # cube A
im[34:50, 20:40, 20:40] = True          # cube B
im[25:35, 27:33, 27:33] = True          # 6x6 voxel bar  -> throat width 60 nm

vox = 10e-9    # porespy works in metres by convention when voxel_size is SI
snow = ps.networks.snow2(im, voxel_size=vox, boundary_width=0)
net = snow.network

print(f"\nregions: {snow.regions.max()} labels")
print("network keys:")
for k in sorted(net.keys()):
    v = net[k]
    if isinstance(v, np.ndarray):
        s = f"shape={v.shape} dtype={v.dtype}"
        if v.size and v.ndim == 1 and np.issubdtype(v.dtype, np.number):
            s += f"  min={v.min():.4g} max={v.max():.4g}"
        print(f"   {k:38s} {s}")
    else:
        print(f"   {k:38s} {type(v)}")

print("\nthroat.conns:\n", net["throat.conns"])
for k in net:
    if "throat" in k and "diameter" in k:
        print(f"{k}: {net[k]}  -> in nm: {net[k]*1e9}")
print("\nEXPECTED throat width ~60 nm (6 voxels at 10 nm), +/- ~1 voxel EDT bias")
