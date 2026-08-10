"""
Three-phase boundary (TPB) length density.

NO STANDARD LIBRARY EXISTS FOR THIS, so the convention is implemented here and
stated in full.  Different conventions in the literature give absolute numbers
that differ by tens of percent, so this choice matters for any comparison
against published values.

CONVENTION CHOSEN: VOXEL-EDGE COUNTING ("classical edge length")
---------------------------------------------------------------
Treat voxel (i,j,k) as the cube spanning corners i..i+1, j..j+1, k..k+1.  A TPB
is a line where the Ni, YSZ and pore phases all meet.  In a labelled voxel grid
such a line runs along VOXEL EDGES.  Each interior voxel edge is shared by
exactly four voxels:

    edge parallel to z  <->  4 voxels differing in (y, x), same z   -> length dz
    edge parallel to y  <->  4 voxels differing in (z, x), same y   -> length dy
    edge parallel to x  <->  4 voxels differing in (z, y), same x   -> length dx

An edge is counted as TPB if and only if the set of phases among its four
surrounding voxels contains ALL THREE of Ni, YSZ and pore.  Its contribution to
the total TPB length is the physical voxel dimension along that edge, so
anisotropic voxels are handled correctly.

    TPB density = (total TPB length) / (volume of the analysed region)
    units: um of TPB per um^3  ==  um^-2

WHY 4-VOXEL EDGE CONFIGURATIONS AND NOT A 6- OR 26-CONNECTIVITY RULE
--------------------------------------------------------------------
"Connectivity" in the 6/26 sense describes which voxels belong to the same
cluster; it is not the right notion for a LINE where three phase VOLUMES meet.
The three-phase line is a geometric feature of the voxel tessellation, and the
four-voxel edge configuration is its exact discrete counterpart.  This is the
convention most commonly called "voxel edge counting" or the "classical" method.

KNOWN BIAS, AND WHY WE DO NOT SILENTLY CORRECT IT
--------------------------------------------------
Voxel-edge counting follows a staircase path along the cubic grid rather than
the smooth true TPB line, so it OVER-ESTIMATES TPB length relative to
centroid-based or marching-cubes-based measures.  The two papers used a third
convention again: "the length of each TPB line is determined based on the
skeletonization of TPB-voxels in each object" (ma8105370, Methods).  We do NOT
apply a literature correction factor, because a factor quoted from another
study is provenance rather than justification.  Instead we report the raw
voxel-edge value and the RATIO to the published value for every sample; if the
convention difference is a roughly constant multiplicative factor then that
ratio will be approximately constant across samples, which is directly
testable, and rankings are unaffected by a constant factor.

STREAMING
---------
Only two adjacent z-slices are ever held, so a full 1.1-gigavoxel stack can be
processed exactly without sub-sampling.
"""

from __future__ import annotations

import numpy as np


def _tpb_counts_from_pair(a: np.ndarray, b: np.ndarray, lab: dict):
    """TPB edge counts using slices z=a and z=b (b is a's successor).

    Returns (n_edges_z, n_edges_y, n_edges_x) for the layer between a and b.
      * z-parallel edges are evaluated on slice `a` alone (4 voxels in y,x)
      * y-parallel edges need a and b (4 voxels in z,x)
      * x-parallel edges need a and b (4 voxels in z,y)
    """
    ni, ysz, po = lab["Ni"], lab["YSZ"], lab["pore"]

    def trio(*arrs):
        """True where the four supplied arrays together contain all 3 phases."""
        hn = np.zeros(arrs[0].shape, dtype=bool)
        hy = np.zeros(arrs[0].shape, dtype=bool)
        hp = np.zeros(arrs[0].shape, dtype=bool)
        for q in arrs:
            hn |= (q == ni)
            hy |= (q == ysz)
            hp |= (q == po)
        return hn & hy & hp

    # --- z-parallel edges: 4 voxels differing in (y, x) within slice a -----
    nz_edges = int(trio(a[:-1, :-1], a[:-1, 1:], a[1:, :-1], a[1:, 1:]).sum())

    # --- y-parallel edges: 4 voxels differing in (z, x) -------------------
    ny_edges = int(trio(a[:, :-1], a[:, 1:], b[:, :-1], b[:, 1:]).sum())

    # --- x-parallel edges: 4 voxels differing in (z, y) -------------------
    nx_edges = int(trio(a[:-1, :], a[1:, :], b[:-1, :], b[1:, :]).sum())

    return nz_edges, ny_edges, nx_edges


def tpb_density_streaming(slice_iter, lab: dict, spacing_nm,
                          shape=None) -> dict:
    """TPB length density from an iterator of 2D label slices.

    spacing_nm : (dz, dy, dx) in nm
    Returns dict with total length (um), volume (um^3) and density (um^-2).
    """
    dz, dy, dx = (float(s) for s in spacing_nm)
    tot_z = tot_y = tot_x = 0
    prev = None
    nz = 0
    ny = nx = None
    for cur in slice_iter:
        nz += 1
        if ny is None:
            ny, nx = cur.shape
        if prev is not None:
            ez, ey, ex = _tpb_counts_from_pair(prev, cur, lab)
            tot_z += ez
            tot_y += ey
            tot_x += ex
        prev = cur
    # the last slice still contributes its z-parallel edges
    if prev is not None:
        ni, ysz, po = lab["Ni"], lab["YSZ"], lab["pore"]

        def trio(*arrs):
            hn = np.zeros(arrs[0].shape, dtype=bool)
            hy = np.zeros(arrs[0].shape, dtype=bool)
            hp = np.zeros(arrs[0].shape, dtype=bool)
            for q in arrs:
                hn |= (q == ni); hy |= (q == ysz); hp |= (q == po)
            return hn & hy & hp
        tot_z += int(trio(prev[:-1, :-1], prev[:-1, 1:],
                          prev[1:, :-1], prev[1:, 1:]).sum())

    length_nm = tot_z * dz + tot_y * dy + tot_x * dx
    vol_nm3 = float(nz) * float(ny) * float(nx) * dz * dy * dx
    length_um = length_nm / 1e3
    vol_um3 = vol_nm3 / 1e9
    return {
        "tpb_edges_z": tot_z, "tpb_edges_y": tot_y, "tpb_edges_x": tot_x,
        "tpb_edges_total": tot_z + tot_y + tot_x,
        "tpb_length_um": length_um,
        "volume_um3": vol_um3,
        "tpb_density_um-2": length_um / vol_um3 if vol_um3 else float("nan"),
        "nz": nz, "ny": ny, "nx": nx,
    }


def tpb_density_volume(vol: np.ndarray, lab: dict, spacing_nm) -> dict:
    """Same measure, for an in-memory 3D label volume."""
    return tpb_density_streaming((vol[i] for i in range(vol.shape[0])),
                                 lab, spacing_nm)
