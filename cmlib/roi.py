"""Non-overlapping ROI tiling in PHYSICAL space.

Why physical and not voxel space: the six stacks have voxel sizes from 17.9 to
29.14 nm.  Comparing graph metrics (algebraic connectivity, min-cut) between
sub-volumes of equal VOXEL count but unequal PHYSICAL size would confound the
comparison, because those metrics are extensive in domain size.  Every ROI is
therefore the same physical cube, and differs in voxel count between samples.
"""

from __future__ import annotations


def tile_rois(nz, ny, nx, vz_nm, vy_nm, vx_nm, side_um, max_rois=None):
    """Tile non-overlapping cubes of physical side `side_um`.

    Returns a list of dicts with voxel bounds and the achieved physical size.
    ROIs are centred within the stack so that leftover margin is split evenly.
    """
    sz = int(round(side_um * 1000.0 / vz_nm))
    sy = int(round(side_um * 1000.0 / vy_nm))
    sx = int(round(side_um * 1000.0 / vx_nm))
    if sz > nz or sy > ny or sx > nx:
        return []

    kz, ky, kx = nz // sz, ny // sy, nx // sx
    oz = (nz - kz * sz) // 2
    oy = (ny - ky * sy) // 2
    ox = (nx - kx * sx) // 2

    out = []
    for iz in range(kz):
        for iy in range(ky):
            for ix in range(kx):
                z0 = oz + iz * sz
                y0 = oy + iy * sy
                x0 = ox + ix * sx
                out.append(dict(
                    roi=f"z{iz}y{iy}x{ix}",
                    z0=z0, z1=z0 + sz, y0=y0, y1=y0 + sy, x0=x0, x1=x0 + sx,
                    nz=sz, ny=sy, nx=sx,
                    nvox=sz * sy * sx,
                    size_um=(sz * vz_nm / 1000, sy * vy_nm / 1000,
                             sx * vx_nm / 1000),
                ))
    if max_rois:
        out = out[:max_rois]
    return out
