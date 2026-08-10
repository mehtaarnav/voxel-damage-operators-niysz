"""Stack discovery + streaming I/O for the Zenodo segmented TIFF stacks.

Design note: the volumes are 0.48-1.11 GIGAvoxels and this machine has ~2.6 GB
of free RAM, so anything that can be done slice-by-slice IS done slice-by-slice.
Only explicitly-requested sub-volumes are ever materialised in memory.
"""

from __future__ import annotations

import os
import re
from typing import Iterator

import numpy as np
import tifffile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXTRACTED = os.path.join(HERE, "data", "extracted")


def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def slice_paths(zenodo_folder: str) -> list[str]:
    """Sorted list of TIFF slice paths for e.g. '3_Rx36_Segmented'."""
    root = os.path.join(EXTRACTED, zenodo_folder)
    if not os.path.isdir(root):
        raise FileNotFoundError(root)
    hits = []
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if fn.lower().endswith((".tif", ".tiff")):
                hits.append(os.path.join(dp, fn))
    hits.sort(key=lambda p: _natural_key(os.path.basename(p)))
    if not hits:
        raise FileNotFoundError(f"no TIFFs under {root}")
    return hits


def iter_slices(zenodo_folder: str) -> Iterator[np.ndarray]:
    for p in slice_paths(zenodo_folder):
        yield tifffile.imread(p)


def stack_shape(zenodo_folder: str) -> tuple[int, int, int]:
    """(nz, ny, nx) where nz = number of slice files."""
    ps = slice_paths(zenodo_folder)
    a = tifffile.imread(ps[0])
    return (len(ps), a.shape[0], a.shape[1])


def label_histogram(zenodo_folder: str, max_slices: int | None = None) -> dict:
    """Exact voxel count per label value, streamed (never holds the volume)."""
    counts = np.zeros(256, dtype=np.int64)
    ps = slice_paths(zenodo_folder)
    if max_slices:
        ps = ps[:max_slices]
    dtypes, shapes = set(), set()
    for p in ps:
        a = tifffile.imread(p)
        dtypes.add(str(a.dtype))
        shapes.add(a.shape)
        if a.dtype != np.uint8:
            a = a.astype(np.uint8)
        counts += np.bincount(a.ravel(), minlength=256)
    return {
        "counts": counts,
        "n_slices": len(ps),
        "dtypes": sorted(dtypes),
        "shapes": sorted(shapes),
        "total": int(counts.sum()),
    }


def load_subvolume(zenodo_folder: str, z0: int, z1: int,
                   y0: int, y1: int, x0: int, x1: int) -> np.ndarray:
    """Materialise [z0:z1, y0:y1, x0:x1] as a uint8 array."""
    ps = slice_paths(zenodo_folder)[z0:z1]
    out = np.empty((len(ps), y1 - y0, x1 - x0), dtype=np.uint8)
    for i, p in enumerate(ps):
        a = tifffile.imread(p)
        out[i] = a[y0:y1, x0:x1]
    return out
