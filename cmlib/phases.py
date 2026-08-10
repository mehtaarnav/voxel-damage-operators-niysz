"""
Phase-label assignment for the six Zenodo segmented stacks.

WHY THIS FILE EXISTS
--------------------
Every stack uses a different label encoding:
    3_Rx36      0 / 100 / 200  (+255)
    4_Rx37      0 / 100 / 200
    5_Rx38      1 /   2 /   3
    6_Rx41-1    0 / 100 / 200  (+255)
    7_Rx41-2   36 / 121 / 194  (+255)
    8_Rx41-3    0 /  76 / 150
and 5_Rx38 additionally has Ni and YSZ SWAPPED in brightness order relative to
the other five stacks.

ASSIGNMENT RULE (stated explicitly; deliberately NOT fitted to the papers)
-------------------------------------------------------------------------
1. The three labels with the largest voxel counts are the three phases.  Any
   other label (in practice only 255, which occupies 0.00-0.25 % of a volume) is
   'unassigned' and is mapped to no phase.
2. Those three labels are sorted ASCENDING BY GRAYSCALE VALUE.
3. They are then named using the brightness ordering stated by the dataset's own
   metadata file (2_3D_Data_Info.xlsx, column AE), which is one of
       "Ni white, YSZ gray, pore black"  ->  ascending = pore, YSZ, Ni
       "Ni gray, YSZ white, pore black"  ->  ascending = pore, Ni, YSZ
   Note that pore = black = lowest value in both cases.

This rule uses ONLY the dataset's own documentation.  It does not consult the
published volume fractions, so comparing the resulting volume fractions against
Table S4 of ma8095265 is an independent check rather than a circular fit.

DENOMINATOR CONVENTION
----------------------
Volume fraction = (voxels of phase) / (ALL voxels in the stack), i.e. the
'unassigned' 255 voxels stay in the denominator.  This matches the papers: for
the fine pristine sample their three volume fractions sum to 0.997, and the 255
class in that stack occupies exactly 0.249 % of the volume.
"""

from __future__ import annotations

import numpy as np

# brightness note -> phase names in ASCENDING label-value order
ORDER_FROM_NOTE = {
    "Ni white, YSZ gray, pore black": ("pore", "YSZ", "Ni"),
    "Ni gray, YSZ white, pore black": ("pore", "Ni", "YSZ"),
}


def assign_labels(counts: np.ndarray, brightness_note: str) -> dict:
    """Map grayscale label values to phase names.

    counts : length-256 voxel-count array (from cmlib.io.label_histogram)
    returns {'pore': value, 'YSZ': value, 'Ni': value, 'unassigned': [values]}
    """
    if brightness_note not in ORDER_FROM_NOTE:
        raise ValueError(f"unrecognised brightness note: {brightness_note!r}")
    present = np.nonzero(counts)[0]
    if len(present) < 3:
        raise ValueError(f"expected >=3 labels, found {list(present)}")

    # three largest by count = the three phases
    top3 = present[np.argsort(counts[present])[::-1][:3]]
    top3 = np.sort(top3)                       # ascending by VALUE
    names = ORDER_FROM_NOTE[brightness_note]

    out = {name: int(v) for name, v in zip(names, top3)}
    out["unassigned"] = [int(v) for v in present if v not in set(top3.tolist())]
    return out


def volume_fractions(counts: np.ndarray, mapping: dict) -> dict:
    total = int(counts.sum())
    return {ph: float(counts[mapping[ph]] / total) for ph in ("Ni", "YSZ", "pore")}
