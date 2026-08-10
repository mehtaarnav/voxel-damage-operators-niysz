"""
Percolation / connected-component utilities for 3D segmented microstructures.

METHODOLOGICAL CONVENTIONS (stated explicitly — these are choices, not defaults):

1. VOXEL ADJACENCY FOR PHASE CONNECTIVITY.
   The default throughout this study is 6-connectivity (face-sharing neighbours
   only) for deciding whether two voxels of the same phase belong to the same
   cluster.  Reasons:
     (a) It is the adjacency rule for which the simple-cubic *site* percolation
         threshold p_c = 0.311608 is defined (Xu et al., Phys. Rev. E 89, 012120
         (2014) give p_c = 0.3116077(2)).  Using 18- or 26-connectivity would
         change the reference threshold to 0.1372 and 0.09755 respectively, so
         the Phase-0 gate is only meaningful with 6-connectivity.
     (b) It is the conservative choice: two voxels touching only at an edge or a
         corner share zero interfacial area, so treating them as electronically
         connected would over-report percolation.  The SOFC microstructure
         literature conventionally uses face connectivity for solid-phase
         transport.
   `connectivity` is exposed as an argument so the choice is always visible at
   the call site, never implicit.

2. SPANNING ("percolating") DEFINITION.
   A cluster percolates along an axis if and only if a single connected
   component contains at least one voxel in the first slice (index 0) and at
   least one voxel in the last slice (index n-1) along that axis.  Free
   (non-periodic) boundaries in all directions.  This is the standard
   "face-to-face spanning" criterion.

3. PERCOLATING FRACTION.
   Fraction of the phase's voxels that belong to spanning clusters (there can
   in principle be more than one; they are unioned).  This is a volume
   fraction of the phase, not of the whole domain, unless stated otherwise.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

# ---------------------------------------------------------------------------
# Reference site-percolation thresholds on the simple cubic lattice, used only
# for validating this module (Phase 0).
#   6  neighbours  (sc,  NN)          p_c = 0.3116077(2)
#   18 neighbours  (sc,  NN+2NN)      p_c = 0.1372(1)
#   26 neighbours  (sc,  NN+2NN+3NN)  p_c = 0.09755(2)
# Source: Xu, Wang, Lv, Deng, Phys. Rev. E 89, 012120 (2014), Table I.
# ---------------------------------------------------------------------------
SC_SITE_THRESHOLDS = {6: 0.3116077, 18: 0.1372, 26: 0.09755}

_CONN_RANK = {6: 1, 18: 2, 26: 3}


def structure_for(connectivity: int) -> np.ndarray:
    """3x3x3 binary structuring element for 6-, 18- or 26-connectivity."""
    if connectivity not in _CONN_RANK:
        raise ValueError(
            f"connectivity must be one of {sorted(_CONN_RANK)}, got {connectivity}"
        )
    return ndi.generate_binary_structure(3, _CONN_RANK[connectivity])


def label_phase(mask: np.ndarray, connectivity: int = 6):
    """Label connected components of a boolean 3D mask.

    Returns (labels, n_features).  Label 0 is background.
    """
    if mask.ndim != 3:
        raise ValueError(f"expected a 3D mask, got shape {mask.shape}")
    return ndi.label(mask, structure=structure_for(connectivity))


def spanning_labels(labels: np.ndarray, axis: int) -> np.ndarray:
    """Labels present on BOTH the first and last slice along `axis`.

    Excludes background (0).
    """
    first = np.take(labels, 0, axis=axis)
    last = np.take(labels, labels.shape[axis] - 1, axis=axis)
    a = np.unique(first)
    b = np.unique(last)
    common = np.intersect1d(a, b, assume_unique=True)
    return common[common != 0]


def percolates(mask: np.ndarray, axis: int = 0, connectivity: int = 6) -> bool:
    """True if the mask has a face-to-face spanning cluster along `axis`."""
    labels, n = label_phase(mask, connectivity)
    if n == 0:
        return False
    return spanning_labels(labels, axis).size > 0


def percolation_report(mask: np.ndarray, axis: int = 0, connectivity: int = 6) -> dict:
    """Full percolation description of one phase mask along one axis.

    Keys
    ----
    n_phase_voxels     : total voxels of the phase
    volume_fraction    : phase voxels / domain voxels
    n_clusters         : number of connected components
    percolates         : bool, face-to-face spanning along `axis`
    n_spanning_clusters: how many components span
    percolating_frac   : spanning voxels / phase voxels  (0.0 if none)
    largest_frac       : largest component / phase voxels
    axis, connectivity : echoed back so results are self-documenting
    """
    labels, n = label_phase(mask, connectivity)
    n_phase = int(mask.sum())
    out = {
        "axis": axis,
        "connectivity": connectivity,
        "n_phase_voxels": n_phase,
        "volume_fraction": n_phase / mask.size,
        "n_clusters": int(n),
        "percolates": False,
        "n_spanning_clusters": 0,
        "percolating_frac": 0.0,
        "largest_frac": 0.0,
    }
    if n == 0 or n_phase == 0:
        return out

    # bincount over labels; index 0 is background
    counts = np.bincount(labels.ravel(), minlength=n + 1)
    out["largest_frac"] = float(counts[1:].max() / n_phase)

    span = spanning_labels(labels, axis)
    if span.size:
        out["percolates"] = True
        out["n_spanning_clusters"] = int(span.size)
        out["percolating_frac"] = float(counts[span].sum() / n_phase)
    return out


def percolating_mask(mask: np.ndarray, axis: int = 0, connectivity: int = 6) -> np.ndarray:
    """Boolean mask of only those voxels in spanning clusters along `axis`."""
    labels, n = label_phase(mask, connectivity)
    if n == 0:
        return np.zeros_like(mask, dtype=bool)
    span = spanning_labels(labels, axis)
    if span.size == 0:
        return np.zeros_like(mask, dtype=bool)
    return np.isin(labels, span)


def random_site_medium(shape, p: float, rng: np.random.Generator) -> np.ndarray:
    """Independent-site random medium: each voxel 'on' with probability p."""
    return rng.random(shape) < p


def percolation_summary(mask: np.ndarray, axis: int = 0, connectivity: int = 6,
                        check_other_axes: bool = True) -> dict:
    """Percolation summary distinguishing P_span from P_reach.

    Moved here (unchanged in definition) from the `analyse()` helper in
    phase5_percolation.py, which computed exactly this on the six real Ni-YSZ
    stacks.  It is a superset of `percolation_report`: that function reports
    only the P_span quantity (there called `percolating_frac`); this adds
    P_reach, which the real-data study needed because the published "percolation
    factor" P (Pecho et al., ma8095265 Sec. 3) is a MIP-PSD reachable-from-a-
    boundary measure, not a two-face-spanning measure, so P_reach is the
    like-for-like quantity to compare against it. Both defintions:

        P_span   = (voxels in cluster(s) touching BOTH faces of `axis`)
                   / (total phase voxels)     -- our strict definition
        P_reach  = (voxels in cluster(s) touching AT LEAST ONE face of `axis`)
                   / (total phase voxels)     -- P_reach >= P_span always

    `percolates` and `n_spanning_clusters` refer to the STRICT (both-face,
    P_span) sense, matching `percolation_report`'s convention.

    check_other_axes=True additionally reports whether the mask spans the other
    two axes (boolean only, not P_span/P_reach for those axes) -- an isotropy
    check on real data, and on synthetic data a diagnostic for anisotropic
    generators.

    Returns keys: axis, connectivity, n_phase_voxels, volume_fraction,
    n_clusters, percolates, n_spanning_clusters, P_span, P_reach, P_largest,
    and percolates_axis{other axis index} for each other axis if requested.
    """
    n_phase = int(mask.sum())
    out = {
        "axis": axis,
        "connectivity": connectivity,
        "n_phase_voxels": n_phase,
        "volume_fraction": n_phase / mask.size if mask.size else 0.0,
        "n_clusters": 0,
        "percolates": False,
        "n_spanning_clusters": 0,
        "P_span": 0.0,
        "P_reach": 0.0,
        "P_largest": 0.0,
    }
    if check_other_axes:
        for ax in range(mask.ndim):
            if ax != axis:
                out[f"percolates_axis{ax}"] = False

    labels, n = label_phase(mask, connectivity)
    out["n_clusters"] = int(n)
    if n == 0 or n_phase == 0:
        return out

    counts = np.bincount(labels.ravel(), minlength=n + 1)
    out["P_largest"] = float(counts[1:].max() / n_phase)

    lo = set(int(v) for v in np.unique(np.take(labels, 0, axis=axis)) if v > 0)
    hi = set(int(v) for v in np.unique(
        np.take(labels, labels.shape[axis] - 1, axis=axis)) if v > 0)
    both = np.array(sorted(lo & hi), dtype=np.int64)
    either = np.array(sorted(lo | hi), dtype=np.int64)

    if both.size:
        out["percolates"] = True
        out["n_spanning_clusters"] = int(both.size)
        out["P_span"] = float(counts[both].sum() / n_phase)
    if either.size:
        out["P_reach"] = float(counts[either].sum() / n_phase)

    if check_other_axes:
        for ax in range(mask.ndim):
            if ax == axis:
                continue
            a = set(int(v) for v in np.unique(np.take(labels, 0, axis=ax)) if v > 0)
            b = set(int(v) for v in np.unique(
                np.take(labels, labels.shape[ax] - 1, axis=ax)) if v > 0)
            out[f"percolates_axis{ax}"] = bool(a & b)

    return out
