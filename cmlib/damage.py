"""
Ternary YSZ/pore placement and the D4 redox-like damage operator.

Promoted to library code from `scripts/spike/e0_vertical_slice.py`, where both
were first implemented and verified. Definitions are UNCHANGED; only the
location moved, so platform-v2 work and the E0 spike share one implementation
and cannot silently diverge.

YSZ/PORE PLACEMENT (minimal, not optimised for TPB)
---------------------------------------------------
Ni voxels are never touched. The non-Ni remainder is split into YSZ and pore
by thresholding a Gaussian-smoothed random field (blob-like domains, not
salt-and-pepper, so TPB is not swamped by placement noise) at the percentile
that hits the target YSZ-fraction-OF-THE-REMAINDER. Because it only labels
non-Ni voxels, initial Ni P_span is unchanged BY CONSTRUCTION -- still
verified explicitly at every call site rather than assumed.

D4 DAMAGE (non-circular; NOT a threshold on measured neck width, which would
be D1)
-------------------------------------------------------------------------
1. Fixed oxidative expansion: Ni dilates by `expand_vox`, restricted to
   voxels that were originally PORE (never claims YSZ).
2. Stochastic surface erosion, `n_rounds` rounds: each round, current Ni
   SURFACE voxels are removed independently with probability `p_erode`. This
   is uniform and purely geometric -- it never inspects or selects on the
   measured neck-p10 variable under test. Thin necks are emergently more
   vulnerable because their entire cross-section is surface every round,
   whereas a thick particle only loses its outer shell per round.
3. Keep only the single largest remaining connected Ni component; smaller
   fragments are electrically isolated islands (matching the real
   literature's description of disconnected Ni islands under redox cycling).

Voxels that were Ni and are not Ni afterwards become PORE (Ni loss leaves
porosity). **YSZ is never modified by damage** -- asserted at the call site.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from .synthvol import LABELS

STRUCT6 = ndi.generate_binary_structure(3, 1)
SMOOTH_SIGMA_VOX_DEFAULT = 3.0


def add_ysz_pore(ni_mask: np.ndarray, seed: int,
                 ysz_frac_of_rest: float,
                 smooth_sigma_vox: float = SMOOTH_SIGMA_VOX_DEFAULT):
    """Label the non-Ni remainder as YSZ/pore. Returns (ternary_vol, ysz_mask).

    `ysz_frac_of_rest` is the YSZ share OF THE NON-Ni REMAINDER, i.e.
    phi_YSZ/(phi_YSZ+phi_pore) of the anode being matched -- not phi_YSZ of
    the whole volume, which is not attainable once phi_Ni is fixed by the Ni
    geometry.
    """
    rng = np.random.default_rng(seed)
    field = ndi.gaussian_filter(
        rng.standard_normal(ni_mask.shape).astype(np.float32),
        sigma=smooth_sigma_vox)
    non_ni = ~ni_mask
    thresh = np.percentile(field[non_ni], 100.0 * (1.0 - ysz_frac_of_rest))
    ysz_mask = non_ni & (field >= thresh)
    vol = np.full(ni_mask.shape, LABELS["pore"], dtype=np.uint8)
    vol[ysz_mask] = LABELS["YSZ"]
    vol[ni_mask] = LABELS["Ni"]
    return vol, ysz_mask


def apply_d4(ni_mask: np.ndarray, ysz_mask: np.ndarray, n_rounds: int,
             p_erode: float, expand_vox: int, seed: int):
    """D4 redox-like damage. Returns (damaged_ni_mask, diagnostics)."""
    rng = np.random.default_rng(seed)
    pore0 = ~ni_mask & ~ysz_mask

    dil = ndi.binary_dilation(ni_mask, structure=STRUCT6, iterations=expand_vox)
    cur = ni_mask | (dil & pore0)
    added_expand = int(cur.sum() - ni_mask.sum())

    removed_erosion = 0
    for _ in range(n_rounds):
        eroded = ndi.binary_erosion(cur, structure=STRUCT6)
        boundary = cur & ~eroded
        remove = boundary & (rng.random(cur.shape) < p_erode)
        removed_erosion += int(remove.sum())
        cur = cur & ~remove

    labels, n = ndi.label(cur, structure=STRUCT6)
    if n == 0:
        final = np.zeros_like(cur)
        removed_islands = int(cur.sum())
    else:
        counts = np.bincount(labels.ravel(), minlength=n + 1)
        counts[0] = 0
        final = labels == int(np.argmax(counts))
        removed_islands = int(cur.sum() - final.sum())

    return final, dict(voxels_added_expand=added_expand,
                       voxels_removed_erosion=removed_erosion,
                       voxels_removed_islands=removed_islands,
                       voxels_pre=int(ni_mask.sum()),
                       voxels_post=int(final.sum()),
                       n_components_before_island_removal=int(n))


def rebuild_ternary(ni_mask: np.ndarray, ysz_mask: np.ndarray) -> np.ndarray:
    """Ternary volume from a (damaged) Ni mask + the UNCHANGED YSZ mask.
    Anything neither Ni nor YSZ becomes pore."""
    return np.where(ysz_mask, LABELS["YSZ"],
                    np.where(ni_mask, LABELS["Ni"],
                             LABELS["pore"])).astype(np.uint8)
