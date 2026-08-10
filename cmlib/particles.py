"""
Ni particle sizing: watershed (SNOW-style) and c-PSD (bulge-size) measures.

`watershed_particles` / `size_stats` are moved here UNCHANGED in definition
from `phase4d_particles.py`, which validated them on the real Ni-YSZ stacks
(overlay gate passed; ordering fine < medium < coarse survived a
`min_distance` in {2,3,4,6,8} sensitivity sweep — see out/phase4/). They are
the PRIMARY particle-size measure throughout this project.

`cpsd_r50max` is NEW, added for the synthetic decoupling study (R1 in
out/next/EXECUTION_SPEC.md): the watershed measure is, by construction,
sensitive to neck width (widening a neck merges watershed regions, inflating
measured particle size even when the underlying particle geometry is
unchanged). c-PSD ("continuous pore size distribution", here applied to the Ni
solid phase rather than the pore phase -- the method is phase-agnostic) instead
measures, at every voxel, the radius of the largest inscribed sphere covering
that voxel. For a voxel deep inside a bulge, this is set by the bulge geometry
and is only weakly perturbed by a distant, narrow neck, so its 50th-percentile
(r50,max in the Pecho et al. terminology -- see references/Table S4 in
cmlib/ground_truth.py) is a neck-INSENSITIVE size measure and the natural
counterpart to report alongside the watershed measure when the study's whole
point is to move neck width independently of "particle size".

CONVENTIONS
-----------
  * watershed markers: sigma=0.4 voxel Gaussian blur on the physical-unit
    distance transform; peak_local_max with a stated min_distance (voxels);
    watershed on the negated blurred map, mask=Ni. Same as phase4d.
  * particle size (watershed): equivalent-sphere diameter from voxel counts,
    reported as plain mean, median, and volume-weighted mean (the one
    comparable to a laser-diffraction d50). Border-truncated regions excluded
    from size statistics (standard unbiased-stereology choice), counted
    separately.
  * particle size (c-PSD): `porespy.filters.local_thickness`, the standard
    library implementation of the "insert a maximally inscribed sphere at
    every voxel" granulometry method -- used per project convention 4
    (prefer existing libraries for the hard parts) rather than a hand-rolled
    morphological-opening sweep.

POREPY UNIT TRAP (found empirically, documented so it is never silently
re-introduced)
------------------------------------------------------------------------------
`porespy.filters.local_thickness` MUST be called on a distance transform in
VOXEL units (spacing=1), not physical units. Verified on an analytic sphere of
true voxel-radius 10 (`probe_local_thickness.py`): calling it with a
physical-unit (nm) distance transform silently returns wrong, `sizes`-
dependent values (160.4, 20.0, or 179.5 "nm" for `sizes` = 25, None, and 50
respectively, for a sphere whose true radius is 200 nm) with NO error or
warning. Calling it with a voxel-unit distance transform recovers the exact
answer (10.0 voxels) for every `sizes` value tested except `sizes=None`, which
is ALSO broken (returns a uniform 1.0 regardless of true thickness) and must
never be used. The fix used throughout this module: always compute the
distance transform at unit voxel spacing, run `local_thickness` with an
explicit integer `sizes`, and rescale the OUTPUT to physical units afterward.
This rescaling is only valid for (near-)isotropic voxels, consistent with
Q7 in the execution spec (synthetic study runs on isotropic voxels).
"""

from __future__ import annotations

import numpy as np
import porespy as ps
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.segmentation import watershed

SIGMA_VOX_DEFAULT = 0.4

# `sizes` must be an int/array; `sizes=None` is buggy (see module docstring).
# RAISED 25 -> 100 (2026-08-10). sizes=25 quantizes the returned c-PSD into
# bins ~6% wide in diameter terms -- comparable to the +-5% c-PSD gate itself,
# which made per-structure deviations unreliable: on the platform-v2
# qualification set, one seed measured -3.13% deviation at sizes=25 (a "pass")
# but -9.18% at sizes=200 (a clear fail), and only 8 distinct c-PSD values
# appeared across 15 structures. Deviations had substantially converged by
# sizes=100. Anything comparing c-PSD against a few-percent tolerance MUST use
# a converged `sizes`; 25 is only safe for coarse/qualitative use.
CPSD_SIZES_DEFAULT = 100


def watershed_particles(mask: np.ndarray, spacing_nm, min_distance: int = 4,
                        sigma_vox: float = SIGMA_VOX_DEFAULT):
    """Watershed-segment `mask` into particles. Returns (labels, edt, n_peaks)."""
    edt = ndi.distance_transform_edt(mask, sampling=spacing_nm)
    edt_s = ndi.gaussian_filter(edt, sigma=sigma_vox)
    coords = peak_local_max(edt_s, min_distance=int(min_distance),
                            labels=mask, exclude_border=False)
    markers = np.zeros(mask.shape, dtype=np.int32)
    if len(coords):
        markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
    labels = watershed(-edt_s, markers, mask=mask)
    return labels, edt, len(coords)


def size_stats(labels: np.ndarray, spacing_nm, exclude_border: bool = True):
    """Equivalent-sphere diameter statistics from a watershed label image.

    Returns (stats_dict, per_region_diameter_array_nm). stats_dict is {} if no
    usable region exists.
    """
    n = int(labels.max())
    if n == 0:
        return {}, np.array([])
    vox_vol_nm3 = float(np.prod(spacing_nm))
    counts = np.bincount(labels.ravel(), minlength=n + 1)[1:]

    border = set()
    for ax in range(3):
        for i in (0, labels.shape[ax] - 1):
            border.update(np.unique(np.take(labels, i, axis=ax)).tolist())
    border.discard(0)
    keep = np.ones(n, dtype=bool)
    if exclude_border:
        for b in border:
            if 1 <= b <= n:
                keep[b - 1] = False

    vols = counts[keep] * vox_vol_nm3
    vols = vols[vols > 0]
    if vols.size == 0:
        return {}, np.array([])
    d = (6.0 * vols / np.pi) ** (1.0 / 3.0)          # nm
    vw = float((d * vols).sum() / vols.sum())        # volume-weighted mean
    return {
        "n_regions_total": n,
        "n_regions_border_excluded": int((~keep).sum()),
        "n_regions_used": int(keep.sum()),
        "d_mean_nm": float(d.mean()),
        "d_median_nm": float(np.median(d)),
        "d_volweighted_nm": vw,
        "d_p90_nm": float(np.percentile(d, 90)),
        "vol_median_nm3": float(np.median(vols)),
    }, d


def cpsd_r50max(mask: np.ndarray, voxel_size_nm: float,
                sizes: int = CPSD_SIZES_DEFAULT) -> dict:
    """Neck-insensitive size measure: 50th percentile of the c-PSD (bulge size).

    `voxel_size_nm` is a SCALAR (isotropic, or the geometric mean of an
    anisotropic spacing -- see cmlib.pnm.geometric_voxel_size_nm for the same
    convention used by SNOW). The distance transform and `local_thickness` call
    are done in voxel units and the result is rescaled to nm afterward -- see
    the module docstring for why this order is mandatory.

    Returns d_cPSD_r50max_nm (DIAMETER, i.e. 2x the r50 radius, for direct
    comparability with `d_volweighted_nm` from `size_stats`), plus the
    r25/r50/r75 radii in nm and the fraction of mask voxels that received a
    nonzero thickness (should be ~1.0; a low value flags a `sizes` resolution
    problem).
    """
    if not mask.any():
        return {}
    dt_vox = ndi.distance_transform_edt(mask)            # voxel units, no spacing
    lt_vox = ps.filters.local_thickness(mask, dt=dt_vox, method="dt", sizes=sizes)
    r_vox = lt_vox[mask]
    r_vox = r_vox[r_vox > 0]
    if r_vox.size == 0:
        return {}
    r_nm = r_vox * float(voxel_size_nm)
    return {
        "cpsd_nonzero_frac": float(r_vox.size / int(mask.sum())),
        "cpsd_r25_nm": float(np.percentile(r_nm, 25)),
        "cpsd_r50_nm": float(np.percentile(r_nm, 50)),
        "cpsd_r75_nm": float(np.percentile(r_nm, 75)),
        "d_cPSD_r50max_nm": float(2.0 * np.percentile(r_nm, 50)),
    }
