"""
Project 2 damage operators O1 / O2 / O3, plus a corrected TPB estimator.

SEPARATE MODULE BY REQUIREMENT. `cmlib/damage.py` and `cmlib/synth.py` are
frozen by Project 1's pre-registration and are neither modified nor deleted
here. O1 deliberately *calls* the frozen `apply_d4` rather than reimplementing
it, so the Ni surface-erosion operator is bit-identical to Project 1's.

MONOTONICITY. The bisection protocol requires damage to be monotone in
`n_rounds`: a contact severed at n must stay severed at n+1. For O2 and O3 the
per-round Bernoulli process is therefore realised as a single uniform draw per
contact, with the contact surviving n rounds iff u < (1-p)^n. This has the exact
marginal survival probability of n independent rounds, is monotone in n by
construction, and makes an evaluation O(1) in n instead of O(n) -- which is what
makes 180 bisections affordable. O1 cannot be collapsed this way (erosion
changes the surface each round) and is iterated.

FROZEN PARAMETERS (Step 2 main arm, set by the advisor, never tuned):
    O1  p_erode = 0.35, expand_vox = 1        (inherited from Project 1)
    O2  p_sever = 0.25, neck_percentile_cut = 25
    O3  p_fracture = 0.25 per round
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from .damage import apply_d4
from .project2 import rasterize_ysz

STRUCT6 = ndi.generate_binary_structure(3, 1)

O1_P_ERODE, O1_EXPAND = 0.35, 1
O2_P_SEVER, O2_PCT = 0.25, 25
O3_P_FRACTURE = 0.25


def survival_mask(n_items, n_rounds, p_per_round, seed):
    """Monotone realisation of `n_rounds` independent Bernoulli(p) trials.

    Returns a boolean array, True = item has SURVIVED all n rounds. Monotone
    non-increasing in `n_rounds` for a fixed seed.
    """
    u = np.random.default_rng(seed).random(n_items)
    return u < (1.0 - p_per_round) ** n_rounds


# ---------------------------------------------------------------- O1
def apply_o1(ni_mask, ysz_mask, n_rounds, seed,
             p_erode=O1_P_ERODE, expand_vox=O1_EXPAND):
    """O1 -- Ni surface erosion / coarsening proxy.

    Thin wrapper over the FROZEN `cmlib.damage.apply_d4`, so this operator is
    identical to the one Project 1 characterised (collapse is surface-area
    mediated; failure voxels are one-voxel-deep surface material). YSZ is not
    touched.
    """
    if n_rounds <= 0:
        return ni_mask.copy(), dict(voxels_pre=int(ni_mask.sum()),
                                    voxels_post=int(ni_mask.sum()))
    return apply_d4(ni_mask, ysz_mask, n_rounds, p_erode, expand_vox, seed)


# ---------------------------------------------------------------- O2
def o2_candidates(throat_conns, throat_diam_nm, pct=O2_PCT):
    """Throats in the lower `pct` percentile of inscribed diameter."""
    if throat_diam_nm.size == 0:
        return np.zeros(0, dtype=bool), np.nan
    thr = float(np.percentile(throat_diam_nm, pct))
    return throat_diam_nm < thr, thr


def apply_o2(ni_mask, regions, throat_conns, throat_diam_nm, n_rounds, seed,
             p_sever=O2_P_SEVER, pct=O2_PCT, region_slices=None):
    """O2 -- Ni contact/neck failure.

    Severs narrow SNOW throats: the watershed interface between the two regions
    is set to pore and dilated by one voxel so the two regions are guaranteed
    6-disconnected. Only throats in the lower `pct` percentile are candidates,
    each severed with probability `p_sever` per round. The largest remaining
    component is kept (isolated Ni is electrically dead), matching O1/D4.
    """
    cand, thr = o2_candidates(throat_conns, throat_diam_nm, pct)
    ni = ni_mask.copy()
    n_sev = 0
    if n_rounds > 0 and cand.any():
        alive = survival_mask(int(cand.sum()), n_rounds, p_sever, seed)
        sev_idx = np.flatnonzero(cand)[~alive]
        n_sev = int(sev_idx.size)
        if n_sev:
            # BOUNDING-BOX LOCALISATION. Doing `regions == label` plus two
            # dilations over the whole domain for every severed throat costs
            # several full-domain passes each; with ~100 severed throats per
            # evaluation and ~180 bisections that would dominate the run. Each
            # throat's interface lies inside the union of its two regions'
            # bounding boxes, so the work is confined to that slab.
            if region_slices is None:
                region_slices = ndi.find_objects(regions.astype(np.int32))
            for t in sev_idx:
                a, b = int(throat_conns[t, 0]) + 1, int(throat_conns[t, 1]) + 1
                sa = region_slices[a - 1] if a - 1 < len(region_slices) else None
                sb = region_slices[b - 1] if b - 1 < len(region_slices) else None
                if sa is None or sb is None:
                    continue
                box = tuple(slice(max(0, min(sa[d].start, sb[d].start) - 2),
                                  min(regions.shape[d],
                                      max(sa[d].stop, sb[d].stop) + 2))
                            for d in range(3))
                sub = regions[box]
                ra, rb = sub == a, sub == b
                if not ra.any() or not rb.any():
                    continue
                iface = (ra & ndi.binary_dilation(rb, STRUCT6)) | (
                    rb & ndi.binary_dilation(ra, STRUCT6))
                ni[box] &= ~ndi.binary_dilation(iface, STRUCT6)
    lab, n = ndi.label(ni, structure=STRUCT6)
    if n == 0:
        final = np.zeros_like(ni)
    else:
        c = np.bincount(lab.ravel())
        c[0] = 0
        final = lab == int(np.argmax(c))
    return final, dict(n_candidates=int(cand.sum()), n_severed=n_sev,
                       neck_threshold_nm=thr,
                       voxels_pre=int(ni_mask.sum()),
                       voxels_post=int(final.sum()))


# ---------------------------------------------------------------- O3
def apply_o3(ysz_centres, ysz_pairs, sintered, r_ysz, w_ysz, shape, ni_mask,
             n_rounds, seed, p_fracture=O3_P_FRACTURE):
    """O3 -- YSZ mechanical fracture of explicit sintered contacts.

    Each intact sintered contact fractures with probability `p_fracture` per
    round. A fractured contact has its NECK VOLUME REMOVED (the mask is rebuilt
    from grains + surviving necks) and is struck from the contact graph, which
    is returned so downstream analysis sees the updated topology.

    NO largest-component pruning is applied to YSZ: disconnected YSZ is still
    physically present and still locally ionically relevant, unlike electrically
    dead Ni. This asymmetry is deliberate and pre-registered (DESIGN_MEMO 3, O3a).
    Ni is not touched.
    """
    sint = np.asarray(sintered, dtype=bool)
    idx = np.flatnonzero(sint)
    intact = sint.copy()
    n_frac = 0
    if n_rounds > 0 and idx.size:
        alive = survival_mask(idx.size, n_rounds, p_fracture, seed)
        broken = idx[~alive]
        n_frac = int(broken.size)
        intact[broken] = False
    ysz = rasterize_ysz(ysz_centres, ysz_pairs, intact, r_ysz, w_ysz, shape)
    ysz &= ~ni_mask
    return ysz, intact, dict(n_sintered_pre=int(sint.sum()),
                             n_fractured=n_frac,
                             n_intact_post=int(intact.sum()))


# ---------------------------------------------------------------- TPB
def tpb_density_um2(ni_mask, ysz_mask, voxel_nm=20.0):
    """Three-phase-boundary line density, in um^-2 (= um of line per um^3).

    Counts voxel-EDGE sites at which all three phases appear among the four
    voxels sharing that edge, for edges along each axis, times the edge length.

    TWO THINGS THAT WERE WRONG IN AN EARLIER SCRATCH VERSION, both fixed and
    unit-tested here:

    1. UNITS. The raw result is nm/nm^3 = 1/nm^2, and 1/nm^2 = 1e6 um^-2.
       (Real anode range, for reference: 1.07-2.65 um^-2.)
    2. PERIODIC WRAP. Using np.roll to gather the four voxels around an edge
       wraps at the domain faces and manufactures triple lines that do not
       exist -- on the analytic single-line test case it over-counted by
       exactly 4x, because the torus closes three extra junctions. Our domains
       are FREE-boundary everywhere else (see cmlib.percolation), so edges are
       now gathered by SLICING and only interior edges are counted.
    """
    pore = ~ni_mask & ~ysz_mask
    shape = ni_mask.shape
    total = 0
    for ax in range(3):
        a1, a2 = [x for x in range(3) if x != ax]
        n1, n2 = shape[a1], shape[a2]
        if n1 < 2 or n2 < 2:
            continue
        acc = None
        for phase in (ni_mask, ysz_mask, pore):
            has = None
            for o1 in (0, 1):
                for o2 in (0, 1):
                    sl = [slice(None)] * 3
                    sl[a1] = slice(o1, n1 - 1 + o1)
                    sl[a2] = slice(o2, n2 - 1 + o2)
                    v = phase[tuple(sl)]
                    has = v.copy() if has is None else (has | v)
            acc = has if acc is None else (acc & has)
        total += int(acc.sum())
    length_nm = total * voxel_nm
    volume_nm3 = float(np.prod(shape)) * voxel_nm ** 3
    return length_nm / volume_nm3 * 1e6
