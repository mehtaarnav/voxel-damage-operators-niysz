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


# ---------------------------------------------------------------- O5
# Volume-conserving Ni agglomeration. Frozen in
# out/project2/PREREG_AMENDMENT_O5.md (commit cb2ca49) BEFORE any run.
#
# WHY IT EXISTS. Step 2 failed on a magnitude, not an ordering: O1 destroyed
# 99.5% of TPB by the time Ni lost percolation, against a real fine anode that
# RETAINS 79.9% of its TPB at its worst Ni-percolation retention. Every Step 2
# operator removes volume, yet real Phi_Ni RISES in the coarse anode (+0.0146).
# The missing mechanism is redistribution, not removal.
#
# FLOW DIRECTION. In the dewetting/agglomeration regime matter moves from high
# mean curvature to low: a neck of radius r has mean curvature 1/(2r) against a
# particle of radius R with 1/R, so for r < R/2 material leaves the neck and
# joins the particle -- the Rayleigh-type instability that thins necks and grows
# bodies. Local Ni density in a W_CURV window is used as the curvature proxy:
# low density at a surface voxel means locally convex, i.e. high curvature.
#
# NO LARGEST-COMPONENT PRUNING, unlike O1/D4 and deliberately so: O1's pruning
# modelled Ni loss; O5 models Ni redistribution, and deleting disconnected Ni
# would destroy the volume conservation that is the entire point. Disconnected
# Ni stays in the mask and is simply not counted by P_span.
O5_MOVE_FRAC, O5_W_CURV = 0.05, 5


def apply_o5(ni_mask, ysz_mask, n_rounds, seed,
             move_frac=O5_MOVE_FRAC, w_curv=O5_W_CURV):
    """O5 -- volume-conserving, curvature-driven Ni agglomeration."""
    ni = ni_mask.copy()
    n0 = int(ni.sum())
    rng = np.random.default_rng(seed)
    moved_total = 0
    for _ in range(max(0, int(n_rounds))):
        non_ni = ~ni
        surface = ni & ndi.binary_dilation(non_ni, STRUCT6)
        pore = non_ni & ~ysz_mask
        front = pore & ndi.binary_dilation(ni, STRUCT6)
        n_surf = int(surface.sum())
        if n_surf == 0:
            break
        k = int(round(move_frac * n_surf))
        if k <= 0:
            break
        dens = ndi.uniform_filter(ni.astype(np.float32), size=w_curv)
        jit = rng.random(ni.shape).astype(np.float32) * 1e-3

        si = np.flatnonzero(surface.ravel())
        sv = dens.ravel()[si] + jit.ravel()[si]
        k_rm = min(k, si.size)
        rm = si[np.argpartition(sv, k_rm - 1)[:k_rm]]        # lowest density

        fi = np.flatnonzero(front.ravel())
        if fi.size == 0:
            break
        fv = dens.ravel()[fi] + jit.ravel()[fi]
        k_add = min(k_rm, fi.size)
        add = fi[np.argpartition(-fv, k_add - 1)[:k_add]]    # highest density
        rm = rm[:k_add]

        flat = ni.ravel()
        flat[rm] = False
        flat[add] = True
        moved_total += k_add
    n1 = int(ni.sum())
    return ni, dict(voxels_pre=n0, voxels_post=n1,
                    voxels_moved=moved_total,
                    volume_error=abs(n1 - n0) / max(n0, 1))


def ni_ysz_interface_area_vox(ni_mask, ysz_mask):
    """Count of Ni/YSZ shared voxel faces -- the perimeter TPB lives on."""
    n = 0
    for ax in range(3):
        for sh in (1, -1):
            n += int((ni_mask & np.roll(ysz_mask, sh, axis=ax)).sum())
    return n


# ---------------------------------------------------------------- O6
def apply_o6(ni_mask, ysz_mask, n_rounds, seed, p_erode=O1_P_ERODE):
    """O6 -- reduction-only Ni surface erosion. O1 with NO oxidative expansion.

    Frozen in out/project2/PREREG_O6.md (commit f980562) before implementation.

    WHY THE EXPANSION STEP IS GONE. O1/D4 dilates Ni by one voxel into pore,
    which models Ni->NiO oxidation (~70% volume expansion; NiO ~11.2 vs Ni
    ~6.59 cm^3/mol). The process modelled here is Ni loss under REDUCING
    conditions -- coarsening, dissolution, electrochemical removal -- which is
    volume-losing or volume-conserving, not expanding. On the synthetic platform
    the step was benign only because gate G1-c enforced pristine P_span=1.000,
    so there was nothing to heal. Measured on real voxels, it raised pristine
    P_span from 0.9821/0.9713/0.8878 to exactly 1.0000 and multiplied TPB by
    7.7-15.2x. That is an operator artifact, not a degradation mechanism.

    O1 in cmlib/damage.py is untouched and is not run on real data.
    """
    if n_rounds <= 0:
        return ni_mask.copy(), dict(voxels_pre=int(ni_mask.sum()),
                                    voxels_post=int(ni_mask.sum()),
                                    voxels_removed_erosion=0)
    rng = np.random.default_rng(seed)
    cur = ni_mask.copy()
    removed = 0
    for _ in range(int(n_rounds)):
        eroded = ndi.binary_erosion(cur, structure=STRUCT6)
        boundary = cur & ~eroded
        rm = boundary & (rng.random(cur.shape) < p_erode)
        removed += int(rm.sum())
        cur &= ~rm
    lab, n = ndi.label(cur, structure=STRUCT6)
    if n == 0:
        final = np.zeros_like(cur)
    else:
        c = np.bincount(lab.ravel())
        c[0] = 0
        final = lab == int(np.argmax(c))
    return final, dict(voxels_pre=int(ni_mask.sum()),
                       voxels_post=int(final.sum()),
                       voxels_removed_erosion=removed,
                       voxels_removed_islands=int(cur.sum() - final.sum()))


# ---------------------------------------------------------------- O5v2
# True volume-conserving curvature-driven agglomeration. Frozen in
# out/project2/PREREG_O5V2.md (commit 304aa8a) before implementation.
# NOT erosion (no p_erode) and NOT O5's density proxy, which roughened instead
# of agglomerating. Matter leaves convex surface (high curvature, high chemical
# potential) and joins concave necks -- the direction that reduces surface area
# at fixed volume.
O5V2_P_COARSEN = 0.03


def _mean_curvature(ni, sites):
    """k = (#pore 6-neighbours) - (#Ni 6-neighbours); + convex, - concave."""
    nb = ndi.convolve(ni.astype(np.int8), STRUCT6.astype(np.int8),
                      mode="constant", cval=0)
    return (6 - 2 * nb)[sites]


def apply_o5v2(ni_mask, ysz_mask, n_rounds, seed,
               p_coarsen=O5V2_P_COARSEN, conn26=False):
    """Remove the most-convex surface voxels, add at the most-concave neck
    sites, exactly conserving Ni volume. `conn26` switches the curvature
    stencil to 26-connectivity -- the pre-registered implementation fix if the
    6-connectivity proxy fails the surface-area-reduction gate."""
    st = ndi.generate_binary_structure(3, 3) if conn26 else STRUCT6
    ni = ni_mask.copy()
    rng = np.random.default_rng(seed)
    n_surf0 = int((ni & ~ndi.binary_erosion(ni, structure=STRUCT6)).sum())
    k_move = int(round(p_coarsen * n_surf0))
    moved = 0
    for _ in range(max(0, int(n_rounds))):
        if k_move <= 0:
            break
        nb = ndi.convolve(ni.astype(np.int16), st.astype(np.int16),
                          mode="constant", cval=0)
        nmax = int(st.sum())
        surf = ni & ~ndi.binary_erosion(ni, structure=STRUCT6)
        pore = (~ni) & (~ysz_mask)
        front = pore & ndi.binary_dilation(ni, structure=STRUCT6)
        si = np.flatnonzero(surf.ravel())
        fi = np.flatnonzero(front.ravel())
        if si.size == 0 or fi.size == 0:
            break
        jit = rng.random(ni.size).astype(np.float32) * 1e-3
        # convexity of a Ni voxel = few Ni neighbours -> nmax - 2*nb large
        conv = (nmax - 2.0 * nb.ravel()[si]) + jit[si]
        # concavity of a pore site = many Ni neighbours -> nb large
        conc = nb.ravel()[fi] + jit[fi]
        k = min(k_move, si.size, fi.size)
        rm = si[np.argpartition(-conv, k - 1)[:k]]
        add = fi[np.argpartition(-conc, k - 1)[:k]]
        flat = ni.ravel()
        flat[rm] = False
        flat[add] = True
        moved += k
    n0, n1 = int(ni_mask.sum()), int(ni.sum())
    return ni, dict(voxels_pre=n0, voxels_post=n1, voxels_moved=moved,
                    volume_error=abs(n1 - n0) / max(n0, 1), k_move=k_move)


# ------------------------------------------------------- O5v2 Option B
# Greedy (zero-temperature) KMC agglomeration. Frozen in
# out/project2/PREREG_O5V2_OPTIONB.md (c414e63) before implementation.
#
# EXACT dA, derived rather than proxied -- this is what Option A got wrong.
# Removing surface voxel a changes exposed Ni faces by 2*nb(a) - 6; adding at
# pore site b by 6 - 2*nb(b), where nb counts Ni 6-neighbours. Hence
#     dA = 2*(nb(a) - nb(b))       and   dA <= 0  <=>  nb(a) <= nb(b).
# Curvature RANK was a proxy for this; the identity is exact. dV = 0 by
# construction (one voxel out, one in). gamma=1.0, lambda enforced structurally,
# kT=0 so no dA>0 move is ever accepted -- gate (ii) cannot fail by construction
# up to the non-adjacency condition, which is enforced below.
O5V2B_GAMMA, O5V2B_KT = 1.0, 0.0


def apply_o5v2b(ni_mask, ysz_mask, n_rounds, seed,
                p_coarsen=O5V2_P_COARSEN):
    ni = ni_mask.copy()
    rng = np.random.default_rng(seed)
    n_surf0 = int((ni & ~ndi.binary_erosion(ni, structure=STRUCT6)).sum())
    k_move = int(round(p_coarsen * n_surf0))
    proposed = accepted = 0
    for _ in range(max(0, int(n_rounds))):
        if k_move <= 0:
            break
        nb = ndi.convolve(ni.astype(np.int16), STRUCT6.astype(np.int16),
                          mode="constant", cval=0)
        surf = ni & ~ndi.binary_erosion(ni, structure=STRUCT6)
        front = (~ni) & (~ysz_mask) & ndi.binary_dilation(ni, structure=STRUCT6)
        si = np.flatnonzero(surf.ravel())
        fi = np.flatnonzero(front.ravel())
        if si.size == 0 or fi.size == 0:
            break
        nba, nbb = nb.ravel()[si], nb.ravel()[fi]
        a_ord = si[np.argsort(nba, kind="stable")]          # fewest Ni nbrs
        b_ord = fi[np.argsort(-nbb, kind="stable")]         # most Ni nbrs
        nba_s = np.sort(nba, kind="stable")
        nbb_s = -np.sort(-nbb, kind="stable")
        k = min(k_move, a_ord.size, b_ord.size)
        flat = ni.ravel()
        shp = ni.shape
        moved_this = 0
        for t in range(k):
            proposed += 1
            if nba_s[t] > nbb_s[t]:      # dA > 0 -> reject (and all later)
                break
            ia, ib = int(a_ord[t]), int(b_ord[t])
            za, ya, xa = np.unravel_index(ia, shp)
            zb, yb, xb = np.unravel_index(ib, shp)
            if abs(int(za) - int(zb)) + abs(int(ya) - int(yb)) + \
               abs(int(xa) - int(xb)) <= 1:
                continue                  # adjacency breaks the dA algebra
            flat[ia] = False
            flat[ib] = True
            accepted += 1
            moved_this += 1
        if moved_this == 0:
            break
    n0, n1 = int(ni_mask.sum()), int(ni.sum())
    return ni, dict(voxels_pre=n0, voxels_post=n1,
                    volume_error=abs(n1 - n0) / max(n0, 1),
                    volume_delta_vox=n1 - n0, k_move=k_move,
                    proposed=proposed, accepted=accepted,
                    acceptance_rate=accepted / max(proposed, 1),
                    gamma=O5V2B_GAMMA, kT=O5V2B_KT)
