"""
Project 2 additions. Kept in a SEPARATE module so that `cmlib/damage.py` and
`cmlib/synth.py` -- frozen by Project 1's pre-registration and required to stay
bit-reproducible -- are not touched.

Contains only the A2.3 YSZ-morphology scaling rule at this stage. The damage
operators O1/O2/O3 are NOT implemented here or anywhere; they are specified in
out/project2/DESIGN_MEMO.md sec 3 and await authorization.
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# A2.3 -- YSZ morphology scaling
#
# WHY. `cmlib.damage.add_ysz_pore` builds YSZ by thresholding a Gaussian-
# smoothed random field. Its `smooth_sigma_vox` already exists as a parameter,
# so NO change to damage.py is required -- what was missing is a *rule* fixing
# sigma per analog. With sigma fixed at 3.0 (its default) the three analogs'
# YSZ phases differ only in volume fraction, not in length scale, and a YSZ
# damage operator would then have no morphological lever to act on.
#
# THE RULE (frozen 2026-08-10, before any Step 1 result was seen): sigma is
# proportional to the analog's own target particle diameter, anchored so the
# medium analog keeps the existing default of 3.0.
#
#     sigma(analog) = SIGMA_MEDIUM * D_particle(analog) / D_particle(medium)
#
# With the Design Memo sec 1.2 diameters (420 / 484 / 560 nm) this gives
# 2.603 / 3.000 / 3.471 voxels for fine / medium / coarse.
#
# HONEST NOTE ON WHAT THIS RULE CAN AND CANNOT DO, recorded before running:
# a sigma sweep at fixed Phi_YSZ (medium analog) measured YSZ component count
# falling monotonically with sigma -- 615 / 246 / 98 / 47 / 38 at sigma =
# 1.5 / 3.0 / 5.0 / 7.0 / 10.0. Real pristine YSZ runs the OTHER way: fragment
# density rises with coarseness (1.03 / 2.46 / 6.07 per Mvoxel for
# fine / medium / coarse). So this rule is expected to reproduce the YSZ
# LENGTH-SCALE ordering (gate G1-h) while working AGAINST the YSZ
# FRAGMENTATION ordering (gate G1-i, second clause). That is a prediction, not
# a hedge: it is recorded here so the gate outcome is a genuine test of the
# placement model rather than a story told afterwards.
# ---------------------------------------------------------------------------

SIGMA_MEDIUM_VOX = 3.0
D_PARTICLE_MEDIUM_NM = 484.0


def ysz_sigma_for_analog(d_particle_nm: float,
                         sigma_medium: float = SIGMA_MEDIUM_VOX,
                         d_medium_nm: float = D_PARTICLE_MEDIUM_NM) -> float:
    """A2.3 sigma for `add_ysz_pore`, proportional to particle diameter.

    `d_particle_nm` is the analog's OWN target particle diameter (2*R*voxel_nm),
    not the real anode's -- per the Step 1 instruction. Returns voxels.
    """
    if d_particle_nm <= 0:
        raise ValueError(f"d_particle_nm must be positive, got {d_particle_nm}")
    return float(sigma_medium) * float(d_particle_nm) / float(d_medium_nm)


# ===========================================================================
# A4 -- PARTICULATE YSZ GENERATOR WITH EXPLICIT CONTACTS AND SINTERING YIELD
#
# Added 2026-08-10 under PREREGISTRATION_V2_1.md (committed e62f30b) sec 3-4.
# Replaces, FOR PROJECT 2 ONLY, the thresholded-random-field placement in
# cmlib.damage.add_ysz_pore -- which is NOT modified or deleted, because
# Project 1 and Step 0/1 results depend on it.
#
# WHY. Step 1 measured that a smoothed random field gets the YSZ fragmentation
# trend BACKWARDS: at fixed Phi, raising sigma rescales the morphology without
# changing its topology, so component count FALLS with coarseness while the
# real anodes' RISES. A random field has no grains and no contacts, which is
# precisely the mechanism that makes a real coarse YSZ backbone fragile.
#
# ARCHITECTURE. YSZ grains sit on a lattice INTERPENETRATING the Ni lattice --
# offset by half a pitch on every axis, i.e. at the body centres of the Ni
# lattice cells. This is the arrangement that lets two particulate networks
# occupy one domain without either being an afterthought, and it reuses the
# Ni generator's own validated rasterizer rather than inventing a second one.
#
#   grains    explicit spheres, radius r_ysz, at recorded centres
#   contacts  explicit nearest-neighbour pairs on the YSZ lattice
#   sintered  contact present with probability p_sinter -> an explicit neck bar
#             is rasterized centre-to-centre, so connection is GUARANTEED
#   unsintered  no bar; the inter-grain gap (pitch - 2*r_ysz) leaves the pair
#             6-disconnected, so disconnection is also GUARANTEED
#
# The gap is what makes requirement A4.6 ("no reliance on accidental raster
# overlap") checkable: at p_sinter = 0 the number of YSZ components must equal
# the number of grains. That is exactly gate K0.
#
# Ni ALWAYS TAKES PRECEDENCE: the YSZ mask is clipped by the Ni mask, never the
# other way round, matching add_ysz_pore's convention that Ni is never touched.
# ===========================================================================

import numpy as np  # noqa: E402


LATTICE_BASIS = {
    "SC":  [(0.5, 0.5, 0.5)],
    "BCC": [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
    "FCC": [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
}
# touching-sphere packing fraction, and nearest-neighbour distance in units of
# the cube edge a
LATTICE_FPACK = {"SC": np.pi / 6, "BCC": np.pi * np.sqrt(3) / 8,
                 "FCC": np.pi / (3 * np.sqrt(2))}
LATTICE_NN = {"SC": 1.0, "BCC": np.sqrt(3) / 2, "FCC": 1 / np.sqrt(2)}


def max_phi_ysz(lattice: str, phi_ni: float) -> float:
    """Analytic ceiling on Phi_YSZ for DISJOINT grains coexisting with Ni.

    Disjoint grains cap out at the lattice touching-packing fraction, and
    lattice-placed grains lose ~Phi_Ni of their raw volume to Ni clipping:

        Phi_YSZ_max ~= f_pack * (1 - Phi_Ni)

    A resolvable (>= 1 voxel) inter-grain gap requires staying STRICTLY below
    this. Measured 2026-08-10: an SC YSZ lattice is INFEASIBLE for the fine
    analog -- it needs 119% of its own cap (0.421 required vs 0.355 available)
    -- which is why the SC architecture was rejected before K0 rather than
    tested and failed. BCC and FCC both clear all three analogs.
    """
    return LATTICE_FPACK[lattice] * (1.0 - phi_ni)


def ysz_lattice_geometry(shape, cube_a_vox, jitter_frac, rng,
                         lattice: str = "FCC"):
    """YSZ grain centres on an SC/BCC/FCC lattice interpenetrating the Ni phase.

    `cube_a_vox` is the CUBE EDGE, not the nearest-neighbour distance; those
    differ by `LATTICE_NN[lattice]`. Contacts are nearest-neighbour pairs,
    identified by DISTANCE rather than index arithmetic, so one code path
    serves all three bases.

    `jitter_frac` is a fraction of the NEAREST-NEIGHBOUR distance and must be
    small enough that no two neighbouring grains touch once rasterized -- worst
    case approach along a bond is 2 * jitter_frac * nn. K0 is the check.

    Returns (centres, pairs) keyed by integer grain id.
    """
    basis = LATTICE_BASIS[lattice]
    nn = LATTICE_NN[lattice] * cube_a_vox
    j = jitter_frac * nn
    # BOUNDARY MARGIN. Grain centres are allowed to fall OUTSIDE the domain by
    # up to one nearest-neighbour distance, and the rasterizer clips them.
    # Without this, no grain covers the extreme slices whenever the lattice is
    # incommensurate with the domain, and the phase CANNOT span by
    # construction -- measured 2026-08-10: at a_ysz = 46 in a 181-deep domain
    # the last grain centre sits at z = 161 and reaches z = 176 against a last
    # slice of 180, giving P_span = 0.0000 at every p_sinter including 1.0.
    # The Ni generator avoids this by pinning its z boundary layers to the
    # faces; a lattice at an unrelated pitch cannot, so it must overhang.
    ncell = [int(np.ceil(s / cube_a_vox)) + 2 for s in shape]
    centres, coords, k = {}, [], 0
    for iz in range(-1, ncell[0]):
        for iy in range(-1, ncell[1]):
            for ix in range(-1, ncell[2]):
                for (bz, by, bx) in basis:
                    c = [(iz + bz) * cube_a_vox, (iy + by) * cube_a_vox,
                         (ix + bx) * cube_a_vox]
                    c = [c[a] + rng.uniform(-j, j) for a in range(3)]
                    if all(-nn <= c[a] < shape[a] + nn for a in range(3)):
                        centres[k] = tuple(c)
                        coords.append(c)
                        k += 1
    if not centres:
        return {}, []
    pts = np.asarray(coords)
    keys = list(centres.keys())
    d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1)
    iu, ju = np.triu_indices(len(pts), k=1)
    sel = d2[iu, ju] <= (nn * 1.02) ** 2
    pairs = [(keys[a], keys[b]) for a, b in zip(iu[sel], ju[sel])]
    return centres, pairs


def draw_sintered(n_pairs: int, p_sinter: float, rng) -> np.ndarray:
    """Independent Bernoulli sintering of each candidate contact."""
    return rng.random(n_pairs) < p_sinter


def build_ysz_mask(centres, pairs, sintered, r_ysz, neck_w_vox, shape,
                   ni_mask):
    """Rasterize YSZ grains + sintered necks, clipped by Ni.

    Reuses `cmlib.synth.rasterize` unmodified: unsintered contacts are passed a
    neck width of 0, which that function skips. Ni always wins the overlap.
    """
    from .synth import rasterize
    widths = np.where(sintered, float(neck_w_vox), 0.0)
    ysz = rasterize(centres, pairs, r_ysz, widths, shape)
    return ysz & ~ni_mask


def solve_r_ysz_for_phi(centres, pairs, sintered, neck_w_vox, shape, ni_mask,
                        phi_target, r_lo=1.0, r_hi=None, max_iter=14):
    """PROTOCOL A: bisect grain core radius to hit a hard Phi_YSZ target.

    Phi_YSZ is measured on the CLIPPED mask against the FULL domain, which is
    how the real per-anode values are defined. Returns
    (r_final, ysz_mask, phi_achieved, n_iter).
    """
    dom = float(np.prod(shape))
    if r_hi is None:
        r_hi = 0.5 * float(np.min([np.ptp([c[a] for c in centres.values()])
                                   for a in range(3)])) or 20.0
    best = None
    for it in range(max_iter):
        r_mid = 0.5 * (r_lo + r_hi)
        m = build_ysz_mask(centres, pairs, sintered, r_mid, neck_w_vox, shape,
                           ni_mask)
        phi = m.sum() / dom
        if best is None or abs(phi - phi_target) < abs(best[2] - phi_target):
            best = (r_mid, m, phi, it + 1)
        if phi > phi_target:
            r_hi = r_mid
        else:
            r_lo = r_mid
    return best


def _add_capsule(out, c0, c1, radius, shape):
    """OR into `out` every voxel within `radius` of the segment c0->c1.

    WHY THIS EXISTS. `cmlib.synth.rasterize` draws a neck as an AXIS-ALIGNED
    bar, which is correct for the Ni generator (6-connected lattice, every
    contact axis-aligned) but silently wrong for BCC/FCC, whose nearest
    neighbours lie along body/face diagonals: the bar advances along one axis
    while staying at c0's coordinates in the others, so it never reaches the
    partner grain. Measured 2026-08-10: reusing the axis-aligned rasterizer for
    an FCC YSZ lattice gave P_span = 0.0000 even at p_sinter = 1.0. This is a
    true distance-to-segment capsule, valid for any contact direction.
    """
    c0 = np.asarray(c0, dtype=float)
    c1 = np.asarray(c1, dtype=float)
    lo = [max(0, int(np.floor(min(c0[a], c1[a]) - radius - 1))) for a in range(3)]
    hi = [min(shape[a], int(np.ceil(max(c0[a], c1[a]) + radius + 2)))
          for a in range(3)]
    if any(hi[a] <= lo[a] for a in range(3)):
        return
    zz, yy, xx = np.ogrid[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    d = c1 - c0
    L2 = float((d * d).sum())
    pz, py, px = zz - c0[0], yy - c0[1], xx - c0[2]
    if L2 <= 0:
        t = 0.0
    else:
        t = np.clip((pz * d[0] + py * d[1] + px * d[2]) / L2, 0.0, 1.0)
    dz, dy, dx = pz - t * d[0], py - t * d[1], px - t * d[2]
    out[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] |= \
        (dz * dz + dy * dy + dx * dx) <= radius * radius


def rasterize_ysz(centres, pairs, sintered, r_ysz, neck_w_vox, shape):
    """YSZ grains (spheres) OR sintered necks (capsules), any contact direction.

    Grains are spheres of radius `r_ysz`; a sintered contact adds a capsule of
    radius `neck_w_vox / 2` spanning centre to centre, so connection is
    guaranteed. An unsintered contact adds nothing, so the inter-grain gap
    (nearest-neighbour distance - 2*r_ysz) leaves the pair 6-disconnected.
    """
    out = np.zeros(shape, dtype=bool)
    for c in centres.values():
        _add_capsule(out, c, c, float(r_ysz), shape)
    rad = 0.5 * float(neck_w_vox)
    for (a, b), s in zip(pairs, np.asarray(sintered)):
        if s:
            _add_capsule(out, centres[a], centres[b], rad, shape)
    return out


def build_ysz_mask_v2(centres, pairs, sintered, r_ysz, neck_w_vox, shape,
                      ni_mask):
    """As `build_ysz_mask` but using the direction-agnostic rasterizer."""
    return rasterize_ysz(centres, pairs, sintered, r_ysz, neck_w_vox,
                         shape) & ~ni_mask


def solve_r_ysz_for_phi_v2(centres, pairs, sintered, neck_w_vox, shape,
                           ni_mask, phi_target, r_lo=1.0, r_hi=None,
                           max_iter=12):
    """PROTOCOL A radius bisection on the direction-agnostic rasterizer."""
    dom = float(np.prod(shape))
    best = None
    for it in range(max_iter):
        r_mid = 0.5 * (r_lo + r_hi)
        m = build_ysz_mask_v2(centres, pairs, sintered, r_mid, neck_w_vox,
                              shape, ni_mask)
        phi = m.sum() / dom
        if best is None or abs(phi - phi_target) < abs(best[2] - phi_target):
            best = (r_mid, m, phi, it + 1)
        if phi > phi_target:
            r_hi = r_mid
        else:
            r_lo = r_mid
    return best
