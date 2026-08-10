"""
T5b generator scaffolding: mass-conservative, percentile-targeted neck
widening on a cubic sphere lattice.

SCOPE: this is NOT the full Family A/B/C generator (no config file, no
YSZ/pore placement, no damage model). It is scoped narrowly to the T5b
coupling-decision experiment, per the 2026-08-10 review decision recorded in
out/next/preregistration.md ("0. Refined primary causal question...") and
out/next/EXECUTION_SPEC.md. Family B proper is only built after T5b's
acceptance criteria (preregistration.md #0b) are checked against real T5b
output.

WHAT THE ORIGINAL T5 SCAFFOLD (scripts/next/t5_coupling_experiment.py) GOT
WRONG, AND THE FIX (full diagnosis in out/next/t5_coupling_decision_report.md)
--------------------------------------------------------------------------------
1. MASS CONSERVATION. Widening a neck by simply adding Ni voxels increases
   total Ni volume. In the original 300-neck-pair lattice this overshot the
   +/-5% Phi_Ni tolerance by 2x-25x even at the mildest tested point (11%-128%
   relative deviation across all 8 tested configurations).

   FIX: every sphere's radius is uniformly shrunk -- as a continuous float,
   then matched by DIRECT VOXEL-COUNT bisection (not an analytic
   surface-area approximation, which would be off by several voxels on a
   discrete grid) -- so that total-Ni-removed (smaller spheres) matches
   total-Ni-added (wider necks) to within a small residual. Every voxel
   change (added / removed / net residual) is computed exactly by rasterizing
   and diffing masks, and is returned in full, never silently absorbed.

2. PERCENTILE TARGETING. "Widen the bottom K necks (a fixed count) to a fixed
   target width" cannot move the MEASURED p10 once the target width exceeds
   the untouched population's range: the widened group simply ranks ABOVE the
   untouched necks and stops occupying low-percentile ranks, so p10 is set by
   whichever untouched values are left at the bottom -- independent of the
   target width chosen. Confirmed both analytically (rank arithmetic: with a
   fixed K=60 out of n=300, the new p10 equals the ORIGINAL (10%+20%)=30th
   percentile of the untouched population, for ANY target width bigger than
   that) and empirically (T5's measured p10 was stuck at exactly 80 nm across
   four different target widths).

   FIX: a MAX-CLIP rule -- every neck narrower than a threshold T is raised to
   EXACTLY T (`new_width = max(old_width, T)`). This keeps the raised group at
   the BOTTOM of the sorted distribution (never pushed above the untouched
   population, because untouched necks that were already >= T are left
   completely alone), so once the fraction of necks originally below T
   exceeds ~10%, the new p10 equals T exactly, for T chosen by BISECTION
   against the ACTUAL MEASURED (SNOW) p10 -- not just the intended
   construction value, because watershed/SNOW segmentation can shift the
   realised value away from the raw geometric intention (as already found for
   watershed particle size under R1).

3. Uniform (whole-distribution) widening is kept as a secondary/exploratory
   POSITIVE CONTROL, per the review decision -- not the primary axis, but
   useful because its measured-p10 response is already well characterised
   (T5: delta = target_w - median(base_widths) predicts the resulting p10
   almost exactly) and does not need iterative bisection.

GEOMETRY: cubic lattice of Ni spheres, nearest-neighbour pairs joined by a
SHORT connecting neck spanning CENTRE TO CENTRE (never a long freestanding
bridge -- see the T5b design-rule finding in
scripts/next/phase0_validate_synthetic_pipeline.py T5b, which showed a long
freestanding bridge can spawn its own spurious watershed region). Because the
neck bar spans exactly from one sphere centre to the other, connectivity
between the two spheres is guaranteed regardless of how much the sphere
radius is later shrunk for mass conservation (the bar always includes both
centres).
"""

from __future__ import annotations

import numpy as np

# ===========================================================================
# FAMILY B DISORDERED PILOT (added 2026-08-10, post-T5b review)
# ===========================================================================
#
# The T5b lattice (perfect cubic grid, plain discrete-uniform base neck
# widths) was accepted as a valid DIAGNOSTIC scaffold, not the main-claim
# structure (preregistration.md #0c amendment D). Two things change here:
#
# 1. GEOMETRY: `jittered_lattice_geometry` breaks positional regularity by
#    randomly displacing sphere centres (a "perturbed lattice", the simpler
#    of the two options the review offered over a full random packing).
#    Topology (which spheres are candidate neighbours) stays the deterministic
#    lattice adjacency -- this is a deliberate scoping choice for a first
#    pilot: it avoids the disconnection/instability risk of a
#    distance-threshold-based random-packing topology, at the cost of capping
#    achievable coordination at 6 (interior) rather than reaching the
#    "approximately 6-8" the review asked for "if possible". Reported
#    honestly via the pilot's own `mean_degree` diagnostic, not silently
#    overstated. The z-axis (percolation axis) boundary layers keep their
#    EXACT on-lattice z-coordinate so the domain-spanning guarantee
#    established in T5/T5b still holds; all other coordinates (interior z,
#    and every sphere's y,x) are independently jittered.
#
# 2. BASE NECK-WIDTH DISTRIBUTION: a plain discrete-uniform draw cannot reach
#    p50/p10 >= 2.5 while also keeping the floor at the required 3-4 voxel
#    resolution within a physically sensible neck-width range (checked
#    numerically: for discrete-uniform[lo,hi], p50/p10 -> 1 as lo grows
#    relative to the range, so lo>=4 forces an unreasonably large hi to reach
#    ratio 3). `mixture_neck_widths` instead draws a MIXTURE: a "normal"
#    neck population plus a genuinely narrow-tail subpopulation -- which also
#    happens to be more physically realistic than a uniform draw (real
#    necks are not uniformly distributed; there is typically a bulk
#    population plus some genuinely weak/narrow necks from imperfect
#    sintering). `draw_valid_base_widths` enforces the pre-registered
#    validity criterion (p50/p10 >= 2.5) BEFORE any widening is attempted,
#    by rejection sampling -- every rejected attempt is logged, never
#    silently discarded (this is exactly the mechanism that would have
#    caught T5b's seed=1 anomaly before it ever reached the widening step).
# ===========================================================================


def jittered_lattice_geometry(nlat: int, pitch_vox: int, r_vox: float,
                              margin: int, jitter_frac: float,
                              rng: np.random.Generator):
    """Perturbed cubic lattice: same topology as `lattice_geometry`, but every
    sphere centre is independently displaced by up to
    +/- jitter_frac * pitch_vox per axis, EXCEPT the z-coordinate of the
    iz=0 and iz=nlat-1 layers, which stays exactly on-lattice so the
    domain-spanning guarantee (see `lattice_geometry` docstring) still holds.

    Returns (centres, pairs, shape) -- same shapes/types as `lattice_geometry`.
    """
    j = jitter_frac * pitch_vox
    centres = {}
    for iz in range(nlat):
        for iy in range(nlat):
            for ix in range(nlat):
                z0 = iz * pitch_vox
                y0 = margin + r_vox + iy * pitch_vox
                x0 = margin + r_vox + ix * pitch_vox
                if 0 < iz < nlat - 1:
                    z0 = z0 + rng.uniform(-j, j)
                y0 = y0 + rng.uniform(-j, j)
                x0 = x0 + rng.uniform(-j, j)
                centres[(iz, iy, ix)] = (z0, y0, x0)

    pairs = []
    for (iz, iy, ix) in centres:
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            nb = (iz + d[0], iy + d[1], ix + d[2])
            if nb in centres:
                pairs.append(((iz, iy, ix), nb))

    nz = (nlat - 1) * pitch_vox + 1
    n_xy = nlat * pitch_vox + 2 * int(r_vox) + 2 * margin
    return centres, pairs, (nz, n_xy, n_xy)


def platform_v2_lattice_geometry(nlat_z: int, nlat_xy: int, pitch_vox: int,
                                 r_vox: float, margin: int, jitter_frac: float,
                                 rng: np.random.Generator):
    """Platform-v2 Ni generator qualification (2026-08-10): same construction
    as `jittered_lattice_geometry`, generalised to a SEPARATE particle count
    along z (`nlat_z`) vs y/x (`nlat_xy`).

    WHY SEPARATE COUNTS. The domain-spanning guarantee requires the z-boundary
    layers (iz=0, iz=nlat_z-1) to sit exactly on the domain faces, so
    nz = (nlat_z-1)*pitch + 1 (no margin term). The y/x extent instead carries
    a margin on both sides: n_xy = nlat_xy*pitch + 2*r + 2*margin. These two
    formulas cannot be satisfied by the SAME particle count while both landing
    in the target domain range (160-192 vox on every axis) at a pitch set by
    the Phi_Ni target -- solved here by decoupling nlat_z from nlat_xy rather
    than distorting pitch or margin to force symmetry. This function does NOT
    change the topology rule (still plain nearest-neighbour/6-connectivity
    adjacency, i.e. `pairs` connects only axis-adjacent lattice sites) -- see
    scripts/platform_v2/design_probe.py for the empirical check of whether
    this plain topology already lands in the coordination target band before
    any topology modification (face-diagonal bonds, etc.) is considered.

    Returns (centres, pairs, shape) -- same shapes/types as
    `jittered_lattice_geometry`.
    """
    j = jitter_frac * pitch_vox
    centres = {}
    for iz in range(nlat_z):
        for iy in range(nlat_xy):
            for ix in range(nlat_xy):
                z0 = iz * pitch_vox
                y0 = margin + r_vox + iy * pitch_vox
                x0 = margin + r_vox + ix * pitch_vox
                if 0 < iz < nlat_z - 1:
                    z0 = z0 + rng.uniform(-j, j)
                y0 = y0 + rng.uniform(-j, j)
                x0 = x0 + rng.uniform(-j, j)
                centres[(iz, iy, ix)] = (z0, y0, x0)

    pairs = []
    for (iz, iy, ix) in centres:
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            nb = (iz + d[0], iy + d[1], ix + d[2])
            if nb in centres:
                pairs.append(((iz, iy, ix), nb))

    nz = (nlat_z - 1) * pitch_vox + 1
    n_xy = nlat_xy * pitch_vox + 2 * int(r_vox) + 2 * margin
    return centres, pairs, (nz, n_xy, n_xy)


def mixture_neck_widths(n_pairs: int, rng: np.random.Generator,
                        frac_weak: float, weak_range: tuple,
                        normal_range: tuple) -> np.ndarray:
    """Base neck widths drawn from a mixture: `frac_weak` of necks from a
    genuinely narrow `weak_range`, the rest from a `normal_range`. Both
    ranges are (lo, hi) inclusive, voxels. See module docstring for why a
    plain uniform draw cannot reach the required p50/p10 spread at the
    required resolution floor simultaneously."""
    is_weak = rng.random(n_pairs) < frac_weak
    weak = rng.integers(weak_range[0], weak_range[1] + 1, n_pairs)
    normal = rng.integers(normal_range[0], normal_range[1] + 1, n_pairs)
    return np.where(is_weak, weak, normal).astype(float)


def base_distribution_stats(widths: np.ndarray) -> dict:
    p10 = float(np.percentile(widths, 10))
    p50 = float(np.percentile(widths, 50))
    return {"p10_vox": p10, "p50_vox": p50,
           "ratio": p50 / p10 if p10 > 0 else float("inf"),
           "min_vox": float(widths.min()), "max_vox": float(widths.max())}


def draw_valid_base_widths(n_pairs: int, seed: int, frac_weak: float,
                           weak_range: tuple, normal_range: tuple,
                           min_ratio: float = 2.5, max_attempts: int = 20):
    """Rejection-sample a base neck-width draw until p50/p10 >= min_ratio.

    Tries `seed`, then `seed*100000 + attempt` for attempt=1,2,...  Every
    attempt (accepted or rejected) is logged with its stats, per
    preregistration.md #0c amendment E ("record all rejected base seeds and
    reasons"). Returns (widths, accepted_attempt_seed, n_attempts, log) where
    `log` is a list of dicts, one per attempt, each with the RNG seed used,
    the resulting p10/p50/ratio, and whether it was accepted.
    """
    log = []
    for attempt in range(max_attempts):
        sub_seed = seed if attempt == 0 else seed * 100_000 + attempt
        rng = np.random.default_rng(sub_seed)
        widths = mixture_neck_widths(n_pairs, rng, frac_weak, weak_range,
                                     normal_range)
        stats = base_distribution_stats(widths)
        accepted = stats["ratio"] >= min_ratio
        log.append({"attempt": attempt, "sub_seed": sub_seed,
                    "accepted": accepted, **stats})
        if accepted:
            return widths, sub_seed, attempt + 1, log
    # exhausted attempts: return the best (highest-ratio) draw found, but the
    # caller MUST check the log / n_attempts and treat this seed as a
    # documented failure, never silently accept it as if it had passed.
    best = max(log, key=lambda r: r["ratio"])
    rng = np.random.default_rng(best["sub_seed"])
    widths = mixture_neck_widths(n_pairs, rng, frac_weak, weak_range,
                                 normal_range)
    return widths, best["sub_seed"], max_attempts, log


def lattice_geometry(nlat: int, pitch_vox: int, r_vox: float, margin: int):
    """Cubic lattice of `nlat`^3 sphere centres.

    The z-axis (axis 0) is sized so the first and last sphere LAYERS touch
    the domain faces exactly -- required for `cmlib.pnm`'s spanning-cluster
    restriction to accept the lattice as a single network (found empirically
    while building the original T5 scaffold: a margined z-axis silently
    yields an EMPTY network because nothing touches the domain boundary).
    y and x get `margin` padding (irrelevant to spanning, avoids boundary
    truncation of border spheres in the particle-size statistics).

    Returns (centres: {(iz,iy,ix): (z,y,x)}, pairs: [((iz,iy,ix),(iz,iy,ix)), ...],
    shape: (nz, ny, nx)).
    """
    centres = {}
    for iz in range(nlat):
        for iy in range(nlat):
            for ix in range(nlat):
                z = iz * pitch_vox
                y = margin + r_vox + iy * pitch_vox
                x = margin + r_vox + ix * pitch_vox
                centres[(iz, iy, ix)] = (z, y, x)

    pairs = []
    for (iz, iy, ix) in centres:
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            nb = (iz + d[0], iy + d[1], ix + d[2])
            if nb in centres:
                pairs.append(((iz, iy, ix), nb))

    nz = (nlat - 1) * pitch_vox + 1
    n_xy = nlat * pitch_vox + 2 * int(r_vox) + 2 * margin
    return centres, pairs, (nz, n_xy, n_xy)


def base_neck_widths(n_pairs: int, rng: np.random.Generator, lo: int = 2,
                     hi: int = 6) -> np.ndarray:
    """Randomized base (pristine) neck widths, voxels, discrete uniform [lo,hi]."""
    return rng.integers(lo, hi + 1, size=n_pairs).astype(float)


def max_clip_widths(base_widths: np.ndarray, threshold: float) -> np.ndarray:
    """Percentile-targeted lower-tail widening: raise everything below
    `threshold` up TO `threshold`; leave everything already >= threshold
    untouched. This is what keeps the raised group at the BOTTOM of the
    distribution -- see module docstring point 2."""
    return np.maximum(base_widths, threshold)


def uniform_shift_widths(base_widths: np.ndarray, delta: float,
                         lo: float, hi: float) -> np.ndarray:
    """Secondary/exploratory positive control: shift every neck by the same
    increment, clipped to [lo, hi] (hi should be < sphere diameter so necks
    never fully engulf a sphere)."""
    return np.clip(base_widths + delta, lo, hi)


def rasterize(centres: dict, pairs: list, r_vox: float,
             neck_widths_vox: np.ndarray, shape) -> np.ndarray:
    """Spheres of radius r_vox at every centre, OR'd with a short connecting
    bar (spanning centre to centre) of the given width for every pair.
    A width of 0 means "no neck" for that pair."""
    nz, ny, nx = shape
    zz, yy, xx = np.ogrid[:nz, :ny, :nx]
    ni = np.zeros(shape, dtype=bool)
    for (z, y, x) in centres.values():
        ni |= (zz - z) ** 2 + (yy - y) ** 2 + (xx - x) ** 2 < r_vox ** 2

    for (p0, p1), w in zip(pairs, neck_widths_vox):
        if w <= 0:
            continue
        c0 = np.array(centres[p0])
        c1 = np.array(centres[p1])
        axis = int(np.argmax(np.abs(c1 - c0)))
        lo = int(min(c0[axis], c1[axis]))
        hi = int(max(c0[axis], c1[axis]))
        half = int(round(w)) // 2
        wv = max(int(round(w)), 1)
        sl = [slice(None)] * 3
        sl[axis] = slice(lo, hi + 1)
        for a in range(3):
            if a != axis:
                sl[a] = slice(int(c0[a]) - half, int(c0[a]) - half + wv)
        ni[tuple(sl)] = True
    return ni


def match_radius_for_mass_conservation(centres: dict, pairs: list,
                                       neck_widths_vox: np.ndarray,
                                       r_base: float, shape,
                                       target_total_voxels: int,
                                       r_lo: float = 1.0,
                                       max_iter: int = 12):
    """Bisect the sphere radius (a float, EXACT voxel-count matched, not an
    analytic surface-area approximation) so that
    rasterize(centres, pairs, r, neck_widths_vox, shape).sum() is as close as
    possible to `target_total_voxels`.

    Returns (r_final, ni_mask_final, actual_total_voxels).
    """
    r_hi = r_base
    lo_mask = rasterize(centres, pairs, r_lo, neck_widths_vox, shape)
    if lo_mask.sum() > target_total_voxels:
        # even the smallest allowed radius overshoots -- cannot conserve mass
        # by radius shrinkage alone at this neck width; return the smallest
        # tried, caller must inspect and report the residual honestly.
        return r_lo, lo_mask, int(lo_mask.sum())

    best_r, best_mask, best_n = r_base, None, None
    for _ in range(max_iter):
        r_mid = 0.5 * (r_lo + r_hi)
        mask = rasterize(centres, pairs, r_mid, neck_widths_vox, shape)
        n = int(mask.sum())
        if best_n is None or abs(n - target_total_voxels) < abs(best_n - target_total_voxels):
            best_r, best_mask, best_n = r_mid, mask, n
        if n > target_total_voxels:
            r_hi = r_mid          # too big -> shrink further
        else:
            r_lo = r_mid          # too small -> grow back
    return best_r, best_mask, best_n


def build_mass_conservative_structure(centres, pairs, shape, r_base_vox,
                                      base_widths_vox, final_widths_vox):
    """Build a widened, mass-conservative structure and report the exact
    voxel-level accounting (no analytic approximation, no silent absorption).

    Returns a dict:
        ni_mask, r_final_vox,
        voxels_base, voxels_added_by_necks, voxels_removed_by_shrink,
        voxels_final, voxels_net_residual, phi_Ni_base, phi_Ni_final
    """
    domain_vox = int(np.prod(shape))
    base_mask = rasterize(centres, pairs, r_base_vox, base_widths_vox, shape)
    n_base = int(base_mask.sum())

    widened_no_compensation = rasterize(centres, pairs, r_base_vox,
                                        final_widths_vox, shape)
    n_widened_uncompensated = int(widened_no_compensation.sum())
    voxels_added_by_necks = n_widened_uncompensated - n_base

    r_final, final_mask, n_final = match_radius_for_mass_conservation(
        centres, pairs, final_widths_vox, r_base_vox, shape,
        target_total_voxels=n_base)
    voxels_removed_by_shrink = n_widened_uncompensated - n_final

    return {
        "ni_mask": final_mask,
        "r_base_vox": r_base_vox,
        "r_final_vox": r_final,
        "voxels_base": n_base,
        "voxels_added_by_necks": voxels_added_by_necks,
        "voxels_removed_by_shrink": voxels_removed_by_shrink,
        "voxels_final": n_final,
        "voxels_net_residual": n_final - n_base,
        "phi_Ni_base": n_base / domain_vox,
        "phi_Ni_final": n_final / domain_vox,
    }
