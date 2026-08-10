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
