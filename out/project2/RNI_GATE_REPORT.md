# R_Ni metric gate — **PASS.** O6 is now usable; bisection not yet run.

Under `PREREG_RNI_METRIC.md` (frozen `bfc5d89`). Data: `c1real_rni_gate.csv`.
`p_erode` unchanged at 0.35. `cmlib/damage.py`, `cmlib/synth.py` frozen.

## Mandatory sanity check — exact

| anode | pristine `P_span` | `R_Ni(0)` | \|diff\| |
|---|---|---|---|
| fine | 0.982135 | 0.982135 | **0.00e+00** |
| medium | 0.971254 | 0.971254 | **0.00e+00** |
| coarse | 0.887775 | 0.887775 | **0.00e+00** |

## R_Ni under O6 — monotone in all three anodes

| anode | R_Ni(0) | n=1 | n=3 | n=5 | n=8 | vol loss n=8 |
|---|---|---|---|---|---|---|
| fine | 0.9821 | 0.9505 | 0.8692 | 0.7540 | **0.6086** | 0.391 |
| medium | 0.9713 | 0.9467 | 0.8867 | 0.8126 | **0.6988** | 0.301 |
| coarse | 0.8878 | 0.8632 | 0.7965 | 0.7212 | **0.5765** | 0.424 |

**The pruning artifact is gone.** R_Ni declines monotonically from a pristine
value that is not forced to 1.0. **The metric fix works and the operator is now
scoreable.**

## Amendment F — answered from source, not asserted

`cmlib/damage2.py:309` — removal sites are
`boundary & (rng.random(cur.shape) < p_erode)`: **uniform over all surface
voxels, no curvature weighting.** The rate argument in this milestone is
therefore about **surface area**, not curvature.

## TPB — reported only, artifact confirmed

TPB rises 7–15× at n=1 then collapses to ~0.1–0.5 µm⁻² by n=8. **O6 is
volume-losing and surface-area-increasing; real coarsening is volume-conserving
and surface-area-reducing.** Quoting the ruling: *"C3-real (TPB retention) is
not addressable with O6. The TPB retention puzzle requires a volume-conserving,
surface-area-reducing coarsening operator, which is a separate milestone."*
No gate, no conclusion.

## An early signal, flagged but NOT a result

At n = 8 the R_Ni ordering is **coarse 0.577 < fine 0.609 < medium 0.699**,
against real retention **fine 0.680 < medium 0.855 < coarse 0.947**. Fine is not
lowest, and coarse — best in reality — is worst here. **This is one ROI per
anode at one intensity and is not the pre-registered test**; the bisection at
both frozen thresholds decides. It is recorded now so it cannot look like a
post-hoc observation later.

Note coarse starts far lower (0.8878 vs 0.9821), which is exactly why
Amendment C's partial correlations controlling for pristine `P_span` are
required.

## Status — what remains

**Not run:** Amendment D (v3 audit re-run with corrected extents, mandatory
**before** the bisection) and the C1-real bisection itself with its Spearman
analysis. Both are specified and unblocked; neither was reached in this pass.
No result is claimed for C1-real.
