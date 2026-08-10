# Pre-registration amendment — R_Ni outcome metric (frozen before any run)

Amends `PREREG_O6.md` (`f980562`). Frozen 2026-08-11.

## Metric

**R_Ni(n) = spanning-cluster Ni voxels at intensity n / PRISTINE Ni voxels.**

Monotone non-increasing under erosion; unaffected by island pruning (deleting
non-spanning clusters changes neither the spanning cluster nor the pristine
denominator); directly comparable to the real retention values
(0.680 / 0.855 / 0.947).

**Sanity check, mandatory before the bisection:** `R_Ni(0)` must equal pristine
`P_span`. If not, the metric is implemented wrongly — STOP and fix.

**Transition** = first intensity at which R_Ni falls below a threshold. **Two
frozen thresholds, 0.50 and 0.10; both reported.**

## Amendment F — erosion site selection, answered from the source

`cmlib/damage2.py::apply_o6` selects removal sites as
`boundary & (rng.random(cur.shape) < p_erode)` — **uniform over all surface
voxels, with no curvature weighting.** Therefore **the rate argument in this
milestone is about SURFACE AREA** (more surface ⇒ more removal sites), not about
curvature. Recorded so the mechanism interpretation cannot drift.

## TPB — reported only

O6 is volume-losing and surface-area-increasing; real Ni coarsening is
volume-conserving and surface-area-reducing. **C3-real is not addressable with
O6.** TPB values are diagnostic only, with the roughening artifact noted. No TPB
gate, no TPB conclusion.

## Amendment D — mandatory ordering

The v3 graph audit must be re-run with corrected axis extents **before** the
bisection, so the min-cut fractions entering the Spearman analysis are final.

## Decision rules (frozen, as instructed)

PRIMARY C1-real: fine reaches the threshold at strictly lower intensity than
medium AND coarse, ≥ 1.0 round, ≥ 3 ROIs per anode, **at both thresholds**; one
threshold only = partial. MECHANISM: raw and partial Spearman (controlling for
pristine `P_span`) vs specific Ni surface area and vs corrected min-cut
fraction; |ρ| ≥ 0.6 with leave-one-out sign consistency; if both correlate,
report both and let the partials decide dominance. C1-real failure with a
validated operator and well-posed metric = definitive scope boundary; stop
building operators of this class.

## Out of scope for this milestone

TPB retention puzzle; YSZ damage mechanism; volume-conserving coarsening
operator. Separate milestones; not to be conflated.
