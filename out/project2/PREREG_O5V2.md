# Pre-registration — O5v2: true volume-conserving curvature-driven agglomeration

**Frozen 2026-08-11, committed BEFORE any O5v2 code is written or run.**
`cmlib/damage.py`, `cmlib/synth.py` frozen; O5v2 goes in `cmlib/damage2.py`.
**No erosion parameter is used — this is not erosion.**

## Hypothesis
Volume-conserving agglomeration consumes critical Ni necks, killing percolation
in fine first (C1: real retention 0.680 < 0.855 < 0.947) while preserving
Ni–YSZ interfacial area, retaining TPB (C3: 0.799 > 0.621 > 0.574).

## Operator (Option A — morphological curvature flow), frozen
Per round: surface = `ni & ~erode(ni, STRUCT6)`; discrete mean curvature
`k = (#pore 6-neighbours) − (#Ni 6-neighbours)` (positive convex, negative
concave); remove the `N_move` most-convex surface voxels; add at exactly
`N_move` most-concave pore sites adjacent to Ni. Volume conserved by
construction. YSZ never touched.

**Rate, frozen:** `N_move = round(p_coarsen × N_surface_pristine)`,
**`p_coarsen = 0.03`**. Not adjusted after the gate.

If (ii) fails, the curvature proxy is insufficient: upgrade to 26-connectivity
curvature or Option B (KMC). **That is an implementation fix; `p_coarsen` and
`N_move` are NOT touched.**

## Validity gate A1v2 — all must hold, every ROI and seed
1 ROI per anode, n = 1, 3, 5, damage seeds 300/301/302.
- (i) |ΔΦ_Ni| ≤ 0.005 absolute for all n
- (ii) `S_spec(1) < S_spec(0)`, strict, monotonic
- (iii) `TPB(n) ≤ TPB(0)` for all n
- (iv) `R_Ni(n) ≤ R_Ni(0)` and non-increasing
- (v) YSZ untouched (array equality)

**Any failure ⇒ STOP, report, adjust nothing, do not bisect.**

## Bisection (only if the gate passes)
Same real ROIs as O6 (fine/medium 8 µm ×3, coarse 12 µm ×3), 3 damage seeds,
integer bisection [1,20] expand-only width ≤ 1, thresholds **0.50 and 0.10** on
`R_Ni = largest spanning Ni cluster / pristine Ni voxels`.

## Decision rules, frozen
**PRIMARY C1-real:** fine reaches threshold at strictly lower intensity than
medium AND coarse by ≥ 1.0 round, both thresholds, ≥ 3 ROIs/anode. One
threshold only = partial.
**MECHANISM:** raw and partial Spearman (controlling pristine `P_span`) vs
specific Ni surface area and vs corrected min-cut fraction; |ρ| ≥ 0.6 with
leave-one-out sign consistency; if both correlate report both.
**SECONDARY C3-real:** TPB retention reported only, no gate.
**If C1 fails with a validated O5v2:** definitive scope boundary for
volume-conserving curvature-driven coarsening. Combined with the O6 boundary for
surface-erosion operators, **both principal geometry-based degradation classes
are then tested and neither reproduces the ordering.** Report and STOP building
operators in these two classes.

Ordinal comparisons only. No fitting.
