# E0 vertical slice — report

> **Correction (2026-08-10, post-review):** this report originally stated
> "8 intensities × 45 = 360 damage runs" in three places. The actual sweep
> used 6 intensities (`{2,5,6,7,8,10}`), giving 270 runs — verified directly
> against `e0_combined_all_intensities.csv`
> (`sorted(df.n_rounds.unique())` = `[2, 5, 6, 7, 8, 10]`, `len(df) = 270`).
> Fixed in place below; flagged by external review, not self-caught. No
> other number in this report was affected — the per-intensity retained-value
> tables, the saturation finding, and the "no differentiation by p10 ratio"
> finding are all independent of this count and unchanged.

Run 2026-08-10. **Exploratory, non-confirmatory pipeline spike**, per the
review decision. This is NOT the main Family B experiment, NOT a valid test
of the scientific hypothesis, NOT a basis for Path B, for tuning the damage
model post-hoc, for Family C, or for claiming a design principle. Code:
`scripts/spike/e0_vertical_slice.py` (pre-registered run, n_rounds∈{2,5,10}),
`scripts/spike/e0b_saturation_bridge.py` (supplementary run, n_rounds∈{6,7,8},
added per the pre-registered saturation rule below — no other parameter
changed). Data: `e0_vertical_slice.csv`, `e0b_saturation_bridge.csv`,
`e0_combined_all_intensities.csv`. Figures:
`e0_retained_pspan_vs_p10.png`, `e0_combined_retained_pspan_vs_intensity.png`.

## What this spike set out to test, and the answer to each

**1. Can minimal YSZ/pore placement be added without corrupting the
validated Ni structures? YES.** Ni voxels are never touched by construction
(YSZ/pore placement only assigns labels to the non-Ni remainder). Verified
explicitly, not just assumed: Ni `P_span` was checked before and after
placement for all 15 structures and found bit-identical in every case.
Achieved composition (e.g. seed 0, base): Φ_Ni=0.330, Φ_YSZ=0.347,
Φ_pore=0.324. **Caveat, exactly as anticipated:** Φ_Ni could not be forced to
the nominal medium-anode value (0.250) because the Ni geometry is fixed and
unchanged from the already-validated Family B pilot (Φ_Ni here is whatever
the pilot's lattice produces, ~0.32–0.33). What *was* targeted and achieved
is the YSZ:pore *ratio* of the non-Ni remainder, matching the medium anode's
own Φ_YSZ/(Φ_YSZ+Φ_pore) = 0.5173.

**2. Can a non-circular D4 model be implemented and produce non-saturating
damage? Partially — see the saturation finding below, which is itself an
answer, not a failure.** D4 (fixed 1-voxel oxidative expansion into pore only,
then stochastic surface erosion at p_erode=0.35 per round, then keep only the
largest remaining Ni component) never inspects or selects on the measured
neck-p10 variable — it is a uniform geometric rule applied identically
regardless of a structure's neck statistics. TPB density itself is
**not** saturating — it degrades smoothly and monotonically across every
intensity tested (mean retained TPB: 4.21 at n=2 → 0.90 at n=5 → 0.34 at n=6
→ 0.11 at n=7 → 0.028 at n=8 → 0.0007 at n=10 — note n=2's retained TPB
>1 is real, not an error: the fixed 1-voxel expansion step adds surface area
faster than 2 rounds of mild erosion remove it, a real, if slightly
counterintuitive, emergent consequence of running expansion and erosion in
sequence at low intensity). What *does* saturate sharply is `P_span`: exactly
1.0 for every one of 45 (5 seeds × 3 damage seeds) structures at n_rounds
2, 5, 6, 7, and 8, and exactly 0.0 for every one of 45 at n_rounds 10 — a
step transition with no intermediate value observed anywhere between n=8 and
n=10 across the whole pre-registered + bridge sweep (6 intensities [2,5,6,7,8,10]
× 45 combinations = 270 damage runs total — corrected 2026-08-10 from an
arithmetic error that had said "8 intensities × 45 = 360"; verified directly
against `e0_combined_all_intensities.csv`: `sorted(df.n_rounds.unique())` =
`[2, 5, 6, 7, 8, 10]`, `len(df)` = `270`). Per the pre-registered rule (D): this
is reported as saturation and the *intensity range* was adjusted (the 6/7/8
bridge), not the damage model or any scientific criterion — p_erode and
expand_vox are unchanged from the original pre-registration throughout.

**3. Can retained percolation and retained TPB be recomputed consistently
after damage? YES.** Every metric (Φ_Ni, P_span, P_reach, TPB density,
largest-component fraction) was computed identically pre- and post-damage
using the same `cmlib.api` functions used throughout this project, for all
270 damage runs, with no failures, no NaN propagation issues, and physically
coherent results (e.g. at the P_span=0 transition, `largest_component_fraction`
stays at 1.0 for every one of the 45 collapsed structures — meaning the Ni
network necks down to a single non-spanning remnant rather than fragmenting
into many disconnected islands, a coherent, interpretable failure mode for
this damage mechanism, not a computational artifact).

**4. Is there a directional hint that achieved neck-p10 ratio affects
retained percolation? NO — and this is reported exactly as the pre-registered
rules require, not softened.** See below.

**5. Pipeline failures found cheaply? None.** All 15 structures reconstructed
bit-identical to the recorded pilot CSV (verified programmatically, not just
assumed); YSZ/pore placement, D4 damage, and metric recomputation ran
without error across 270 damage evaluations; total wall time for both runs
combined was under 6 minutes.

## The directional-signal result, read against the pre-registered rules (§F)

Retained P_span by nominal group, at every intensity tested:

| n_rounds | nominal 1.0× (achieved ~1.00×) | nominal 1.5× (achieved ~1.39×) | nominal 2.0× (achieved ~1.99×) |
|---|---|---|---|
| 2 | 1.000 | 1.000 | 1.000 |
| 5 | 1.000 | 1.000 | 1.000 |
| 6 | 1.000 | 1.000 | 1.000 |
| 7 | 1.000 | 1.000 | 1.000 |
| 8 | 1.000 | 1.000 | 1.000 |
| 10 | 0.000 | 0.000 | 0.000 |

**Every column is identical at every row.** Retained TPB tells the same
story (values in the combined-run table agree to 3 significant figures
across all three nominal groups at every intensity — e.g. at n_rounds=6:
0.342, 0.345, 0.345). The combined figure
(`e0_combined_retained_pspan_vs_intensity.png`) shows the three group curves
exactly superimposed — only one color is visible because the other two are
drawn directly underneath it.

**Per the pre-registered rule (F): "If there is no effect, call this
inconclusive on the hypothesis and use it to debug the pipeline."** That is
the finding. There is no effect, at any intensity tested, at any resolution
this spike could distinguish. This is called inconclusive, not negative —
the pre-registered rules explicitly distinguish "no effect" from "Path B,"
and this spike is barred from declaring either Path A or Path B regardless.

## Why "no effect visible here" is expected, and what it does and does not imply

This is a 64-particle system with the geometry's own connectivity graph
"neck-limited" by everything, not just the widened necks specifically — a
domain-spanning path in a 4×4×4 lattice with mean coordination ~4.3 has very
few alternative routes, so ANY sufficiently deep erosion is more likely to
sever the *whole* structure via its single most exposed weak point (wherever
that happens to be) than to reveal a graded, neck-distribution-dependent
retention curve. This is exactly the kind of resolution limit that
platform qualification (higher coordination, larger domain, more seeds) is
supposed to fix — it is a reason FOR that qualification step, not a
finding that supersedes it. **This spike does not, and is not permitted to,
argue that coordination/scale-up qualification can be skipped.**

## Explicit compliance with the "this spike is not" list

- Not the main Family B experiment — confirmed, no claim is made about the
  final hypothesis.
- Not a valid test of the final scientific hypothesis — confirmed, "no
  effect visible in this spike" is reported as inconclusive-and-pipeline-
  informative only.
- Not a reason to skip platform v2 — confirmed above.
- D4's parameters (p_erode=0.35, expand_vox=1) were never adjusted after
  seeing outcomes; only the intensity GRID (n_rounds) was extended, exactly
  as the pre-registered saturation rule (D) permits.
- Not a basis for Family C — nothing beyond this spike has been built.
- Not a basis for claiming a design principle — none is claimed.

## Stop condition

Per instruction H: **stopping here.** Not built: platform v2 (higher
coordination, larger domain), any larger structure, Family C, or any
calibration against the real Holzer/Pecho dataset. Awaiting review.
