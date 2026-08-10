# Pre-registration — neck-size decoupling study

Frozen 2026-08-10, before Phase 3 (or any generator/damage-sweep) is run.
Changes after this point must be logged with a reason, not silently made.
This document is the single source of truth for what counts as a positive
result; `out/next/EXECUTION_SPEC.md` describes *how* the work is done, this
document describes *what would make it count*.

> **AMENDMENT 2026-08-10 (post-T5).** The T5 coupling-decision experiment
> (`out/next/t5_coupling_decision_report.md`) found the scaffold's naive
> construction failed for two specific, fixable reasons (added Ni mass not
> compensated; a percentile-definition trap in "widen bottom 20%"), not
> because decoupling is infeasible. Per the review decision, the primary
> causal question and the particle-size measurement hierarchy are refined
> below (§0), and a corrected experiment, T5b, is required before any Family
> B / full generator work. Everything below §0 is the ORIGINAL pre-
> registration and remains in force except where §0 explicitly narrows it.

## 0. Refined primary causal question and measurement hierarchy (frozen after T5)

**Primary causal question (supersedes the informal framing used through T5):**

> At fixed Ni volume fraction and fixed neck-insensitive particle-size
> measure, does changing the lower tail of the neck-width distribution
> independently improve retained Ni percolation under redox-like damage?

**Measurement hierarchy for "particle size," in order of authority:**

1. **Generator-known primary particle geometry**, where available (e.g. the
   exact sphere radius used to rasterize a structure) — the ground-truth
   control, immune to segmentation artifacts.
2. **Neck-insensitive c-PSD particle-size measure** (`cmlib.particles.cpsd_r50max`)
   — the primary MEASURED control when generator ground truth is not directly
   comparable (e.g. once particles are no longer literal spheres).
3. **Watershed/SNOW chamber size** — reported always, but **diagnostic only,
   never the sole gate**. It is expected to move under neck widening (this
   is the R1 mechanism the whole study is about); a gate built on it would be
   circular.

Any gate phrased as "Φ_Ni within X% and particle size within Y%" (§7-8 below,
and the original EXECUTION_SPEC §1 gates) is read, from this amendment
forward, as requiring (1) or (2) to hold — never (3) alone.

## 0a. T5b — the corrected experiment (required before Family B)

Full spec in `out/next/EXECUTION_SPEC.md` and implementation in
`cmlib/synth.py` / `scripts/next/t5b_coupling_experiment.py`. Summary of what
changed from T5:

- **Mass-conservative widening.** Neck widening no longer just adds Ni
  voxels. Total Ni volume is held near-fixed by compensating elsewhere in the
  Ni phase (sphere-radius shrinkage in this lattice scaffold). Every voxel
  change (added / removed / net residual) is recorded, never silently
  absorbed.
- **Percentile-targeted widening.** Replaces "widen the bottom 20%" (which
  cannot move measured p10 by construction — see the T5 report) with a
  max-clip rule (`new_width = max(old_width, T)`) whose threshold T is found
  by bisection against the ACTUAL MEASURED (SNOW) p10, not the intended
  construction value.
- **Uniform compensated widening kept as a secondary/exploratory positive
  control**, not the primary axis.

## 0b. Acceptance criteria for T5b → Family B (frozen, supersedes §7-8 below for this checkpoint)

For **at least one** lower-tail (percentile-targeted, mass-conservative)
design point:

- measured neck p10 increases by **≥1.5×**, preferably 2×;
- **Φ_Ni deviation ≤2% relative if possible, and never more than 5% relative**
  (this is a hard ceiling, not a target — do not silently relax it; a
  deviation beyond it is reported as a failed point, not rounded down);
- **c-PSD mean particle-size deviation within ±5%** (this, not watershed, is
  the size gate per §0's hierarchy);
- watershed size deviation is reported but is **not** a pass/fail condition;
- **p50 remains much more stable than p10** (the whole point of "lower-tail"
  engineering — if p50 moves nearly as much as p10, the design is not
  actually tail-selective and should be reported as such);
- initial `P_span` remains intact (no collapse of pristine percolation);
- no large number of disconnected Ni fragments is created (`n_clusters`
  tracked, reported).

**Decision tree (frozen):**
- All criteria met on ≥1 lower-tail point → proceed to Family B.
- Lower-tail fails but uniform compensated widening meets the equivalent
  criteria → **stop and report.** This is not itself a green light for Family
  B — it requires an explicit amendment reframing the primary axis from
  lower-tail engineering to whole-distribution shifting before any further
  code is written.
- Neither works → stop, prepare a Path-B-style limitation memo. Do not claim
  physical impossibility; report "not testable within this synthetic
  framework, tested two constructions, both failed for stated reasons."

No damage models, Family C, or large sweeps are built until T5b's outcome is
reviewed against this checkpoint.

## 0c. Amendments after T5b review (frozen 2026-08-10) — mandatory before Family B

T5b is accepted as a real partial success. Decision-tree branch confirmed:
proceed toward Family B, WITH the following mandatory amendments. Family C
and damage models remain out of scope until the Family B disordered pilot
(§0d) is reviewed.

**A. Primary target envelope.** Use base, ~1.5×, and ~2.0× measured neck-p10
ratio as the primary Family B contrast. 2.5× is a stress/limit case only —
never part of the main causal claim (T5b showed the compensation mechanism
breaks down there for some base-distribution shapes; see §0b amendment
below).

**B. Per-seed gating — no aggregate-mean pass criterion is allowed.**
(T5b's own aggregate-mean check silently passed a point that failed on
per-seed inspection — see the corrected T5b report. This is now
structurally forbidden, not just avoided by care.) Each seed is evaluated
individually against ALL of:

1. achieved measured p10 ratio ≥ the target minimum;
2. Φ_Ni deviation ≤2% relative (target), ≤5% relative (hard ceiling);
3. c-PSD size deviation within ±5%;
4. p50 remains much more stable than p10 — **p50 ratio ≤1.15** (a concrete
   number, not the qualitative "much more stable" used through T5b);
5. `n_nodes` ≥95% of the base value;
6. initial Ni percolation (`P_span`) remains intact;
7. no large disconnected-Ni-fragment population is created.

A target ratio is **feasible** if ≥4 of 5 seeds pass all seven. The per-seed
table is always reported in full, pass or fail — an aggregate summary is
never presented without it standing behind it.

**C. Measured achieved p10 ratio is the analysis variable, not the nominal
target.** The nominal target is a generation label only. All plots and
causal contrasts use the achieved value; nominal and achieved are always
reported side by side.

**D. Family B moves beyond the cubic lattice.** T5b's lattice was a valid
diagnostic scaffold, not the main-claim structure. Family B uses a disordered
packing or perturbed lattice, target coordination ~6–8 (report the achieved
value honestly — a scoping choice may fall short, see §0d), initial spanning
Ni connectivity required, particle diameters ~24–32 voxels, lower-tail neck
widths resolved to ≥3 voxels (preferably 4). Pilot at 128³–160³-class before
any 192³-class run.

**E. Base neck-width distribution is designed deliberately, not left to a
plain uniform draw.** Target population p50/p10 ≈3.0–4.0 (broadly consistent
with the real anodes' own neck/particle-size ratios). **Validity criterion,
applied BEFORE widening, never after seeing outcomes: p50/p10 ≥2.5 for that
seed's specific draw.** Seeds failing this are rejected and REDRAWN
(rejection sampling), with every rejected attempt and its reason logged, not
discarded silently. A plain discrete-uniform distribution cannot satisfy both
the ≥2.5 spread and the ≥3-4 voxel floor simultaneously within a physically
sensible neck-width range (checked numerically — see §0d); the base
distribution is therefore a MIXTURE (a "normal" neck population plus a
genuinely narrow-tail subpopulation), which is also more physically
realistic than a uniform draw. An optional narrow-spread STRESS distribution
(deliberately low p50/p10, to reproduce the T5b seed=1-style failure mode on
purpose) may be kept for robustness demonstration but is never used for the
primary causal claim.

## 0d. Family B disordered pilot (required before Family C / damage models)

Implementation: `cmlib/synth.py` (extended, not rewritten — T5b's
mass-conservation and percentile-targeting machinery is geometry-agnostic and
is reused as-is), `scripts/next/familyB_pilot.py`. Scope: 5 seeds, {base,
1.5×, 2.0×}, disordered/jittered structures, full per-seed diagnostic table
(§0b criteria), pass/fail per seed and per target ratio. **Stop and report
this pilot before building ternary YSZ/pore placement or any damage model.**

> **PILOT RESULT (2026-08-10):** run, reported in
> `out/next/familyB_pilot_report.md`. Base-distribution validity check
> passed on the first attempt for all 5 seeds (no T5b-style anomaly
> recurred). Structural quality across all 10 widened structures was
> essentially exact: Φ_Ni deviation −0.0% to +0.03%, c-PSD deviation −0.78%
> to 0.00%, p50 ratio **exactly 1.00 at every single point**. Under the
> strict per-seed gate: nominal 1.5× is NOT feasible (1/5 — a quantization
> effect from the smaller 144-pair population; every failing seed is
> otherwise near-perfect), nominal 2.0× IS feasible (4/5, one seed narrowly
> short at 1.94×, reported as a strict fail per the frozen rule, not rounded
> up). Achieved coordination (mean_degree ~4.3) fell short of the 6-8
> target, reported honestly — a scoping choice (lattice-adjacency topology
> over distance-threshold random packing) traded off for lower disconnection
> risk in a first pilot. No YSZ/pore placement, damage model, Family C, or
> larger sweep has been built; awaiting review.

> **T5b RESULT (2026-08-10):** run, reported in
> `out/next/t5b_coupling_decision_report.md`. Corrected verdict (the script's
> own aggregate-mean check over-reported one point as passing — see the
> report for the sign-cancellation diagnosis): **lower-tail passes cleanly at
> target ratio 1.5×** (4/5 seeds essentially exact: 0.0% c-PSD deviation,
> <0.3% Φ_Ni deviation, p50 completely flat; 1/5 seed marginally over the
> c-PSD ceiling, fully diagnosed, not unexplained noise). Ratio 2.0 fails
> (one seed's c-PSD swings -31%). Ratio 2.5's apparent pass does not survive
> per-seed inspection (4 seeds mildly over ceiling, 1 seed severely so,
> averaging to a false pass). Uniform (whole-distribution) compensated
> widening fails mass conservation at every tested ratio, increasingly
> badly (+7% to +32% mean, up to +88% individual). This is the "lower-tail
> criteria met at ≥1 point" branch of §0b's decision tree — indicated next
> step is Family B — but NOT yet started; awaiting review, per the two
> caveats in the T5b report (achievable envelope is ~1.5-2.0x not 2.5x; one
> seed's failure mode is fully diagnosed and traces to the base neck-width
> distribution's own p50/p10 spread, which a Family-B generator should
> address by design, not rediscover by accident).

## 0e. Amendments after E0 spike review (frozen 2026-08-10) — mandatory before platform v2

E0 (`out/spike/e0_vertical_slice_report.md`) ran cleanly (270 damage runs, not
360 — an arithmetic error in the report, caught on external review and
corrected in place) and returned an honest **inconclusive** result: retained
P_span showed zero differentiation by achieved p10 ratio at any of 6
intensities tested (values were always exactly 1.0 or exactly 0.0, never
intermediate). Per E0's own rules this is not Path A, not Path B, and not
grounds to skip platform qualification. Review of the platform-v2 targets
that followed E0 surfaced real errors, corrected here.

**Coordination target corrected — it was ungrounded and backwards.** §0d's
"target coordination ~6–8" had no derivation attached anywhere in this
document. Pulling the actual measured values
(`out/phase4/phase4c_metrics_per_anode_8.0um.csv`, `mean_degree_mean`) from
the real Holzer/Pecho SNOW-extracted Ni networks:

| | fine | medium | coarse |
|---|---|---|---|
| measured mean_degree | 3.65 ± 0.24 | 3.41 ± 0.28 | 3.55 ± 0.59 |

The real anodes' own coordination is **~3.4–3.65**, which the Family B pilot
(achieved 4.3) already exceeds. **§0d's "~6–8" target is superseded: platform
v2 targets mean_degree ≈3.5–4.5**, matching the real data, not a generic
random-packing prior.

**p50/p10 target tightened, now precisely grounded.** Measured p50/p10 for
the real anodes: fine 4.28, medium 3.29, coarse 3.08. §0c amendment E's
"≈3.0–4.0" is superseded: **target 3.0–4.3**, matching the measured range
exactly.

**Domain-size / particle-size caveat, stated explicitly (not previously
flagged).** p10/diameter ratio varies ~2× across the real anodes: fine
0.0606, medium 0.1086, coarse 0.1197 (computed from the real neck-p10 and
watershed particle-diameter figures already in `REPORT.md`). Resolving a
p10-width neck at ≥3 voxels requires particle diameter ≥3/ratio: **~50
voxels for fine-like ratios, ~26–28 for medium/coarse-like.** §0d's "~24–32
voxel" particle-size target under-resolves fine-like neck geometry by
roughly 2×. Platform v2 must either (a) state explicitly that it is scoped
to medium/coarse-like neck geometry only, or (b) size particles for the
fine-like case (~50–66 voxels), which pushes the domain well past the
192³-class ceiling in §0d. This choice must be made explicitly before
platform v2 is built, not left implicit.

**Working hypothesis, stated falsifiably (E0's "single-bottleneck" reading
was descriptive, not a registered, falsifiable claim — fixed here).**

> H: in a low-coordination (~4), small (64-particle) network, retained
> percolation is dominated by a single structure-wide bottleneck, so damage
> severs the domain as an abrupt, all-or-nothing event rather than as a
> function of the lower-tail neck distribution specifically.

**Falsification condition:** run the same design (multiple intensities
bracketing the transition, multiple seeds) on platform v2. If **any**
intensity shows retained P_span values that are neither uniformly ~1.0 nor
uniformly ~0.0 across the p10-ratio groups, H is falsified and E0's null is
attributed to low coordination, not to the underlying mechanism. **If the
null repeats at platform v2, that does NOT confirm H** — it only means
platform v2 was insufficient too. Resolving that ambiguity, if it arises,
requires tracking transition sharpness across ≥2 platform scales and testing
whether it narrows systematically (the same finite-size-scaling logic
already used successfully on the site-percolation threshold in Phase 0 of
the real-data study) — not declaring H true by default.

**Damage-intensity bisection procedure for the platform-v2 calibration
pilot, frozen now, before it is run.** E0's fixed shared grid
({2,5,6,7,8,10}) is replaced by a per-(structure, damage-seed) bisection on
integer `n_rounds` in bracket `[1, 20]`, mirroring the max-clip threshold
search already validated in Family B: evaluate retained P_span at the
midpoint, narrow the bracket by whether it is still 1.0 or already 0.0, stop
when bracket width ≤1 round. Record `(n_lo, n_hi)` per combination. The
analysis is a regression of bracket midpoint against achieved p10 ratio
across structures — this directly tests the causal question, rather than
hoping a shared coarse grid straddles every structure's transition
simultaneously (which E0 showed is not guaranteed: platform v2's
coordination differs from E0's, so the {2,5,6,7,8,10} grid does not carry
over and must not be reused as-is).

**New pre-registered interpretation branch, added to §F (T5b/Family-B
acceptance) and to be carried into the platform-v2 decoupling experiment's
own pre-registration:** E0's rules named "clean improvement," "no effect,"
and "unstable/intensity-dependent," but never named the mechanistically
plausible opposite. **If retained P_span decreases monotonically (or shows a
clear negative trend) with achieved p10 ratio, report this as a genuine
INVERSE finding — not "no effect," not "unstable."** Leading candidate
mechanism to name, not assume confirmed: mass-conservative widening
compensates by shrinking particle radius, which thins particle *bodies*,
which could make erosion cut through particles rather than necks —
inverting the naive expectation. One pilot cannot confirm this mechanism;
name it as the hypothesis to test, not the conclusion.

**D4 re-validation on platform v2 — required, not assumed to carry over.**
The same three checks E0 ran (bit-identical reconstruction against whatever
platform-v2 generator output is recorded, Ni-untouched-by-YSZ/pore-placement
check, metric-coherence/no-NaN check) must be re-run on the new geometry.
The erosion/expansion mechanism's behaviour is geometry-dependent (more
alternative paths at higher coordination), so E0's intensity range does not
predict platform v2's — this is exactly why the bisection-per-structure
procedure above replaces a fixed shared grid rather than porting E0's
numbers forward.

## 0f. Platform v2 Ni generator qualification (2026-08-10) — RESULT

Run per §0e's corrected targets. Full report:
`out/platform_v2/qualification_report.md`, design memo:
`out/platform_v2/design_memo.md`.

> **Gate P2-A (composition/topology): PASS.** 5 base seeds, Φ_Ni mean=0.2502
> (range 0.2474–0.2561), SNOW mean_degree mean=4.193 (range 4.184–4.206, all
> in [3.5,4.5]), single connected cluster and intact P_span for all 5. **The
> open question is resolved: plain 6-connectivity lattice adjacency already
> satisfies the coordination target — no topology modification (face-diagonal
> bonds etc.) was built or needed.** New generator function:
> `cmlib.synth.platform_v2_lattice_geometry` (asymmetric nlat_z=6/nlat_xy=4,
> R=12.1vox, pitch=32vox — decoupling Φ_Ni control (R/pitch) from topology,
> per §0c's finding that these were conflated in the Family B pilot).
>
> **Gate P2-B (base neck distribution): validity floor PASS, target range
> mean MISS.** All 5 seeds clear p50/p10≥2.5 (range 2.50–3.00), zero
> rejections. Mean (2.90) falls short of the 3.0–4.3 target — reported, not
> rounded up or retuned.
>
> **Gate P2-C (lower-tail decoupling): NOT FEASIBLE at either tested ratio
> (0/5 at intermediate ~1.45x nominal, 1/5 at high 2.0x), for a traced,
> non-arbitrary reason.** Mass conservation and p50 stability are excellent
> (Φ_Ni dev ≤0.02%, p50 ratio exactly 1.00, no node loss) — the sole failure
> mode is c-PSD, which fails by −9% to −12% (vs the ±5% ceiling) for 4/5 seeds
> at the high ratio. **Traced with data, not asserted:** the mass-conservative
> radius shrink needed at high ratio (12.10→11.3–11.5 voxels, ~6-7%) is
> harmless at the Family B pilot's larger base radius (R=14) but large enough
> relative to platform v2's smaller base radius (R=12.1, driven by the
> Φ_Ni=0.250 target) to bleed directly into c-PSD (base 414–452nm → high-ratio
> 396–402nm, matching the radius shrink in direction and magnitude). This is
> the item-6 headroom risk, confirmed empirically — a MILDER, more systematic
> manifestation (4/5 seeds, not one outlier; no seed near the r_lo floor) than
> T5b's catastrophic seed=1 failure, not a recurrence of the same failure
> mode.
>
> **No parameter was retuned after seeing this result** (radius, neck range,
> c-PSD tolerance all unchanged from their frozen values). YSZ/pore placement,
> D4 damage on platform v2, Family C, and real-dataset calibration remain
> unbuilt, per the explicit stop condition. Gate P2-C's failure is a decision
> point for the next review, not something resolved here.

## 1. Claim boundaries (what this study can and cannot conclude)

- **Scope: redox-like damage to Ni/YSZ microstructures only.** No claim about
  long-term thermal/electrochemical Ni migration, no claim this covers all
  SOEC degradation modes, no universal SOEC lifetime claim.
- **No claim of physical impossibility.** If the synthetic generator cannot
  decouple neck size from particle size (gate G1 fails), the conclusion is
  "not testable within this synthetic framework," never "physically
  impossible."
- **No p-values or correlation coefficients from the n=3 real anodes.** That
  restriction, inherited from the original falsification study
  ([[soec-connectivity-margin-study]]), still applies and is never lifted.
  The real data is a QUALITATIVE anchor only (§4).
- **Seed-level statistics on SYNTHETIC structures are legitimate** (effect
  sizes, standard deviations, CIs across seeds) — this is genuine replication,
  unlike the n=3 real-anode case. The two are never pooled into one statistic.
- **Medium vs coarse remains unresolved in the real data** (published
  retention 0.916 vs 0.924, an 0.8-point gap smaller than measurement
  variability) and is never treated as resolved or used as a calibration
  target for anything finer than "roughly comparable."
- **The one robust real-data fact used for calibration** (§4) is that the fine
  (finest-grained, smallest-particle, narrowest-neck) real anode retains Ni
  percolation WORST. This is the only real-data ordering fact this study
  leans on.
- **Damage models are phenomenological, not first-principles.** No claim that
  D1-D4 mechanistically reproduce the physics of NiO/Ni redox volume change;
  they are calibrated to reproduce one qualitative real-data fact (above) and
  are validated on that basis alone.

## 2. Primary and secondary outcome metrics

- **Primary:** retained `P_span` (degraded / pristine), the strict
  face-to-face spanning definition, matching the metric the original
  falsification study used as its headline outcome.
- **Secondary:** retained `P_reach` (reachable-from-either-face, the
  like-for-like analogue of the published *P*), retained TPB density.
- All three computed via `cmlib.api.compute_percolation` /
  `compute_tpb`, the same functions used throughout Phase 0-6 here, so
  synthetic and real-data numbers are never computed by two different code
  paths.

## 3. Damage-model hierarchy (D4/D2/D1), frozen before Phase 2 is built

| model | status | role |
|---|---|---|
| **D4** (morphological redox surrogate: dilation/erosion cycle + removal of disconnected Ni) | **REQUIRED** | primary evidence model — neck sensitivity EMERGES from geometry, not asserted |
| **D2** (random Ni removal, matched total loss) | **REQUIRED** | baseline / null-effect control |
| **D1** (neck-width-targeted removal) | **SANITY / UPPER-BOUND ONLY** | not evidence on its own (see rule below) |
| D3 (betweenness/current-proxy-targeted damage) | optional, if implementable | exploratory only |

**Rule (frozen):** an effect that appears under D1 but NOT under D4 and NOT
under D2 is scored **NEGATIVE**, regardless of its size under D1. Rationale:
D1 ("preferentially destroy narrow necks, then observe narrow-neck-poor
structures survive better") is close to tautological — of the four models it
is the one most likely to manufacture the hypothesized effect by construction.
D4 is the physically-motivated model where the effect, if real, has to emerge
from geometry rather than be assumed into the damage rule; D2 is the null
control that must NOT show the effect if a genuine geometric predictor is
present but WOULD show it if the effect is really just "structures with less
total Ni survive better regardless of geometry" (a volume-fraction confound,
not a neck-geometry effect).

## 4. How the real Holzer/Pecho data is used (qualitative anchor, not a fit target)

- Used ONLY to check that at least one damage model reproduces the ordering
  "fine-like structures lose Ni percolation worse than coarse-like structures"
  (gate G2 in the execution spec).
- NOT used to fit exact retention values, redox-cycle counts, or damage
  magnitudes — fitting three numbers with an arbitrary damage-intensity
  parameter is not a meaningful validation and is explicitly disallowed.
- If a damage model also happens to reproduce "fine retains TPB relatively
  better than coarse" (the second, TPB-side real-data fact from the original
  study), that is reported as a bonus consistency check, not required for
  gate G2.

## 5. Damage-intensity selection procedure (frozen BEFORE Phase 3 results are seen)

Damage intensity can manufacture or erase any effect: mild damage saturates
retained `P_span` near 1.0 for every structure; severe damage saturates it
near 0 for every structure. Either extreme erases the between-structure
signal the experiment is looking for.

**Procedure:**
1. During Phase 2 calibration (Family A structures only, never Family B),
   sweep damage intensity for each required model (D4, D2) and record mean
   retained `P_span` per structure family as a function of intensity.
2. Select the SINGLE intensity, per damage model, at which the **Family A
   mean retained `P_span` falls in [0.5, 0.8]**. If no single intensity value
   achieves this for all three Family-A analogues simultaneously, use the
   intensity that minimizes the spread of family-mean retained `P_span`
   around the [0.5, 0.8] window (documented, not silently chosen).
3. **Freeze that intensity** before generating or damaging any Family B
   structure. It is the PRIMARY comparison point for gate G3.
4. Report the full intensity sweep regardless (never only the frozen point),
   so a reader can see whether the conclusion is an artifact of the specific
   intensity chosen.

## 6. Noise floor (frozen definition)

**Noise floor** = pooled within-structure standard deviation of retained
`P_span` across **≥5 damage seeds**, at FIXED generator parameters and FIXED
(frozen, per §5) damage intensity. Computed separately per damage model.

## 7. Positive-effect criteria for gate G3 (Family B decoupling experiment)

ALL of the following, at the frozen intensity from §5, must hold for a
**positive (Path A)** result:

1. **Effect size:** difference in retained `P_span` between the
   low-neck-p10 and high-neck-p10 ends of Family B ≥
   `max(0.05 absolute, 3 × noise_floor)` (§6).
2. **Monotonicity:** the trend across Family B (ordered by pristine neck p10)
   is monotonic or near-monotonic (no more than one non-monotonic adjacent
   pair, and that pair's reversal must be smaller than the noise floor).
3. **Survives an alternative damage model:** the effect (criteria 1-2) must
   hold under **at least one of {D4, D2}** — NOT D1 alone (§3's frozen rule).
4. **TPB cost bound:** the neck-engineered (or backbone, Family C) design must
   not reduce INITIAL (pristine) TPB density by more than 20% relative to the
   best-TPB baseline at the same Ni volume fraction, unless the result is
   explicitly framed as a tradeoff curve rather than a strict improvement.

Failing any of 1-4 at the frozen intensity → **negative (Path B)** result.
Passing all four → proceed to Phase 4 robustness checks (seeds, intensity
sweep, alternative damage model, ROI/domain size, SNOW `r_max` sensitivity)
before any Path-A claim is written up; a candidate effect that does not
survive Phase 4 is downgraded to "suggestive, not robust," never reported as
established.

## 8. What would make the Family B design itself untestable (gate G1)

Per `out/next/EXECUTION_SPEC.md` §8 Q1/Q2, decided by the T5 coupling
experiment (run immediately after this document is frozen, before any
generator code is written):

- If selective lower-tail neck widening cannot move neck p10 by ≥1.5×
  (target 2×) while holding `d_watershed_volwt` within ±5%, but CAN do so
  while holding `d_cPSD_r50max` within ±5%, both are reported and the primary
  claim is scoped to "neck-vs-cPSD-size decoupling," not
  "neck-vs-watershed-size decoupling" — a narrower but still valid claim.
- If neither size measure can be held within ±5% at the required neck p10
  swing, this is reported as **"the decoupling hypothesis cannot be tested
  within this synthetic framework"** (never "physically impossible") and the
  study proceeds directly to a Path-B-style limitation memo without a Family
  B experiment.

## 9. Sign-off

This document is treated as frozen from this point forward. Any deviation
during Phase 2/3 execution must be logged in
`out/next/phase2_calibration_report.md` / `out/next/phase3_experiment_report.md`
with an explicit "deviation from preregistration" note and reason — not
silently applied.

## 0f. Amendments after platform-v2 qualification review (frozen 2026-08-10)

**1. Gate P2-C is split into two independently-tracked gates.**

- **P2-C1 — mass-conservation / generator self-consistency: PASS.**
  Φ_Ni deviation ≤0.02%; net voxel residuals ±44 on ~30,000 moved; no seed
  near the radius floor (>92% headroom); radius shrink −5.17% to −6.79% at
  high ratio; p50 ratio exactly 1.00; node counts stable (97–101 vs base
  97–98); P_span intact (1.000) throughout.
- **P2-C2 — measured image-based size comparability: DEFERRED** to the
  real-data calibration phase. `cpsd_r50max` failed convergence
  (non-monotone; no resolution with both seeds <0.5 pp); raw EDT is not
  equivalent to either local thickness or generator radius (~3× compressed);
  local thickness is itself unconverged; opening granulometry would require
  the same 5-point ladder cost; and there is **no Platform-v2 consumer** for
  this metric. Carried forward as a named open obligation, not dropped.

**2. `generator_radius_deviation` language, qualified.** It is exact ground
truth **for the input sphere-radius parameter only** — not for what an
image-based thickness measurement of the final rasterized structure would
read. Neck material near junctions changes what a thickness metric sees,
especially as necks widen. See `out/platform_v2/design_memo.md`.

**3. Internal generator-radius guardrail (model-sanity limit, NOT derived
from real ROI variability).**

| | threshold | ≈ sphere-volume loss |
|---|---|---|
| target | shrink ≤ **7%** | ~20% |
| hard review | shrink ≤ **10%** | ~27% |

Beyond that the intervention risks becoming *particle-body thinning* rather
than *neck-tail widening*. Current structures (−5.17% to −6.79%) pass; this
guardrail **has not yet constrained anything**. If a future design point
needs to exceed 7%, that is a trigger to revisit the justification — not a
checkbox to tick.

**4. Causal-interpretation obligation, pre-registered NOW (before any
retention result exists).** The high-ratio intervention couples two changes:
lower-tail neck widening **and** ~5–7% primary sphere shrinkage. Any future
retention-benefit finding tied to achieved p10 ratio **must not be attributed
to neck widening alone** until the radius-shrink confound is addressed by one
of:

- sensitivity analysis across achieved p10 ratio *and* radius deviation;
- a **matched-shrink control** structure with no lower-tail widening, Φ_Ni
  preserved by placing the compensated Ni volume outside the lower tail;
- a widened structure with reduced radius shrink, if feasible.

Sensitivity analysis alone is **necessary but may not be sufficient** if
radius deviation and achieved p10 ratio remain collinear. If the effect
cannot be separated from radius shrinkage, the claim is downgraded to:
*"mass-conservative lower-tail neck widening with accompanying particle-body
shrinkage affects retention."*

**5. Scientific-question language updated.** The intervention is described as
*"mass-conservative lower-tail neck widening at fixed Ni loading, accompanied
by a known ~5–7% primary sphere shrinkage."* The phrase **"fixed particle
size" is not to be used** unless a validated size condition backs it — which,
per P2-C2, it currently does not.

## 0g. Frozen before any widened structure is inspected under damage (2026-08-10)

Committed **before** the p10-group damage experiment is run. Nothing in the
p10-group experiment may proceed until this section is in the repository.

### 1. Damage-seed averaging requirement

The base-only bisection found **within-structure damage-seed variance
comparable to across-structure variance** — transition midpoints of
8.5 / 9.5 / 8.5 occurred *within a single structure*. The damage process is
itself stochastic, so a retention comparison built on one damage seed per
structure would conflate structural effect with damage-process noise.

Therefore, binding on the p10-group experiment:

- **Minimum 3 damage seeds per structure; 5 if compute remains cheap.**
- **If a branch decision would depend on a group difference smaller than 1.0
  damage round, the comparison must be re-run or completed to 5 damage seeds
  before it is interpreted.**
- **No single-damage-seed comparison per structure is permitted.**
- Damage seeds must be **independent of structure seeds** and recorded
  explicitly.
- Required reporting:
  - per-damage-seed transition midpoints;
  - per-structure mean transition midpoint;
  - within-structure damage-seed variance;
  - across-structure variance;
  - group means with seed-level spread.

### 2. Status of D4 parameters `p_erode = 0.35`, `expand_vox = 1`

Recorded explicitly so their provenance is never mistaken for derivation:

- Both **originated in the E0 spike**.
- Neither was **re-derived from first principles** for Platform v2.
- Both are now **deliberately frozen** as part of the D4 operator definition
  for the Platform-v2 damage experiment.
- **The intensity variable expected to shift with geometry is `n_rounds`**,
  not `p_erode` or `expand_vox`.
- The Platform-v2 base-only bisection result — transition at **8.77 ± 0.59
  with no bracket expansion required on any of 15 runs** — is accepted as
  evidence that this frozen parameterization produces a **resolvable**
  transition on the new geometry.
- **Changing `p_erode` or `expand_vox` defines a different damage model and
  requires a new amendment.**

**Additional obligation.** If the p10-group experiment yields a **positive or
weak-positive** branch, a **damage-parameter sensitivity check is required
before any causal claim**: at minimum one alternative `p_erode`, or one
alternative expansion/erosion budget, demonstrating the qualitative result is
not specific to the exact E0 parameterization. A positive result under a
single frozen parameterization is not a causal claim.

### 3. TPB magnitude caveat — corrected

The earlier wording "8–10× real" **understated the low end of the real
range** and is superseded by:

> Platform v2 TPB density after minimal YSZ/pore placement is
> **19.4–20.5 µm⁻²**. This is roughly **7.3–19.2× the real-anode TPB range of
> 1.07–2.65 µm⁻²** (ratio bounds computed as 19.4/2.65 and 20.5/1.07),
> depending on which real anode and which TPB convention is used. This is
> expected from the smoothed random placement field with a small correlation
> length, and is **acceptable for internal percolation/damage work**.
> **TPB magnitude is not yet a real-data-comparable quantity.** Real-data TPB
> comparability belongs with the deferred calibration-phase obligations
> (see P2-C2, §0f).
