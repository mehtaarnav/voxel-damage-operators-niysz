# Path B memo — Ni-YSZ electrode connectivity margin: falsification of the
# lower-tail neck hypothesis, and a mechanism-specific null

**Status: final project memo.** Written after the p10-group causal experiment
and its null autopsy were accepted as closed. No new experiments were run to
produce this document; every number below is cross-referenced to a committed
CSV or report (see `REPRODUCIBILITY_MANIFEST.md`).

**Date:** 2026-08-10. **Scope of closure:** the p10-group causal question is
closed for Platform v2 + frozen D4. Deliberately not run, and not part of this
project: further generator variants, alternative damage models, Family C,
real-data calibration, opening granulometry, a third image-based size metric,
and instruction item 6.

---

## Executive summary

The project asked whether a *pristine-state connectivity margin* — the width of
the narrow tail of the Ni neck-size distribution — predicts, and then causally
governs, how much Ni percolation an SOEC/SOFC electrode retains after
redox-like degradation.

It produced two results, both negative, and the second one with a mechanism.

1. **On real data (Holzer/Pecho, three anodes), pristine connectivity-margin
   metrics do not beat coarseness.** The metrics that reproduce the outcome
   ordering are the ones that are monotone in coarseness and therefore
   indistinguishable from mean particle size at n = 3; the metrics that are
   genuinely graph-theoretic (min-cut conductance, effective conductance) match
   no outcome ordering at all. Separately, Ni volume fraction and pristine TPB
   density rank the anodes **exactly backwards** for retained Ni percolation.

2. **On a purpose-built synthetic platform where neck statistics are decoupled
   from coarseness, widening the lower tail has no effect on the percolation-loss
   threshold.** Structures differing in 67,803 voxels, with a 2× difference in
   measured neck p10, fail percolation at *identical* damage intensity under all
   five damage seeds. The null autopsy explains why: under this damage model
   collapse is **surface-area-mediated global thinning**, and the voxels whose
   removal breaks spanning are one-voxel-deep surface material — not the
   lower-tail neck population. The intervention is **bypassed, not erased**.

The design implication is the general one: lower-tail neck widening is not a
universal percolation-retention strategy, and a degradation-mitigation design
rule is only meaningful once matched to the degradation mechanism that actually
operates.

---

## A. Original predictor falsification on the Holzer/Pecho dataset

**Source:** `REPORT.md`, `out/phase6/phase6_comparison_table.csv`,
`out/phase6/phase6_rankings.csv`, `out/phase5/phase5_retention.csv`.

### A.1 The outcome

Retained Ni percolation on full stacks (0.48–1.11 Gvoxel each, no sub-sampling),
degraded ÷ pristine:

| anode | P_span retained | P_reach retained | published *P* retained |
|---|---|---|---|
| fine | **0.680** | **0.857** | **0.821** |
| medium | 0.855 | 0.979 | 0.916 |
| coarse | 0.947 | 0.942 | 0.924 |

Our absolute values track the published values closely (`P_reach` within
0.1–7.4 % on all six stacks).

**The single robust outcome fact is that the fine anode retains Ni percolation
worst.** Medium vs coarse flips with the definition — `P_span` puts coarse
first, `P_reach` puts medium first, and the published gap is 0.8 points.
**Medium vs coarse is unresolved and is treated as such throughout.**

This bounds how much any predictor "match" is worth: if the only reliable
outcome fact is *fine is last*, a predictor need only put fine last, and 2 of
the 6 possible orderings do that by chance (33 %).

### A.2 The predictors split into coarseness proxies and failed graph metrics

| predictor | ordering (best first) | matches retained Ni percolation? |
|---|---|---|
| 10th-percentile neck width | coarse > medium > fine | yes |
| λ₂ raw (weighted Laplacian) | coarse > medium > fine | yes |
| λ₂ normalised | coarse > medium > fine | yes |
| min-cut, face-to-face | medium > fine > coarse | **no** |
| effective conductance | medium > fine > coarse | **no** |

Three qualifications carry the falsification:

1. **λ₂ raw is measuring node count, not connectivity** — measured, not
   asserted (`phase4e_lambda2_scaling.py`). Pooled over 21 ROIs, λ₂ ∝ N^(−0.867).
   Multiplying by node count collapses a **7.65×** between-anode spread to
   **1.24×** (λ₂·N = 232 / 200 / 247). Normalising removes the separation
   entirely: 0.00094 / 0.00210 / 0.00230 with overlapping error bars.
2. **Only neck-width p10 separates the anodes with non-overlapping error bars**
   (69.5 ± 4.6 / 156.9 ± 43.0 / 205.3 ± 45.0 nm) — and it produces *exactly* the
   mean-particle-size ordering (1148 / 1445 / 1715 nm, watershed,
   volume-weighted) and the SNOW chamber-diameter ordering (671 / 901 / 1082 nm).
   Neck width and particle size are both monotone in coarseness; **on three
   samples they are indistinguishable by construction.**
3. The ranking of that one well-separated metric is **not fully stable** to the
   watershed marker parameter: at `r_max = 2`, medium and coarse swap (184 vs
   179 nm). The swap occurs in exactly the pair the outcome also fails to
   resolve.

**Conclusion.** The hypothesis's claim is that a connectivity margin predicts
retention *better than* mean particle size. This dataset cannot support it: the
connectivity metric that works reproduces the particle-size ranking, and the
connectivity metrics that are *not* reducible to coarseness are the ones that
fail. **Fine-is-worst is the only robust outcome fact; medium vs coarse is
unresolved in predictor and outcome alike.**

### A.3 The clean side result: volume fraction and TPB are anti-predictive

Ni volume fraction orders the anodes fine (0.322) > medium (0.250) > coarse
(0.229) — i.e. predicts the fine anode is *most* robust. It is in fact *least*
robust. Pristine TPB density (3.624 / 2.109 / 1.473 µm⁻²) does the same. Both
are **anti-correlated** with retained Ni percolation, not merely uninformative.

**Retained TPB runs the opposite way to retained percolation:**

| anode | TPB retained (this work) | TPB retained (published) |
|---|---|---|
| fine | **0.799** | **0.743** |
| medium | 0.746 | 0.586 |
| coarse | 0.590 | 0.608 |

Fine is *best* at retaining TPB and *worst* at retaining Ni percolation. Any
claim about "degradation resistance" must therefore name which property it
means. None of the connectivity-margin metrics match the retained-TPB ordering;
Ni volume fraction and pristine TPB density both match it exactly.

### A.4 Statistical standing

n = 3 is a directional check. **No correlation coefficient and no p-value is
reported anywhere in this work, and none should be computed from three points**
— the Phase 6 script enforces this.

---

## B. Successor causal question

**Source:** `out/next/preregistration.md`, `out/platform_v2/design_memo.md`,
`out/platform_v2/qualification_report.md`,
`out/platform_v2/ternary_d4_bisection_report.md`.

### B.1 Why a synthetic platform

The Phase-6 failure mode is a **confound, not a null**: coarseness varies
monotonically across the three real anodes and drags every candidate metric with
it. Separating a connectivity margin from a coarseness proxy requires structures
where the two **disagree** — matched on mean particle size and Ni loading,
differing only in the lower tail of the neck distribution. Three samples that
differ in everything at once cannot do this however carefully the graph metrics
are computed.

### B.2 Platform v2 design

Jittered lattice of overlapping spheres with explicit cylindrical necks, so neck
width is a *generator input* rather than an inferred measurement.

| parameter | value |
|---|---|
| domain | 161 × 168 × 168 voxels @ 20 nm = 3.22 × 3.36 × 3.36 µm |
| lattice | 6 (z) × 4 × 4 = 96 particles, pitch 32 vox, jitter 0.15·pitch |
| sphere radius (base) | R = 12.1 vox (242 nm) |
| bonds | 224 nearest-neighbour pairs (6-connectivity adjacency only) |
| neck-width mixture | 20 % weak in [4,6] vox, 80 % normal in [12,20] vox, min p50/p10 ratio 2.5 |
| Φ_Ni target | 0.250 (medium anode); YSZ share of remainder 0.5173 (medium anode's own value) |
| geometry seed | 999; structure seeds 0–4 |

### B.3 Qualification

- **P2-A composition/topology — PASS.** Φ_Ni mean 0.2502 (range 0.2474–0.2561);
  SNOW mean degree 4.184–4.206 (target band 3.5–4.5); single connected Ni
  cluster and P_span = 1.000 on all 5 seeds. **Plain nearest-neighbour lattice
  adjacency already meets the coordination target — no topology modification was
  built.** The earlier "coordination too low, needs face-diagonal bonds" framing
  was an artifact of an ungrounded 6–8 target and was discarded rather than
  built upon.
- **P2-B base neck distribution — PASS with a note.** p50/p10 mean 2.90 (range
  2.50–3.00): every seed clears the pre-registered 2.5 floor, but the mean sits
  *below* the stated 3.0–4.3 target. Reported, not rounded up. Base p10 resolved
  at 5.2 vox (preferred ≥4).
- **YSZ/pore placement — PASS.** Ni phase untouched (asserted by array equality);
  Φ_Ni shifts by ≤3.4 × 10⁻⁵; Φ_YSZ ≈ 0.385–0.389 vs the medium anode's 0.388;
  P_span unchanged to 1 × 10⁻¹².
- **D4 re-validation on the new geometry — ALL CHECKS PASS** (reconstruction
  bit-identical 5/5; YSZ untouched 5/5; removed Ni → pore 5/5; all metrics
  finite 5/5).

### B.4 The intervention: mass-conservative lower-tail neck widening

Necks below a threshold T are max-clipped up to T; sphere radius is then solved
so total Ni volume is conserved. This is genuinely tail-selective and genuinely
mass-conservative:

| mode | intended T (vox) | achieved p10 ratio | neck p10 (nm) | neck p50 (nm) |
|---|---|---|---|---|
| base | — | 1.000 | 120.0 | 320.0 |
| lower-tail | 8.5 | 1.333 | 160.0 | 320.0 |
| lower-tail | 11.0 | 2.000 | 240.0 | 320.0 |

**p10 moves 120 → 160 → 240 nm while p50 stays pinned at 320 nm** — the
tail-selective signature. Mass conservation residuals are ±44 voxels on ~30,000
moved (Φ_Ni deviation ≤0.02 %).

**Achieved-ratio rungs, not nominal targets.** The nominal 1.45× rung is not
achievable: attainable p10 values are set by the base mixture's own discrete
value grid, and all 5 seeds landed on 1.33×. Achieved ratio is the scientific
variable, so the rung was retargeted rather than forced.

### B.5 Known radius-shrink confound (recorded, pre-registered before any outcome)

Mass conservation is achieved by shrinking sphere radius **5.17–6.79 % at 2.00×
(1.70–2.13 % at 1.33×)**. The intervention therefore couples two changes:
lower-tail neck widening **and** ~5–7 % primary sphere shrinkage. Pre-registered
§0f/4 forbids attributing any retention benefit to neck widening alone without a
matched-shrink control or a shrink/ratio sensitivity analysis. An internal
guardrail (target ≤7 % shrink, hard review at ≤10 %) was set; it has never
constrained anything.

The phrase **"fixed particle size" is not used**, because no validated
image-based size condition backs it.

### B.6 Image-based size comparability: deferred, not dropped

Gate P2-C was split. **P2-C1 (mass-conservation / generator self-consistency):
PASS.** **P2-C2 (measured image-based size comparability): DEFERRED** to a
calibration phase that this project does not enter. Reasons, all measured:
`cpsd_r50max` failed a convergence test (non-monotone, no resolution with both
seeds <0.5 pp); raw EDT is not equivalent to local thickness or generator radius
(~3× compressed); local thickness is itself unconverged; opening granulometry
would incur the same 5-point ladder cost; and **there is no Platform-v2 consumer
for the metric.** `generator_radius_deviation` is exact ground truth **for the
input sphere-radius parameter only**, not for what an image-based thickness
measurement of the rasterized structure would read.

---

## C. p10-group damage experiment

**Source:** `out/platform_v2/p10_group_report.md`,
`out/platform_v2/p10_group_experiment.csv`,
`scripts/platform_v2/p10_group_experiment.py`.

### C.1 Frozen D4 model

Morphological redox surrogate, promoted unchanged from the E0 spike
(`cmlib/damage.py`): per round, isotropic Ni dilation by `expand_vox`
(6-connectivity, restricted so YSZ is never overwritten), then independent
removal of *surface* voxels with probability `p_erode`, then removal of Ni
disconnected from the spanning backbone.

- **`p_erode = 0.35`, `expand_vox = 1` — FROZEN.** Provenance recorded
  explicitly (§0g/2): both originated in the E0 spike and **neither was
  re-derived from first principles for Platform v2.** They are frozen as part of
  the operator definition. **The intensity variable is `n_rounds`.** Changing
  `p_erode` or `expand_vox` defines a different damage model and requires an
  amendment — not taken.
- Resolvability was verified before use: the base-only bisection transitions at
  **8.77 ± 0.59** with no bracket expansion on any of 15 runs, and sits *inside*
  E0's unresolved 8→10 gap — vindicating per-structure bisection over a shared
  fixed grid.

### C.2 Design

15 structures (5 structure seeds × {base, 1.33×, 2.00× achieved p10}) × **5
damage seeds (200–204, independent of structure seeds)** = **75 bisections**,
each narrowed to a bracket of width 1. Integer `n_rounds`, initial bracket
[1,20], expand-only.

**Primary outcome: the bisection transition midpoint** — the damage intensity at
which the Ni network loses spanning percolation.

Damage-seed averaging was pre-registered (§0g/1) *before* any widened structure
was inspected: minimum 3 seeds, 5 if cheap; no single-seed comparison permitted;
any branch decision resting on a group difference below 1.0 round must be
completed to 5 seeds first.

### C.3 Result: no resolvable effect

| group | achieved p10 | group mean transition | spread of structure means | n structures |
|---|---|---|---|---|
| base | 1.00× | **8.54** | 0.089 | 5 |
| lower-tail | 1.33× | **8.50** | 0.000 | 5 |
| lower-tail | 2.00× | **8.50** | 0.000 | 5 |

Differences vs base: **−0.040 rounds at 1.33×, −0.040 rounds at 2.00×.** Both
are ~25× below the pre-registered 1.0-round interpretability threshold, both are
**smaller than one standard error** (≈0.09 rounds at sd ≈ 0.2, n = 5), and both
are **negative** — not even a sub-threshold hint of the hypothesised benefit.

**74 of 75 bisections returned a midpoint of exactly 8.5.** The single exception
is base structure seed 0, one damage seed, returning 9.5. **All 10 widened
structures returned 8.5 with zero variance across both groups and all 5 damage
seeds.**

### C.4 Classification of the result (frozen wording)

- **No resolvable effect.** Under the pre-registered E0 rules this is the
  "no effect" branch: **inconclusive as a universal scientific hypothesis**, and
  a **clean model-specific null** under Platform v2 + frozen D4.
- **Not Path A.** **Not a weak positive.** **Not a pipeline artifact.**
- The §0g/2 damage-parameter sensitivity obligation applies to a positive or
  weak-positive branch and is **not triggered**.
- The §0f/4 radius-shrink confound is **moot** — there is no effect to attribute
  — but is recorded, and the matched-shrink control is correspondingly not
  triggered.

### C.5 Two integrity checks that the null had to survive

**(i) Damage-seed variance did not reproduce.** §0g/1 was written because the
base-only bisection showed within-structure damage-seed variance comparable to
across-structure variance (8.5/9.5/8.5 within one structure).

| run | damage seeds | n | midpoints | mean | sd |
|---|---|---|---|---|---|
| base-only | 100–102 | 15 | 7.5×1, 8.5×9, 9.5×5 | 8.767 | **0.594** |
| p10-group, base | 200–204 | 25 | 8.5×24, 9.5×1 | 8.540 | **0.200** |

Same structures, same operator, same frozen parameters — only the seeds differ.
The earlier sd of 0.594 was driven by the particular seeds 100–102. The
requirement was still correct to impose (cheap, and it is what exposed this),
but its stated motivation is weaker than the base-only run suggested. **The net
effect is to sharpen the null**, not weaken it: with sd ≈ 0.2 the effect is
absent at the resolution this design can reach, not marginally missed.

**(ii) Provenance of the widened structures was verified before the null was
accepted.** This was a live risk: the script rebuilds each structure from an
`intended_T_vox` looked up in `qualification_run.csv`, so a lookup or
seed-indexing slip would silently produce base geometry while still labelling
the row `lower_tail`.

- Thresholds non-null and distinct (8.5 vs 11.0); `achieved_p10_ratio` has
  `nunique = 1` per group at 1.000 / 1.333 / 2.000 — no base row leaked in.
- Raw neck arrays, structure seed 0: `array_equal` **False**, **46 of 224 necks
  changed**, min 4 → 11 vox while p50 stays 15.0 vox — the max-clip acted on the
  tail only.
- Raw Ni masks: `array_equal` **False**, **67,803 voxels differ** (1.49 % of the
  domain), base-only 33,909 vs high-only 33,894 — the mass-conservation
  signature (~34k removed from bodies, ~34k added at necks, net −15).

**Both checks pass. The zero-variance pattern is a real property of the damage
response, not an absence of signal in the pipeline.**

### C.6 What the design could not have seen

The bisection returns integer-bracketed midpoints, so the finest attainable
resolution is 1.0 round. **An effect smaller than one erosion round is invisible
to this design by construction.** Detecting a sub-round effect would need a
finer-grained outcome (fractional damage intensity via `p_erode`, or retained
conductance rather than binary percolation loss) — which defines a different
damage model and is not taken here.

---

## D. Null autopsy

**Source:** `out/platform_v2/null_autopsy_report.md`, `null_autopsy.csv`,
`null_autopsy_localization.csv`, `scripts/platform_v2/null_autopsy.py`.

Secondary and explanatory only; no D4 parameter was changed and no new primary
outcome was introduced. Diagnostic question: *at n_rounds = 8 (the last spanning
state), do the widened structures still retain a measurable lower-tail geometric
advantage over base?*

### D1 — Remaining Ni volume: no advantage

| group | achieved p10 | vox @ n=0 | vox @ n=8 (sd) | retained @ n=8 (sd) | retained @ n=9 |
|---|---|---|---|---|---|
| base | 1.000× | 1,137,007 | 476,038 (14,674) | 0.4186 (0.0082) | 0.3454 |
| lower-tail | 1.333× | 1,137,025 | 471,078 (14,431) | 0.4142 (0.0079) | 0.3401 |
| lower-tail | 2.000× | 1,137,013 | 470,650 (11,539) | 0.4139 (0.0049) | 0.3373 |

Starting volumes match to ~18 voxels in 1.14 M. By n = 8 all three groups have
lost ~58 % of their Ni and sit within **0.5 percentage points** of each other,
with group differences (0.0044, 0.0047) **smaller than the within-group sd**
(0.005–0.008). The widened groups retain very slightly *less* volume.

### D2 — Lower-tail thickness proxy: **degenerate, excluded from evidence**

The EDT percentile over *all* Ni voxels returns **20.0 nm p10 for every group at
every damage state, including n = 0**, where the structures are known to differ
in 67,803 voxels with neck p10 of 120 vs 240 nm.

20.0 nm is exactly **one voxel**. In any solid, well over 10 % of voxels lie
within one voxel of a surface, so the 10th percentile pins to one voxel
regardless of internal geometry; p25 shows the same quantization one step up
(40 nm at n = 0) and collapses to one voxel after damage because erosion raises
the surface-to-volume ratio.

**Recorded status: degenerate / non-informative. D2 must not be read as
"the intervention was erased" and must not be used as evidence for or against
the intervention in either direction. It will not be repaired or re-run as part
of this project.** A genuine lower-tail diagnostic would have to be restricted
to the *neck* population (e.g. SNOW throat inscribed diameters on the damaged
network) — that is a new metric, not run here.

### D3 — Backbone size: no advantage

At n = 8, **P_span = 1.000 for every group** — the entire remaining Ni network
is a single spanning cluster in all cases, so there is no backbone-size
advantage available to detect. Widened backbones are marginally *smaller*,
within noise. At n = 9, base retains a residual P_span of 0.04 (one seed still
spanning) while both widened groups sit at exactly 0.00 — the same −0.04-round
difference from the primary result, seen from the other side, and in the
**unfavourable** direction.

### D4loc — Failure-step localization: the decisive diagnostic

Voxels removed between n = 8 and n = 9, EDT measured in the n = 8 mask
(structure seed 0, damage seed 200):

| group | n removed | frac of n=8 | removed p10 | p50 | p90 | removed mean | *all* n=8 p50 | *all* n=8 p90 | *all* n=8 mean |
|---|---|---|---|---|---|---|---|---|---|
| base | 83,020 | 0.17 | 20.0 | 20.0 | 20.0 | **20.00** | 28.28 | 74.83 | 38.36 |
| high 2.00× | 86,174 | 0.18 | 20.0 | 20.0 | 20.0 | **20.02** | 28.28 | 72.11 | 36.95 |

The voxels whose removal breaks spanning are **essentially exclusively surface
voxels at exactly one voxel depth** (p10 = p50 = p90 = 20.0 nm) against a
background distribution with median 28.28 nm and p90 ≈ 72–75 nm. The removed set
is statistically indistinguishable between base and 2.00× (20.00 vs 20.02 nm),
and the removed *fraction* is nearly identical (0.17 vs 0.18).

That is D4's erosion operator by construction: each round strips surface voxels
with probability `p_erode`. **Collapse occurs when cumulative uniform surface
stripping has thinned the network globally — not when a specific narrow neck is
severed.**

### Conclusion of the autopsy

Under the pre-registered rules this is **branch A — "the intervention does not
survive to the collapse boundary under D4."** That wording is recorded as the
formal branch. The refinement the rules did not anticipate, and which the
diagnostics support directly:

> **The intervention is bypassed, not erased.** D1 shows all groups converge to
> within 0.5 pp of the same retained volume by n = 8; D4loc shows why. The
> failure mode is uniform surface erosion, which operates on total surface area
> and is indifferent to *where* material sits in the neck-width distribution. A
> lower-tail intervention has nothing to act on if collapse is not mediated by
> the lower tail.

**Therefore, under this D4 parameterization, lower-tail neck widening has no
mechanistic lever on the percolation-loss threshold.** This explains the null
without concluding that lower-tail necks are irrelevant in all damage models —
it identifies the specific property of *this* model responsible.

**Branch C (hidden continuous benefit) is NOT supported.** Every continuous
diagnostic — retained volume, spanning-cluster size, spanning fraction of
original — points marginally the *wrong* way, all within noise. There is no
sub-round benefit hiding beneath the integer outcome. This is recorded
explicitly so that the null is not later re-read as an under-powered positive.

The conclusion rests on **D1, D3 and D4loc**, which are mutually consistent.

---

## E. Design implication

1. **Lower-tail neck widening is not a universal percolation-retention
   strategy.** Doubling the p10 neck width at conserved Ni loading changed the
   percolation-loss threshold by less than the experiment's one-round
   resolution, and by less than one standard error, in the unfavourable
   direction.

2. **Degradation-mechanism match is essential.** The value of a microstructural
   intervention is conditional on the mechanism that actually destroys the
   network. A lower-tail intervention presumes collapse is bottleneck-mediated.
   Under surface-area-mediated thinning that premise is false, and the
   intervention has no lever regardless of how large it is made.

3. **Under surface-area-mediated thinning, the actionable levers are the ones
   that act on surface area and total path redundancy** — total Ni loading, the
   surface-to-volume ratio of the conducting phase, and backbone architecture —
   **not the lower tail of the neck distribution.** This is a statement about
   which lever the mechanism exposes, not a measurement of those alternative
   levers, none of which were tested here.

4. **The corollary for predictor work (Part A):** a pristine-state geometric
   margin can only predict retention if the degradation operator is sensitive to
   that margin. The real-data falsification and the synthetic null are the same
   finding seen from two directions.

---

## F. Limitations

1. **D4 is a phenomenological damage model.** Dilation + stochastic surface
   erosion + island removal is a redox *surrogate*, not a validated redox
   simulation. Its parameters `p_erode = 0.35` and `expand_vox = 1` originated in
   the E0 spike and were **never re-derived from first principles**; they were
   frozen, not justified. The null is specific to this operator at this
   parameterization.
2. **Outcome resolution is one damage round.** Integer-bracketed bisection cannot
   see sub-round effects by construction (§C.6).
3. **Synthetic TPB magnitude is not real-data comparable.** Platform v2 TPB after
   minimal YSZ/pore placement is 19.4–20.5 µm⁻², roughly **7.3–19.2×** the real
   anode range of 1.07–2.65 µm⁻² (bounds 19.4/2.65 and 20.5/1.07). This follows
   from a smoothed random placement field with a 3-voxel correlation length and
   is acceptable for internal percolation/damage work only.
4. **Image-based size comparability is deferred (P2-C2).** No converged
   image-based particle-size metric backs the claim that base and widened
   structures are size-matched as an imaging measurement would read them. The
   generator radius is exact ground truth for the *input parameter* only.
5. **The radius-shrink confound is moot but recorded.** 2.00× structures carry a
   5.17–6.79 % sphere-radius shrink (1.70–2.13 % at 1.33×). With no effect to
   attribute, the confound does not arise and the matched-shrink control was not
   triggered — but any future positive result on this platform must address it
   first.
6. **D2 is invalid for this question and contributes nothing.** One of the three
   requested diagnostics returned no usable information; the conclusion rests on
   the other three.
7. **No real redox calibration was performed.** The damage model is not tuned or
   validated against the Holzer/Pecho degraded stacks; the real data was used
   only as a qualitative anchor, never as a fit target.
8. **No Family C / hierarchical architecture test.** Only a jittered-lattice
   single-scale topology was tested. Whether a hierarchical or bimodal backbone
   behaves differently under the same operator is untested.
9. **Real-data limitations carry through from Part A**, unchanged: the coarse
   anode is sub-REV for network metrics; medium vs coarse retention is within
   noise; TPB ground truth is digitized from a bar chart (self-consistent to
   0.7 pp); post-redox medium and coarse stacks are 40 % voxel-anisotropic; the
   conductance weighting (throat area / throat length) is a modelling choice.
10. **A `porespy.snow2` chunked-watershed bug was active** during the Phase 3/4
    SNOW extractions. It was found, fixed for new work (`cmlib/pnm.py` now
    defaults to serial extraction), and impact-assessed — **no conclusion
    changes** (`IMPACT_NOTE_porespy_parallel_bug.md`).
11. **n = 3 on the real data.** Directional check only; no correlation
    coefficients or p-values are computed anywhere, by design.

---

## G. Future work — explicitly outside this project

None of the following was run, and none is required to support any conclusion
above. Listed so that the closure is a scoping decision on the record rather
than an omission.

1. **Neck-selective damage models.** The direct test of the mechanism-match
   thesis: an operator whose removal probability depends on local thickness or
   on throat identity rather than on surface membership. The prediction this
   project generates is that lower-tail widening *should* show an effect there —
   and the pre-registered D1 rule already flags that such a result would be
   near-tautological on its own, so it would need pairing with an independent
   mechanism.
2. **Hierarchical / bimodal backbone architectures (Family C).** Whether a
   redundant coarse backbone with a fine secondary network resists
   surface-mediated thinning better than a single-scale network at equal loading.
3. **Real-data calibration of damage and size metrics.** Tuning D4 (or a
   successor) against the Holzer/Pecho pristine→degraded pairs, which would also
   discharge the deferred P2-C2 size-comparability obligation.
4. **TPB realism.** Replacing the smoothed random YSZ/pore field with a
   correlated placement at micron-scale domain size to bring synthetic TPB into
   the real range.
5. **Targeted bottleneck interventions.** Widening necks selected by betweenness
   or current-carrying importance rather than by width percentile — a different
   intervention on a different feature of the network.
6. **A finer-grained outcome variable.** Retained effective conductance, or
   fractional damage intensity, to resolve sub-round effects. Requires a damage
   model amendment.

---

*End of memo. Companion documents: `PAPER_OUTLINE.md`,
`REPRODUCIBILITY_MANIFEST.md`.*
