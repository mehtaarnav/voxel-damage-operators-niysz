# CSFM — **UNTESTABLE.** Both pre-registered untestable conditions fire.

Run 2026-08-11 under `PREREG_CSFM.md` (frozen `ffcb066` **before**
implementation). Frozen and unadjusted: ε₀ = 0.0065, σ = 3 vox, r = 0.25,
τ = 0.010. Deterministic operator — no damage seed. Data: `csfm_results.csv`.

## Verdict

The pre-registration declared CSFM **UNTESTABLE** if *"the YSZ throat graph
cannot be extracted"* or *"no cycle count in [1, 20] produces a resolvable
transition."* **Both conditions fire.** No support and no falsification verdict
is issued.

| anode | YSZ throats | pristine P_span | broken at n=20 | P_span at n=20 | retention | **measured target** |
|---|---|---|---|---|---|---|
| fine | 2,401 | 0.9993 | 25 (1.0 %) | 0.9988 | **0.9995** | 0.958 |
| medium | 506 | 0.9709 | 6 (1.2 %) | 0.9704 | **0.9996** | 0.865 |
| **coarse** | — | — | — | — | — | **0.000** |

**Condition 1 — coarse yields no YSZ network.** The 12 µm coarse region produced
no extractable YSZ throat graph: the YSZ phase does not span that region along
the transport axis, so the watershed extraction returns nothing. **Coarse is the
anode the hypothesis is about**, and its measured YSZ retention of 0.000 is the
entire target signature. Without it there is no hypothesis to test.

**Condition 2 — no resolvable transition.** In fine and medium the operator is
very nearly inert across the whole intensity range. At the maximum cycle count
it breaks 1.0–1.2 % of throats and moves YSZ retention to 0.9995 and 0.9996,
against measured targets of 0.958 and 0.865. The transition the bisection was
built to locate does not exist inside [1, 20].

## Diagnosis — why the frozen operator is inert

Reported as diagnosis, **not as grounds for adjustment**. The strain field is the
Gaussian-smoothed local Ni fraction times ε₀, so it is bounded above by ε₀ =
0.0065 and, after smoothing, the per-throat values peak far lower:

| anode | per-throat strain: min / median / max |
|---|---|
| fine | 0.00000 / 0.00046 / **0.00413** |
| medium | 0.00000 / 0.00041 / **0.00210** |

The frozen threshold **τ = 0.010 exceeds the maximum single-cycle strain by
2.4× (fine) and 4.8× (medium).** Only the cumulative multiplier
1 + r(n − 1) — reaching 5.75 at n = 20 — lifts the most-strained throats over τ
at all, and then only a handful. The median throat never comes within an order
of magnitude of the threshold.

τ, ε₀, σ and r were frozen before implementation and are **not adjusted**, per
the explicit constraint that τ is *"never fitted"*. Re-picking τ now — after
seeing that it sits above the strain distribution — would be fitting the
threshold to produce damage, which is exactly what the pre-registration forbids.

## What this does and does not establish

**It does not falsify the hypothesis.** Coarse-anode YSZ fracture at
pristine-weak necks may well be the mechanism; this run could not examine it,
because the target anode produced no graph and the frozen parameterisation
produced no transition in the two anodes that did.

**It does establish two facts about testability**, both of which would have to be
resolved before any future attempt:

1. **A 12 µm coarse region is too small for the YSZ phase to span.** Coarse YSZ
   is the most fragmented phase in the dataset (pristine full-stack P_span
   0.9246, the lowest of six), and at this region size it does not percolate at
   all. Any CSFM test on coarse needs a substantially larger region, or the full
   stack.
2. **A strain proxy bounded by ε₀ cannot exceed a threshold set at 1.5 ε₀**
   without heavy cumulative amplification. Threshold and strain scale must be
   dimensionally reconciled *before* freezing, not after — which requires
   measuring the strain distribution on one region first and pre-registering τ
   relative to it, as a quantile rather than an absolute.

The second point is a genuine pre-registration design lesson: **we froze an
absolute threshold against a quantity whose scale we had not yet measured.**
That is the same class of error as the saturated secondary-outcome rule recorded
earlier in this project — a rule correct in intent, wrong in the value it
selected.

## Limitations

Single region per anode; deterministic operator, so no seed averaging applies;
strain computed from the pristine Ni distribution while real strain acted on
evolving material; YSZ throat p10 ≈ 2 voxels, below this project's 4-voxel
resolution floor, so throat-targeted fracture rests partly on sub-resolution
features; and the C3 controls were not reached, because they presuppose a
measurable effect to control for.
