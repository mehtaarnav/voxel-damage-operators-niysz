# Pre-registration — CSFM: cumulative-strain intergranular fracture of YSZ

**Frozen 2026-08-11 before implementation.** Ordinal only; no fitting of any
parameter to match the data. `cmlib/damage.py`, `cmlib/synth.py` frozen.

## Hypothesis
Coarse-anode YSZ percolation loss is caused by redox-cycling intergranular
fracture concentrated at YSZ necks already structurally weak in the pristine
state, inherited from incomplete sintering of the coarse starting powder.

---

## Three corrections to the brief, made BEFORE committing

### C1 — Step 7 as written is impossible, not merely hard *(fatal)*

The brief asks to *"compare the set of fractured throats to the set of throats
that are disconnected in the degraded coarse stack."*

**Pristine and degraded are different specimens** (Rx38 vs Rx41-3). There is no
throat-to-throat correspondence between them: they are different pieces of
material, sectioned independently. A per-throat overlap set cannot be formed,
and decision rule (a) as written cannot be evaluated on any dataset we hold.
This is the same confound that rendered H1 untestable (`0bccf20`).

**Re-specified step 7, frozen:** compare **distributions and aggregate
outcomes**, never throat identities.
- **(7a) Aggregate:** does CSFM applied to pristine coarse drive YSZ `P_span`
  to the degraded value (**0.000**), and pristine fine/medium toward theirs
  (**0.958 / 0.865**), at a *common* cycle count?
- **(7b) Distributional:** does the predicted-fractured throat population match
  the *size distribution* of the excess disconnected YSZ in the degraded stack
  (two-sample comparison of inscribed-diameter distributions, reported
  ordinally)?
- **(7c) Positional:** is predicted fracture spatially concentrated, as observed
  YSZ loss is, rather than uniform?

### C2 — Cycle count must be an intensity variable, not a fixed 8

Fixing 8 cycles makes a failure uninterpretable: wrong mechanism, or wrong
cycle count? We do not know how many cycles these specimens saw.
**Frozen:** cycle count is the bisection intensity over [1, 20], expand-only,
width ≤ 1 — identical to every other operator in this project. "8 cycles" is
recorded as the brief's nominal value and reported alongside, not used as the
test.

### C3 — The matched-fragility control cannot equalise structure on real data

The brief proposes *"scaling the strain threshold per anode so that pristine
fragility is matched."* Scaling the threshold changes **the operator**, not the
structure, so it does not reproduce the R3 control's logic (which equalised
pristine fragility *in the material* and inverted a +2.87-round result to
−0.27). Real structures cannot be restructured.

**Re-specified control, frozen — two parts, both required:**
- **(C3a) Within-anode:** across regions of interest *within* one anode, does
  predicted fracture concentrate on throats that were already weak pristine?
  This tests the mechanism without any cross-anode comparison.
- **(C3b) Covariate:** across anodes, partial Spearman of predicted YSZ damage
  against coarseness **controlling for pristine YSZ fragility** (1 − P_span:
  0.0011 / 0.0120 / 0.0754). If the association vanishes under the partial, the
  effect is pristine-loading.

---

## Operator (frozen)

1. Local Ni volume fraction field on the **pristine** stack.
2. Ni→NiO expansion applied as uniform volumetric strain, **linear strain
   ε₀ = 0.0065** (mid-point of the 0.55–0.8 % range in the brief) ⚠[unverified].
3. Gaussian smoothing of the strain field, **σ = 3 voxels** (frozen now).
4. Accumulation over `n` cycles with **permanent residual factor r = 0.25** per
   cycle (frozen now): cumulative ε(n) = ε₀ · [1 + r·(n−1)].
5. Per YSZ throat, cumulative strain = mean of the smoothed field over the
   throat's two adjacent grains.
6. Fracture if cumulative strain > **threshold τ = 0.010** (frozen now,
   **never fitted**).
7. Score per C1 above.

**Recorded limitation:** the strain field is computed from the *pristine* Ni
distribution, whereas real strain acted on material that was itself evolving.
This is a deliberate phenomenological simplification, stated up front.

## Decision rules

- **SUPPORT** if (7a) reproduces the aggregate ordering with coarse driven
  furthest, **and** C3a shows concentration on pristine-weak throats,
  **and** C3b's partial association survives.
- **FALSIFY** if predicted damage shows no aggregate ordering, or if C3b's
  association vanishes under the partial (pristine-loading).
- **UNTESTABLE** if the YSZ throat graph cannot be extracted on the degraded
  stacks, or if no cycle count in [1, 20] produces a resolvable transition.

## Sub-resolution caveat (stated, not worked around)

Real YSZ throat p10 ≈ 40 nm = **2 voxels**, below the ≥4-voxel resolution floor
this project holds itself to. Throat-targeted fracture therefore rests partly on
sub-resolution features, and any positive result inherits that limitation.

## Constraints
Ordinal only. **τ, σ, ε₀ and r are frozen above and are never adjusted to
improve agreement.** Report support and falsification outcomes alike.
