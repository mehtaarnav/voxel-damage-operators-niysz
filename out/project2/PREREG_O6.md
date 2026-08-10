# Pre-registration amendment — O6 (reduction-only Ni surface erosion)

**Frozen 2026-08-11, before O6 is written or run.** Amends `PREREG_C1_REAL.md`
(`49b1059`). `cmlib/damage.py` and `cmlib/synth.py` remain frozen; **O1 is
unchanged and is NOT run on real data.** O6 lives in `cmlib/damage2.py`.

## Definition (frozen)

**O6 = O1 with `expand_vox = 0`.** Pure stochastic surface erosion at
`p_erode = 0.35`; largest-component pruning retained; YSZ never touched; no
oxidative dilation step.

**Physics ruling (advisor, adopted):** the dilation step models Ni→NiO
oxidation, a ~70 % volume expansion (NiO ≈ 11.2 vs Ni ≈ 6.59 cm³/mol). The
process being modelled is Ni loss under *reducing* conditions — coarsening,
dissolution, electrochemical removal — which is volume-losing or
volume-conserving, not expanding. On the synthetic platform the step was benign
only because gate G1-c enforced pristine `P_span = 1.000`, leaving nothing to
heal. On real data it bridges genuine disconnections and manufactures TPB. It is
an operator artifact, not a mechanism.

## A — Validity gate, corrected to compare against PRISTINE

All of the following, at every tested `n`:

1. `P_span(n) ≤ P_span(pristine)` — erosion must not increase connectivity.
2. Ni volume loss monotone increasing in `n`.
3. `TPB(n) ≤ TPB(pristine)` — erosion must not manufacture TPB.
4. YSZ untouched (array equality).

**Any failure ⇒ STOP. `p_erode` is NOT adjusted.** The Step-2/A1 check compared
only across damage intensities and could not detect the failure it was written
for; this version compares against pristine, which is where the violation lives.

## B — Pristine disconnection recorded as a standalone finding

| anode | pristine Ni `P_span` | disconnected |
|---|---|---|
| fine | 0.9821 | 1.8 % |
| medium | 0.9713 | 2.9 % |
| coarse | **0.8878** | **11.2 %** |

Coarse is most disconnected, consistent with the YSZ fragility ordering. This is
reported independently of any damage operator, and it establishes that the
synthetic platform's **G1-c (`P_span = 1.000`) was over-constrained relative to
real electrodes.**

## C — Partial correlations in the mechanism test

Because the three classes start from different pristine `P_span`, report **both**:
raw Spearman ρ of transition intensity vs (a) pristine specific Ni surface area
and (b) corrected pristine min-cut fraction; **and partial Spearman for each,
controlling for pristine `P_span`.** Threshold unchanged: |ρ| ≥ 0.6 with
consistent sign under leave-one-out. If both correlate, report both. This
prevents a spurious "rate" result actually driven by coarse starting closer to
failure.

## D — v3 audit re-run with corrected extents BEFORE the bisection

The reversed z/y/x extents make v3's min-cut fractions provisional. The
corrected re-run must complete **before** the C1-real bisection so the
correlation analysis uses final values.

## E — TPB estimator verification on real voxels

Verify the corrected (non-wrapping) estimator on one real ROI and report the
pristine value. If it remains far from the literature (~3.6 µm⁻² fine), record
the discrepancy and **do not gate C3-real on absolute TPB — gate only on the
retention ratio.**

## Unchanged

Hypothesis, decision rules (PRIMARY C1-real, MECHANISM, SECONDARY C3-real),
protocol (3 ROIs/anode min, fine+medium 8 µm, coarse 12 µm, 3 damage seeds,
bisection [1,20] width ≤ 1), and all anti-tuning constraints.
