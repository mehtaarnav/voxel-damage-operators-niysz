# C1-real — **FAIL. Definitive scope boundary.**

Run 2026-08-11 under `PREREG_C1_REAL.md` (`49b1059`), `PREREG_O6.md` (`f980562`),
`PREREG_RNI_METRIC.md` (`bfc5d89`). Code: `scripts/project2/c1real_run.py`.
Data: `c1real_results.csv`. O6 frozen (`p_erode = 0.35`, no expansion, uniform
surface-site removal). `cmlib/damage.py`, `cmlib/synth.py` untouched. Amendment D
completed **before** the bisection.

## Result

**Fine is the LAST anode to lose Ni percolation, not the first.**

| threshold | fine | medium | coarse | fine − medium | fine − coarse | verdict |
|---|---|---|---|---|---|---|
| R_Ni < 0.50 | **9.72** | 9.06 | 9.39 | **+0.67** | **+0.33** | **FAIL** |
| R_Ni < 0.10 | **9.72** | 9.28 | 9.39 | **+0.44** | **+0.33** | **FAIL** |

C1-real required fine strictly **lower** than medium and coarse by ≥ 1.0 round at
both thresholds. Fine is **higher** at both, and every separation is also below
1.0 round, so the ordering is simultaneously **wrong in sign and unresolved in
magnitude**.

Against reality this is a complete reversal: real retention runs
**fine 0.680 < medium 0.855 < coarse 0.947** (worst → best). O6 gives
**medium 9.06 < coarse 9.39 < fine 9.72**.

## Mechanism test — neither variable explains anything

| threshold | variable | raw ρ | partial ρ (controlling pristine P_span) | LOO range | signs consistent |
|---|---|---|---|---|---|
| 0.50 | specific Ni surface area | +0.018 | −0.283 | [−0.110, +0.412] | no |
| 0.50 | min-cut fraction | +0.028 | −0.185 | [−0.124, +0.316] | no |
| 0.10 | specific Ni surface area | +0.089 | −0.079 | [−0.026, +0.509] | no |
| 0.10 | min-cut fraction | +0.213 | +0.109 | [+0.026, +0.509] | yes |

**Nothing reaches |ρ| ≥ 0.6.** The largest raw correlation is +0.213. **Neither
the rate hypothesis (surface area) nor the topology hypothesis (min-cut) is
supported.** The question the milestone was built to decide is not decided in
favour of either alternative — both are rejected.

## Amendment D — corrected min-cut fractions (final, not provisional)

| anode | ROI 0 | ROI 1 | ROI 2 | nodes |
|---|---|---|---|---|
| fine | 0.0196 | 0.0174 | 0.0141 | 439–574 |
| medium | 0.0302 | 0.0242 | 0.0221 | 141–186 |
| coarse | **0.0055** | **0.0063** | **0.0067** | **258–299** |

With coarse properly sized at 12 µm (258–299 nodes, versus v3's 48–72), **coarse
is by far the most topologically fragile** — 3–5× lower min-cut fraction than
fine. Yet coarse is the *best* retainer in reality and mid-pack under O6.
**Pristine topological fragility does not predict real degradation.**

## The early signal — fate reported explicitly

The flagged single-ROI n=8 ordering was **coarse 0.577 < fine 0.609 < medium
0.699**. The bisection gives **medium < coarse < fine**. So:

- **The specific ordering was a single-ROI artifact** — coarse is not last.
- **Its qualitative content held**: fine is not first, which is what mattered
  for C1.

## Pristine disconnection (reported regardless of outcome)

Pristine Ni `P_span` per ROI: fine 0.9821 / 0.9754 / 0.9877; medium 0.9713 /
0.9446 / 0.9554; coarse 0.8878 / 0.9528 / 0.9157. **Real electrodes carry
1.2–11.2 % disconnected Ni when pristine.** The synthetic platform's gate G1-c
demanded exactly 1.0000 and was therefore over-constrained relative to every
real electrode measured.

## TPB — diagnostic only, no conclusion

Pristine TPB 4.48–4.66 (fine), 1.87–2.37 (medium), 1.47–1.55 (coarse) µm⁻²;
against literature ~3.62 / 2.11 / 1.47, i.e. **1.0–1.3×** — the corrected
estimator is sound on real data. O6's TPB values are not interpretable because
the operator is volume-losing and surface-area-increasing while real coarsening
is volume-conserving and surface-area-reducing. **C3-real is not addressable
with O6.** No gate, no conclusion.

## SCOPE BOUNDARY

**Voxel-scale stochastic surface-erosion operators of this class do not
reproduce the real Ni percolation-loss ordering.** This was obtained with a
validated operator, a well-posed pruning-invariant metric, on real disordered
microstructures, at properly sized ROIs, with both the rate and topology
predictors measured and neither correlating.

Per the frozen clause: **the stopping is the result.** No replacement operator
is proposed here. The failure is not attributed to O6 — O6 passed its validity
gate, and the two candidate mechanisms were tested directly and both rejected.

The correct next question is **what class of mechanism does** reproduce the
ordering — not which operator to try next.

## Publication position — either outcome was publishable, and this one is

The project has now systematically eliminated: lower-tail neck widening
(Project 1), narrow-neck severing, YSZ-fracture-as-mechanism (pristine-loaded,
shown by the R3 control), volume-conserving redistribution as implemented, and
now surface erosion on real microstructures. Alongside sits a methods
contribution of four artifact classes, each invisible on synthetic structures
and each capable of silently corrupting a published degradation study:
pruning-dependence of connectivity metrics, TPB manufacture by voxel-scale
erosion, over-constrained pristine connectivity, and lattice min-cut planarity.

Claim scope remains bound by Amendment A5 and requires source verification.

## Limitations

3 ROIs per anode (class-level claims are correspondingly weak, though the sign
reversal is consistent across all 9 ROIs); 3 damage seeds; the two thresholds
are not independent because the spanning cluster vanishes discontinuously —
R_Ni crosses 0.50 and 0.10 within one round in most ROIs, which is itself a
finding about how percolation fails; single 8/12 µm ROI scale; graph metrics use
frozen SNOW parameters.
