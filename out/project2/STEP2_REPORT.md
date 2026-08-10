# Project 2 — Step 2: **FAILS**. C1 unresolved, C2 pristine-state loaded, C3 fails.

Run 2026-08-11. Under `PREREGISTRATION_V2_1.md` (`e62f30b`), Step 1
re-qualification approved (`52aa519`). Operators in `cmlib/damage2.py`;
`cmlib/damage.py` and `cmlib/synth.py` unmodified. **O4 not implemented.**

Data: `step2_r1_scope.csv`, `step2_n_secondary.json`, `step2_r2_tpb_baseline.csv`,
`step2_O1_main.csv`, `step2_O3_main.csv`, `step2_O2_ceiling_check.csv`,
`step2_p_control.json`, `step2_O3_control.csv`.

---

## Verdict

| criterion | result |
|---|---|
| **C1** — fine loses Ni percolation first, ≥1.0 round | **FAIL — UNRESOLVED** (0.07 and 0.00 rounds) |
| **C2** — coarse loses YSZ percolation before fine, ≥1.0 round | **PASS in the main arm (+2.87 rounds)** — but the R3 control shows it is **primarily pristine-state loaded** |
| **C3** — C1 + C2 together, plus retained TPB ordered fine > medium > coarse | **FAIL** (requires C1; TPB retention is also ~164× too low) |
| **O2** | **No transition within the bracket.** Stop-and-report per the boundary rule. |

**Step 2 fails.** Per the instruction this is reported as a negative result, with
the missing feature identified in §7. **No operator was added and no parameter
was tuned in response to any result.**

---

## 1. Deliverable 1 — operators and unit tests

`cmlib/damage2.py`, tested by `scripts/project2/test_operators.py` — **all tests
pass**. Highlights: O1 is asserted **bit-identical to the frozen `apply_d4`**;
O2 and O3 are verified monotone in `n_rounds` (required for bisection to be
well-posed); O3 is verified to remove neck volume, update the contact graph, and
leave Ni untouched; the TPB estimator is verified against an analytic
single-triple-line case.

**Two bugs found and fixed by those tests, before any science ran:**

1. **TPB periodic wrap.** Gathering the four voxels around an edge with
   `np.roll` wraps at the domain faces and manufactures triple lines. On the
   analytic case it over-counted by exactly **4×**. Fixed by slicing (interior
   edges only). *Numerically this barely moves our structures (fine 27.8 → 27.7),
   because the boundary contribution is small in a genuine three-phase volume —
   but the earlier figure I quoted was produced by the wrapping version.*
2. **O2 cost.** Per-throat full-domain dilations would have dominated the run;
   localised to per-region bounding boxes.

Monotonicity is realised by a single uniform draw per contact with survival
`u < (1-p)^n` — exact marginals, monotone by construction, and O(1) in `n`.

## 2. Deliverable 2 — R1 scoping and frozen `n_secondary`

Frozen from the scoping sweep by the stated rule: **O1 = 10, O2 = 5, O3 = 2**.

**The frozen rule produced a saturating value for O1 and this must be reported,
not fixed.** O1's transitions sit at ≈ 8.6, so retained `P_span` at n = 10 is
**0.0000 for all 45 structure × damage-seed combinations in all three classes**.
The secondary outcome is therefore **uninformative for O1** — it cannot
distinguish "no effect" from "resolution-limited", which was its entire purpose.
The rule selected the midpoint of the shallowest observed bracket (8, 13) → 10,
which lies past the transition. A value of 8 would have discriminated.

**I did not re-pick it.** It was pre-registered before the runs and re-picking
after seeing saturation is exactly the tuning the protocol forbids. The
consequence is that C1's null is reported as a **primary null with no usable
secondary evidence**, rather than as a "resolution-limited trend".

## 3. Deliverable 3 — R2 pristine TPB baseline (frozen, µm⁻²)

| analog | median | seed range | real | ratio |
|---|---|---|---|---|
| fine | **27.708** | 27.107–28.020 | 3.624 | 7.6× |
| medium | **17.228** | 17.167–17.335 | 2.109 | 8.2× |
| coarse | **11.294** | 11.039–11.438 | 1.473 | 7.7× |

Ordering fine > medium > coarse holds; **fine/coarse = 2.45 against a real
2.46**. Absolute magnitude is ~8× high, uniformly. TPB is computed and stable,
so **C3 is not blocked** on measurement grounds.

## 4. Deliverable 5–6 — main arm results, per-seed

### C1 — O1 (Ni surface erosion), 45 bisections, all bracketed cleanly

Per-structure means over 3 damage seeds:

| analog | s0 | s1 | s2 | s3 | s4 | **class mean** | sd |
|---|---|---|---|---|---|---|---|
| fine | 8.50 | 8.83 | 8.17 | 8.83 | 8.83 | **8.633** | 0.516 |
| medium | 8.83 | 8.50 | 8.83 | 8.83 | 8.50 | **8.700** | 0.414 |
| coarse | 8.50 | 8.83 | 8.83 | 8.50 | 8.50 | **8.633** | 0.352 |

Separations: fine vs medium **0.067**, fine vs coarse **0.000**. Both far below
1.0. **C1 is UNRESOLVED — reported as a null, not a trend.**

**Why the null happens — the informative part.** Ni volume loss at the
transition state differs strongly and in the predicted direction:

| analog | Ni volume lost at transition |
|---|---|
| fine | **74.4 %** |
| medium | 70.5 % |
| coarse | **62.8 %** |

**Fine does erode faster** — the 1.26× specific-surface-area advantage is real
and visible. But fine's network **tolerates proportionally more volume loss
before losing spanning**, and the two effects cancel almost exactly. Higher
specific surface area and higher connectivity redundancy are both consequences
of fineness, and they oppose each other in this outcome variable.

### C2 — O3 (YSZ fracture), 45 bisections, all bracketed cleanly

| analog | s0 | s1 | s2 | s3 | s4 | **class mean** | sd |
|---|---|---|---|---|---|---|---|
| fine | 6.50 | 5.83 | 5.83 | 6.17 | 6.17 | **6.100** | 0.507 |
| medium | 4.17 | 3.50 | 3.83 | 3.50 | 3.50 | **3.700** | 0.561 |
| coarse | 3.17 | 3.17 | 2.83 | 3.83 | 3.17 | **3.233** | 0.704 |

**C2 passes in the main arm: fine − coarse = +2.87 rounds**, with non-overlapping
seed ranges (coarse max 3.83 < fine min 5.83).

**Ambiguity assessment (deliverable 7).** medium vs coarse is **0.47 rounds** —
below the 1.0 threshold — and their seed ranges overlap (coarse 2.83–3.83,
medium 3.50–4.17). **The three-level ordering is therefore AMBIGUOUS and is not
claimed.** Only the two-level C2 comparison is resolved. This mirrors the real
data, where medium vs coarse is itself unresolved.

### O2 — no transition within the bracket

P_span = 1.0000 at every scoped intensity, and at n = 20 in **15 of 15
structures**, with **100 % of candidate throats severed** and only **1.8 % Ni
volume loss**. Per the boundary rule: stop and report, bracket not expanded.

**This is a substantive result, not an absence of signal.** Severing *every*
neck in the lower quartile of the throat-size distribution does not disconnect
the Ni network. It is the exact converse of Project 1, which found that
*widening* the lower tail gave no retention benefit because collapse is
surface-mediated rather than bottleneck-mediated. **Two independent experiments,
opposite interventions, same conclusion: the narrow-neck population is not
load-bearing for Ni percolation in this architecture.**

## 5. Deliverable 8 — pristine Q vs transition, and the R3 control

Main arm, pristine YSZ fragility rank-correlates perfectly with the transition:

| analog | pristine Q_YSZ | O3 transition |
|---|---|---|
| fine | 0.0009 | 6.10 |
| medium | 0.0110 | 3.70 |
| coarse | 0.0647 | 3.23 |

**R3 matched-fragility control** (p_control = 0.672 / 0.542 / 0.642, targeting
Q ≈ 0.010; achieved 0.0054 / 0.0105 / 0.0076 — within ±0.004 for 2 of 3 classes,
so the Q-matched arm stands and the fallback was not triggered):

| analog | **main** midpoint | **control** midpoint |
|---|---|---|
| fine | 6.100 | 4.433 |
| medium | 3.700 | 4.300 |
| coarse | 3.233 | **4.700** |

**fine − coarse: main +2.87 rounds → control −0.27 rounds.**

**Coarse-worst vanishes, and in fact inverts.** Per the pre-registered rule:
**C2 in the main arm is primarily pristine-state loaded, not mechanistic.** With
pristine fragility equalised, all three classes fail within 0.4 rounds of each
other — i.e. grain size *per se* confers no YSZ damage tolerance under O3.

This is exactly what R3 was added to detect, and without it the main-arm C2 pass
would have been an over-claim.

## 6. C3 — fails, and the TPB diagnosis is the most informative result

C3 requires C1, which failed. But the TPB data is worth stating on its own:

| analog | pristine TPB | TPB at O1 transition | **synthetic retention** | **real retention** |
|---|---|---|---|---|
| fine | 27.708 | 0.135 | **0.0049** | **0.799** |
| medium | 17.228 | 0.082 | 0.0047 | 0.746 |
| coarse | 11.294 | 0.064 | 0.0057 | 0.590 |

**Synthetic TPB retention is ~164× lower than real.** O1 does not merely fail to
reproduce the TPB divergence — it destroys essentially all TPB (99.5 %) by the
time Ni percolation is lost, whereas the real fine anode *retains 80 % of its
TPB* at the point where it has lost the most Ni percolation. The ordering is
preserved but the magnitude is nowhere close, and the ordering is nearly
meaningless at these retention levels.

## 7. Which feature is missing — the diagnosis required by the FAIL branch

Comparing the synthetic damage progression against the real Step 0 data
identifies one decisive mismatch.

**All three operators conserve or reduce phase volume by *removal*. The real
data says Ni is not simply removed.** From `phase2_volume_fractions.csv`:
Φ_Ni changes by **−0.1006 (fine), −0.0165 (medium), +0.0146 (coarse)** — it
*increases* in the coarse anode. A removal-only operator cannot reproduce a
volume increase, and more importantly cannot reproduce **losing connectivity
while retaining Ni–YSZ contact area**, which is what 80 % TPB retention at 68 %
Ni percolation retention means.

**The missing mechanism is volume-conserving Ni redistribution — agglomeration
/ dewetting — rather than Ni removal.** Under reducing conditions Ni dewets YSZ
and coarsens by surface diffusion: Ni retracts into larger, more compact
particles. That breaks *long-range* connectivity (percolation) while the
surviving Ni–YSZ contact perimeter, where TPB lives, is comparatively preserved.
Uniform surface erosion, by contrast, thins everything at once and eliminates
TPB and percolation together — which is precisely the 164× discrepancy above.

*Literature anchors from memory, requiring verification from source before any
manuscript use: Simwonis, Tietz & Stöver (Solid State Ionics 132, 2000) on Ni
particle growth in Ni-YSZ; Sarantaridis & Atkinson (Fuel Cells 7, 2007) on redox
cycling; NiO→Ni molar volumes ≈ 11.2 vs 6.59 cm³/mol.*

This diagnosis was foreshadowed in `DESIGN_MEMO` §0.2 item 2 ("Any operator that
destroys Ni by removing Ni voxels will move Φ_Ni monotonically down and can
never reproduce this") and is now confirmed quantitatively.

## 8. Deliverable 9 — Go/No-Go recommendation

**No-Go on the independent-operator model as specified. No-Go on proceeding
directly to O4.**

- O4 couples Ni strain to YSZ fracture. But §5 shows YSZ damage tolerance is
  **not** grain-size dependent once pristine fragility is matched, and §6 shows
  the Ni half already fails on TPB. **O4 would inherit both defects.** Running
  it now would produce a C2 that is again pristine-loaded and a C3 that again
  fails on TPB.
- **Recommended next step: a volume-conserving Ni agglomeration operator (O5)**,
  tested against the one signature no operator has yet reproduced — TPB retained
  while Ni percolation is lost. This is the highest-leverage move because it is
  the only one that attacks the actual, measured mismatch.

## 9. Limitations

1. `n_secondary` for O1 saturated (§2) — C1's null has no secondary support.
2. Analogs remain ~3× size-compressed; the recorded bias direction is *against*
   C1, so C1's null is partly confounded with the compression.
3. O3 transitions are early (2–6 rounds), so 1.0-round resolution is a large
   fraction of the dynamic range; the two-level C2 margin (2.87) is comfortable,
   the three-level one is not.
4. Real pre/post are different specimens; all real retention values carry that
   confound.
5. TPB is ~8× real in absolute magnitude in every class.
6. One YSZ placement seed per structure seed.
7. O1's cross-check that YSZ is untouched shows YSZ P_span 0.9745 at the Ni
   transition (unchanged from pristine); O3's shows Ni P_span exactly 1.0 — the
   operators are confirmed phase-isolated.
