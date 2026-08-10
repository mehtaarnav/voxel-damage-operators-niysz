# Project 2 — Kill test result: **PASS**

Run 2026-08-10 under `PREREGISTRATION_V2_1.md` (committed `e62f30b`) §5.
Code: `cmlib/project2.py` (generator), `scripts/project2/kill_test.py`.
Data: `killtest_scope.csv`, `killtest_k3k4.csv`, `p_sinter_calibrated.json`.

**No damage operator was implemented or applied.** O1/O2/O3 do not exist in the
codebase. `cmlib/damage.py` and `cmlib/synth.py` are unmodified.

---

## Verdict

| gate | criterion | result |
|---|---|---|
| **K0** | p=0 must not percolate; p=1 must reach P_span ≥ 0.98 | **PASS** (all 3 analogs) |
| **K1** | \|ΔΦ_YSZ\|/Φ_YSZ ≤ 2 % | **PASS** — worst 0.025 % |
| **K2** | Q_low/Q_high ≥ 10, Q_low ≥ 0.05, Q_high ≤ 0.005 | **PASS** — 324× / 1555× / 1569× |
| **K3** | coarse median P_span ∈ [0.90, 0.95] | **PASS** — 0.9282, 3/3 seeds in bracket |
| **K4** | fine median P_span ≥ 0.995 | **PASS** — 0.9991, 3/3 seeds |
| **K5** | filtered diagnostics (never gating) | reported, direction matches real |
| — | K2 grain-diameter drift ≤ 5 % | **exceeded for medium (6.9 %) and coarse (6.7 %)** — see §6 |

**Option A is viable.** The particulate YSZ generator decouples grain size from
sintering yield, which is precisely what the random-field placement could not do.

---

## 1. Two implementation bugs found and fixed before any result was accepted

Both were mine, both produced result-shaped failures, and both are recorded
because either one would have been reported as a false kill.

**Bug 1 — axis-aligned neck rasterizer.** `cmlib.synth.rasterize` draws a neck
as an axis-aligned bar. That is correct for the Ni generator (6-connected
lattice, every contact axis-aligned) but silently wrong for FCC, whose nearest
neighbours lie along face diagonals: the bar advances along one axis while
staying at the source grain's coordinates in the others, never reaching its
partner. **Symptom: P_span = 0.0000 even at p_sinter = 1.0.** Fixed by
`_add_capsule` / `rasterize_ysz`, a true distance-to-segment capsule valid for
any contact direction. Unit-checked: two diagonal grains give 2 components
unsintered, 1 sintered.

**Bug 2 — boundary filter excluded overhanging grains.** Grain centres outside
the domain were discarded, so no grain covered the extreme slices whenever the
YSZ lattice was incommensurate with the domain depth. At a_ysz = 46 in a
181-deep domain the last centre sat at z = 161 and reached z = 176 against a
last slice of 180 — **spanning was impossible by construction, at every
p_sinter.** Fixed by allowing centres to overhang by one nearest-neighbour
distance and letting the rasterizer clip.

**Consequence for an earlier claim of mine, corrected.** I reported mid-run that
the fine analog failed K0 "by a packing bound". **That was wrong** — it was
Bug 2. With the boundary filled, fine hits Φ_YSZ = 0.4210 (deviation 0.00 %) at
a gap of **1.60 voxels**, against my analytic prediction of 1.61. The analytic
model was right all along; the earlier saturation was an unfilled-boundary
artifact.

**What does survive from that analysis:** the **SC** lattice is genuinely
infeasible for the fine analog, by closed form rather than experiment —
Φ_YSZ_max ≈ f_pack·(1−Φ_Ni) = 0.5236 × 0.678 = **0.355** against a required
**0.421**, i.e. 119 % of its own ceiling. SC was rejected on the derivation and
never tested. FCC (cap 0.502) was selected because at equal gap it yields the
most grains per domain.

---

## 2. Generator as built

| | fine | medium | coarse |
|---|---|---|---|
| lattice | FCC | FCC | FCC |
| cube edge a_ysz (vox) | 40 | 46 | 53 |
| neck width w_ysz (vox) | 8 | 10 | 11 |
| YSZ jitter (frac of nn) | 0.02 | 0.02 | 0.02 |
| grains | 792 | 665 | 500 |
| candidate contacts | 3,536 | 2,942 | 2,152 |
| nominal coordination | 8.93 | 8.85 | 8.61 |
| inter-grain gap at p=0 (vox) | **1.60** | **3.74** | **4.66** |

Every A4 requirement is satisfied and verified rather than asserted: grains are
explicit spheres at recorded centres; contacts are an explicit pair list;
sintered contacts are explicit capsules (connection guaranteed); unsintered
contacts leave a ≥1.6-voxel gap (disconnection guaranteed, confirmed by K0 at
p=0); and Φ_YSZ is held to 0.025 % while topology changes.

---

## 3. K0 — topology sanity

| analog | p=0: P_span | clusters/grains | gap | p=1: P_span |
|---|---|---|---|---|
| fine | **0.0000** | 738 / 792 | 1.60 vox | **0.9992** |
| medium | **0.0000** | 579 / 665 | 3.74 vox | **0.9999** |
| coarse | **0.0000** | 385 / 500 | 4.66 vox | **0.9995** |

At p=0 no analog percolates — connectivity is not coming from accidental raster
overlap. At p=1 all exceed 0.98. **PASS.**

(Cluster counts at p=0 are below grain counts because overhang grains are
clipped to the domain and some merge across the periodic-free boundary region;
the decisive quantity, P_span = 0, is unambiguous.)

## 4. K2 — the lever, and the true kill criterion

| analog | low state | Q_low | high state | Q_high | **ratio** |
|---|---|---|---|---|---|
| fine | p=0.35, P_span 0.8967 | 0.2485 | p=1.00, P_span 0.9992 | 0.00077 | **324×** |
| medium | p=0.35, P_span 0.9411 | 0.1969 | p=1.00, P_span 0.9999 | 0.00013 | **1555×** |
| coarse | p=0.35, P_span 0.8544 | 0.7904 | p=1.00, P_span 0.9995 | 0.00050 | **1569×** |

Required: ≥ 10×, Q_low ≥ 0.05, Q_high ≤ 0.005. **All satisfied with two to three
orders of magnitude of margin.** Grain size is fixed by construction across each
sweep and Φ_YSZ is held to 0.025 %, so the Q movement is topological.

## 5. A5 recalibration, K3 and K4

The generator graph is FCC with measured mean coordination ≈ 8.5–8.8, **not SC
bond percolation**, so the SC-derived priors are superseded exactly as A5
anticipates. One-shot recalibration from the measured Q(p) curve, done **before**
any damage and frozen in `p_sinter_calibrated.json`:

| analog | real Q target | **calibrated p_sinter** | SC prior |
|---|---|---|---|
| fine | 0.00109 | **0.955** | 0.70 |
| medium | 0.01197 | **0.523** | 0.56 |
| coarse | 0.07542 | **0.416** | 0.42 |

Medium and coarse land close to the SC priors; fine is much higher, because at
coordination ~8.8 a far higher yield is needed to reach Q ≈ 0.001.

**K3 — coarse, 3 seeds at p = 0.416:** P_span = 0.92819 / 0.93232 / 0.90815,
median **0.92819**, inside [0.90, 0.95], **3/3 seeds in bracket**. **PASS.**
(Real coarse is 0.9246 — reproduced without being fitted to.)

**K4 — fine, 3 seeds at p = 0.955:** P_span = 0.99869 / 0.99908 / 0.99923,
median **0.99908** ≥ 0.995, **3/3 seeds**. **PASS.**

## 6. The one deviation: K2 grain-diameter drift

Protocol A allows ≤ 5 % mean grain-diameter drift during K2 "unless the Q change
is clearly topological".

| analog | low-state d | high-state d | drift |
|---|---|---|---|
| fine | 530 nm | 520 nm | **1.90 %** ✓ |
| medium | 564 nm | 526 nm | **6.91 %** ✗ |
| coarse | 645 nm | 603 nm | **6.66 %** ✗ |

Cause: holding Φ_YSZ fixed while adding ~2,900 (medium) necks of width 10 voxels
requires the core radius to give volume back.

**I believe the escape clause applies, and I am flagging it for your judgement
rather than assuming it.** The Q change is topological on three independent
counts: (i) at p=0 the phase is fully non-spanning with grains disconnected by a
3.7-voxel gap, and no volume change could produce that; (ii) cluster count falls
monotonically 579 → 167 as contacts are added; (iii) Φ_YSZ is constant to
0.025 % throughout, so no volume is being gained or lost — only redistributed
from cores to necks.

**Note also that drift is a sweep-only artifact.** In use each analog sits at a
single frozen p_sinter, where grain diameter is a single well-defined value
(521 / 559 / 641 nm). Nothing downstream inherits the drift.

**I did not reduce w_ysz to force this under 5 %.** That would be tuning a
parameter to pass a gate, which is the pattern A5 exists to prevent.

## 7. K5 — filtered cluster diagnostics (never gating)

Median filtered cluster count (> 0.1 grain volume) at calibrated p_sinter:
**0 / 6 / 21** for fine / medium / coarse. Rises with coarseness, matching the
real direction (0.005 / 0.098 / 0.266 per Mvoxel). Reported as a diagnostic
only, per A1.

## 8. Forward look — not gates, and not fitted

At the calibrated operating points, the quantities the full Step 1
re-qualification will gate:

| | fine | medium | coarse | required | real |
|---|---|---|---|---|---|
| **Q_YSZ median (G1-i)** | 0.00092 | 0.00997 | 0.07181 | rising, ratio ≥ 10 | 0.00109 / 0.01197 / 0.07542 |
| ratio coarse/fine | — | — | **78×** | ≥ 10 | 69× |
| **mean YSZ EDT (G1-h)** | 57.4 | 64.7 | 73.6 nm | rising | — |
| grain diameter | 521 | 559 | 641 nm | — | — |

**The ordering that Step 1 failed on with the random field is now reproduced.**
This is not a fit: p_sinter was calibrated to Q on a single-seed curve, then the
values above were measured on three fresh seeds — but it is also **not yet the
gate**, which requires 5 seeds per class and the anti-outlier clause.

---

## 9. Limitations

1. **G1-i is not yet passed** — §8 is a 3-seed preview at kill-test settings,
   not the 5-seed qualification.
2. **p_sinter was calibrated against the real Q targets**, which is
   pre-authorized (A5) and done before damage, but it means Q agreement in §8 is
   *targeted*, not predicted. The non-trivial content is that the generator
   *can* hit all three simultaneously at fixed Φ — which the random field could
   not do at any σ.
3. **YSZ grain size is not validated against real data.** No YSZ grain-size
   measurement exists in this study; a_ysz was set to scale with the Ni particle
   diameter and to give a ≥1.5-voxel gap for the fine analog. Only the *ordering*
   is claimed.
4. **YSZ jitter is 0.02**, far more regular than the Ni lattice's 0.15, forced
   by the need to keep the fine analog's 1.60-voxel gap open. The YSZ network is
   therefore more crystalline than real.
5. **The frozen sub-grain limitation still applies verbatim:** "The synthetic
   YSZ generator will model grain-scale backbone connectivity and pristine
   percolation fragility. It will not reproduce the measured sub-grain fragment
   population, fragment-size distribution, or fragment surface area."
6. **Single YSZ placement seed per structure seed**; grain geometry and
   sintering draws are not independently replicated.
7. K2 low/high states were read off a **1-seed** scope sweep; K3/K4 used 3 seeds.

---

## 10. Status

- Kill test: **PASS** on K0, K1, K2, K3, K4; K5 diagnostic consistent.
- One deviation recorded: K2 diameter drift 6.9 % / 6.7 % for medium / coarse
  (§6), escape clause argued but **not self-granted**.
- `p_sinter` frozen at 0.955 / 0.523 / 0.416, recalibrated once, before damage,
  logged in `p_sinter_calibrated.json`. **No further changes permitted.**
- **Not run:** O1, O2, O3; full Step 1 re-qualification (5 seeds/class).

**Awaiting your review before Step 1 re-qualification or any damage modelling.**
