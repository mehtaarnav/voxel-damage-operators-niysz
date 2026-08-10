# Project 2 — Step 1: analog qualification — **STOP**

Run 2026-08-10, 15 structures, ~4 min total. Code:
`scripts/project2/step1_analog_qualification.py`, σ rule in `cmlib/project2.py`.
Data: `step1_analog_qualification.csv` (all 15 structures, every generation
parameter and seed), `step1_gates.csv`, `step1_class_means.csv`.

**`cmlib/damage.py` and `cmlib/synth.py` were not modified.** A2.3 needed no
edit to either: `add_ysz_pore` already exposes `smooth_sigma_vox`, so the σ rule
went into a new module, `cmlib/project2.py`. No damage operator was
implemented.

---

## Verdict

**STOP. Three gates failed, and one of them is the pre-registered stop
condition.**

| gate | scope | result |
|---|---|---|
| G1-a Φ_Ni ±2 % | per-structure | **14/15 — FAIL** |
| G1-b Φ_YSZ ±2 % | per-structure | 15/15 PASS |
| G1-c Ni percolates, 1 cluster, P_span = 1 | per-structure | 15/15 PASS |
| G1-d YSZ percolates | per-structure | 15/15 PASS |
| G1-e particles ≥ 30 | per-structure | 15/15 PASS |
| G1-f neck p10 ≥ 4 vox | per-structure | 15/15 PASS |
| G1-g Ni size ordering | class | **PASS** (463 → 536 → 635 nm) |
| **G1-h YSZ length-scale ordering** | class | **FAIL** (tie, see §3) |
| **G1-i YSZ pristine ordering** | class | **FAIL** (inverted, see §4) |

**Per the Step 1 instruction, σ, jitter and every other placement parameter have
been left exactly as specified. Nothing was tuned, re-solved or re-seeded in
response to a gate outcome.**

The G1-i failure is the important one, it was **predicted in writing before the
run** (`cmlib/project2.py` docstring, committed with the code), and it is the
"vital finding" the instruction anticipated: **the smoothed-random-field YSZ
placement model cannot reproduce the real morphological trend, and no choice of
σ can fix it.**

---

## 1. A2.3 σ-scaling, as implemented

Rule frozen before the run: σ ∝ the analog's own target particle diameter,
anchored at medium = 3.0.

| analog | D_particle | σ (vox) |
|---|---|---|
| fine | 420 nm | 2.603 |
| medium | 484 nm | 3.000 |
| coarse | 560 nm | 3.471 |

## 2. What passed, and passed well

Per-structure reproducibility is excellent. Class means over 5 seeds:

| | fine | medium | coarse |
|---|---|---|---|
| Φ_Ni (real target) | 0.3205 (0.322) | 0.2521 (0.250) | 0.2300 (0.229) |
| Φ_YSZ (real target) | 0.4238 (0.421) | 0.3869 (0.388) | 0.3835 (0.384) |
| Ni P_span | 1.0000 | 1.0000 | 1.0000 |
| Ni clusters | 1 | 1 | 1 |
| **SNOW node d_volwt** | **463 nm** | **536 nm** | **635 nm** |
| YSZ P_span | 0.99953 | 0.99887 | 0.99862 |

**G1-g passes cleanly** — measured SNOW node size orders fine < medium < coarse
with seed-to-seed spread of ≤ 5 nm, so the Ni side of the platform is sound.
Every analog percolates in Ni as a single cluster, and Φ_YSZ is within 0.97 % of
the real anode value on all 15 structures.

---

## 3. G1-a and G1-h: two failures that are *not* mechanism failures

### 3.1 G1-a — one structure of fifteen, from freezing `neck_scale`

| analog | seed | Φ_Ni | deviation |
|---|---|---|---|
| coarse | 4 | 0.23651 | **+3.28 %** |

Every other structure is inside ±1.5 %. Cause is plain: `neck_scale` was frozen
per analog at the **seed-0 solved value** (design memo §1.2), so seeds whose
neck-width draw runs fat carry that error straight into Φ_Ni. Per-analog spread
of the deviation:

| analog | min | max | sd |
|---|---|---|---|
| fine | −0.98 % | −0.20 % | 0.30 |
| medium | −0.02 % | +1.42 % | 0.59 |
| coarse | −0.12 % | **+3.28 %** | 1.33 |

Coarse is most exposed because it has the fewest necks (224 pairs vs 732 for
fine), so its width draw has the largest sampling variance.

**Fix, available but NOT applied:** solve `neck_scale` per structure seed rather
than per analog, exactly as Φ_Ni is the quantity we actually want matched. This
is a ~9-iteration bisection per structure, a few seconds each. **I have not done
it** — it changes the frozen generation protocol and needs your approval, and it
would be improper to apply it while other gates are failing.

### 3.2 G1-h — the statistic is quantized, the direction is right

| | fine | medium | coarse |
|---|---|---|---|
| **YSZ EDT p50 (gate statistic)** | **40.00** | **40.00** | **44.72 nm** |
| YSZ EDT mean (supplementary) | 42.41 | 43.72 | 48.65 nm |

The gate requires *strictly* increasing; fine and medium **tie at exactly
40.0 nm**, so it fails. But 40.0 nm is exactly 2 voxels: on a 20 nm grid the
Euclidean distance transform takes discrete values (20, 28.28, 40, 44.72, …) and
a median over all phase voxels snaps to one of them. The **mean** EDT, which is
not snapped, orders 42.41 < 43.72 < 48.65 — **the length-scale trend is
present and correctly ordered; the p50 statistic cannot resolve it.**

**This is the same failure mode as Project 1's D2 diagnostic** — a percentile of
an EDT taken over all voxels of a phase, dominated by voxel quantization and the
surface shell, returning identical values for structures known to differ. I
flagged the risk before running and am reporting it as a failure rather than
quietly switching to the mean, which would be choosing the statistic after
seeing the result.

**Proposed re-specification, for your approval, NOT applied:** score G1-h on
mean YSZ EDT (or a sub-voxel-resolved length scale such as the two-point
correlation length). On the evidence above this would pass. **This is a change
to a pre-registered criterion and is yours to authorize.**

---

## 4. G1-i — the vital finding: the YSZ placement model has the trend backwards

### 4.1 The measurement

| | fine | medium | coarse | required |
|---|---|---|---|---|
| YSZ P_span | 0.99953 | 0.99887 | 0.99862 | decreasing — **PASS** |
| **YSZ n_clusters (raw)** | **211** | **220** | **142** | increasing — **FAIL** |
| **YSZ n_clusters per Mvoxel** | **38.6** | **33.7** | **22.2** | increasing — **FAIL** |

The first clause passes. **The second fails, and not marginally — it is
monotonically inverted.** Real pristine YSZ fragments *more* as the anode
coarsens; the synthetic YSZ fragments *less*.

### 4.2 Against the real data, side by side

| analog | synthetic clusters/Mvox | **real** clusters/Mvox | over-fragmentation | synthetic P_span | **real** P_span |
|---|---|---|---|---|---|
| fine | 38.6 | **1.03** | **37.5×** | 0.9995 | 0.9989 |
| medium | 33.7 | **2.46** | **13.7×** | 0.9989 | 0.9880 |
| coarse | 22.2 | **6.07** | **3.7×** | 0.9986 | **0.9246** |

Two distinct defects, both visible here:

1. **Inverted trend** (the gate failure): real 1.03 → 2.46 → 6.07 rising;
   synthetic 38.6 → 33.7 → 22.2 falling.
2. **Uniform over-fragmentation** (not gated, but disqualifying on its own):
   the synthetic YSZ is **4–38× more fragmented per unit volume than real at
   every coarseness**, and the discrepancy is worst exactly where the model is
   supposed to be most benign. The synthetic YSZ P_span is also pinned at
   ~0.999 for all three classes while the real range is 0.925–0.999 — the model
   has almost no dynamic range in the variable O3 is meant to move.

### 4.3 Why σ cannot fix this — a structural argument, not a tuning problem

`add_ysz_pore` thresholds a Gaussian-smoothed Gaussian random field at the
percentile that hits the target volume fraction. At **fixed Φ**, changing σ
**rescales the morphology without changing its topology**: the excursion set of
a smoothed random field is statistically self-similar under σ, so the number of
connected components per unit volume falls roughly as σ⁻³ while the
component-size distribution simply stretches. My pre-run σ sweep measured
exactly this (medium analog, Φ_YSZ fixed at 0.388):

| σ (vox) | 1.5 | 3.0 | 5.0 | 7.0 | 10.0 |
|---|---|---|---|---|---|
| YSZ components | 615 | 246 | 98 | 47 | 38 |

**Monotonically decreasing.** So "make the coarse analog's YSZ coarser" and
"make the coarse analog's YSZ more fragmented" are *opposing* requests in this
model, and no σ satisfies both. Increasing σ for coarse — which G1-h needs —
necessarily worsens G1-i. **The two gates are in direct conflict under this
placement model.** That is the finding.

Physically, the reason is that a thresholded random field has **no notion of
particles or of contacts between them**. Real YSZ in a coarse anode is a
partially-sintered skeleton of large grains with comparatively few, small
inter-grain contacts — so coarsening reduces contact redundancy and *raises*
fragmentation. A random field coarsens by inflating one continuous blob
structure, which *lowers* fragmentation. The model is missing the very mechanism
that makes the real coarse anode's YSZ backbone fragile — and that mechanism is
precisely the thing O3 was designed to attack.

### 4.4 What this means for O3 and for the pilot

**O3 has no lever, and this is now measured rather than suspected.** The design
memo (§2.3, §5.2) named this as risk #2 and made G1-i the test. The test
returned the bad answer. Building O3a on this YSZ would mean applying a
backbone-failure operator to a phase whose backbone is (a) far more robust than
real in absolute terms, (b) *most* robust in exactly the class where reality is
*least* robust, and (c) essentially identical across the three classes
(P_span 0.9995 / 0.9989 / 0.9986). Any YSZ result obtained on it would be an
artifact of the placement field, not a mechanism finding.

**The Ni half is unaffected.** G1-c, G1-e, G1-f and G1-g all pass; the Ni
analogs are sound and O1/O2 could in principle proceed on them. But C2 and C3 —
two of the three Go criteria — are YSZ criteria, so the pilot as designed cannot
reach a GO verdict without a YSZ phase worth damaging.

---

## 5. Options, none of them implemented

Presented for your decision, with a recommendation. **I have not started any of
these.**

**Option A — replace the YSZ placement model with a particle-based generator.**
Place YSZ as a second population of overlapping spheres with explicit contacts,
the same construction already used and validated for Ni, with contact size and
coordination scaling with analog coarseness. This directly supplies the missing
mechanism (contacts), gives O3 a designed throat population instead of a
random-field artifact (fixing the §2.3 sub-resolution-throat objection at the
same time), and is the only option that can produce *both* G1-h and G1-i.
Cost: a real generator-development effort — the Ni generator took Platform v2's
whole qualification cycle — plus re-qualification. It also changes the Ni/YSZ
interface, so TPB and Φ would need re-checking.

**Option B — proceed Ni-only.** Run O1 and O2 against C1 alone, report the YSZ
half as blocked by the placement model. Cheap and honest, but it abandons the
divergence (C3), which is the scientific point of Project 2 and the one part of
the signature Step 0 showed to be strong and real.

**Option C — re-specify the YSZ outcome away from fragmentation.** Score O3 on
YSZ P_span alone and drop the cluster-count clause. **I recommend against this.**
The synthetic P_span range across classes is 0.9995–0.9986 — three parts in a
thousand — against a real range of 0.925–0.999. There is no dynamic range to
detect anything, so this would not rescue the pilot; it would only hide the
failure.

**Option D — stop Project 2 here** and report the placement-model insufficiency
as the result, together with the Step 0 measurement (which stands on its own as
a first-ever measurement of YSZ percolation on this dataset, with a strong
divergent signature).

**My recommendation: Option A, with Option B as the fallback if the generator
effort is not worth it.** Reasoning: Step 0 established that the divergence is
real and large (YSZ retention 0.958 / 0.865 / 0.000), which makes it worth a
proper YSZ generator; and the specific defect is now precisely characterized, so
Option A has a concrete target to hit rather than an open-ended remit — it must
reproduce rising fragmentation with coarseness and a P_span range of order 0.9,
not 0.999. Option C should be refused outright.

---

## 6. Recorded limitations of Step 1

1. **Class gates are evaluated on 5-seed means with no significance test**, per
   the project-wide no-p-values rule. G1-i's inversion is far larger than the
   seed spread (coarse 111–166 clusters vs fine 176–242, non-overlapping), so
   the conclusion does not rest on the mean alone; G1-h's tie is exact, not
   marginal.
2. **`n_clusters` counts include single-voxel speckle.** In both synthetic and
   real data P_largest ≈ 0.99, so the component counts are dominated by tiny
   fragments. The comparison is like-for-like (same metric, same connectivity)
   but it is a speckle-sensitive statistic, and the real stacks carry
   segmentation speckle of their own. This weakens the *absolute* 4–38×
   over-fragmentation claim more than it weakens the *inverted-trend* claim,
   which is what the gate turns on.
3. **Domain sizes differ between analogs** (5.48 / 6.53 / 6.40 Mvox), so raw
   cluster counts are not directly comparable; the per-Mvoxel column is the
   fair one, and both fail identically.
4. **σ was set from the analog's own particle diameter**, per instruction. Using
   the *real* anode diameters (1148/1445/1715 nm) would give σ = 2.38 / 3.00 /
   3.56 — a slightly wider spread, but the σ sweep shows this makes G1-i
   *worse*, not better. The failure is not sensitive to that choice.
5. **One YSZ placement seed per structure.** Placement seeds are recorded
   (`ysz_seed` = 5000 + 100·class + struct_seed) but not replicated
   independently of structure seed.

---

## 7. Status

- Step 1: **complete, FAILED at G1-h and G1-i**, plus G1-a at 14/15.
- **STOP condition honoured.** No placement parameter was tuned; no gate
  criterion was redefined after seeing results.
- **Not run:** O1, O2, O3 — unimplemented. No damage operator exists in the
  codebase. `cmlib/damage.py`, `cmlib/synth.py` unmodified.
- Step 2 is **not** authorized to begin and should not begin until the YSZ
  placement question is decided.

**Awaiting your decision between Options A–D.**
