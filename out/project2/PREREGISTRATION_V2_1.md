# Project 2 — Pre-registration v2.1 (FROZEN)

Approved 2026-08-10. Supersedes `DESIGN_MEMO.md` §4.2, `STEP1_REPORT.md` §5, and
`AUDIT_RESPONSE_AND_PREREG_V2.md` §6–7 wherever they conflict.

**Committed before any YSZ generator code is written.** Nothing in this document
may be changed after a damage operator is run.

---

## 0. Changes from v2 to v2.1

| # | v2 | **v2.1** |
|---|---|---|
| 1 | Volume control left as a choice between Protocol A and B | **Protocol A frozen** — fixed Φ_YSZ, grain-size compensation ≤ 5 % (§4.2) |
| 2 | K2 failure consequence unstated | **K2 failure ⇒ Option A not viable with this architecture.** Report achieved Q range, fall back to Option B. **No rescue** via grain geometry, domain size, or contact model (§5.2) |
| 3 | Seed counts not distinguished | **Kill test 3 seeds/condition; full Step 1 re-qualification 5 seeds/class** (§6) |
| 4 | G1-i: P_span primary, filtered count secondary-but-gating-ish | Filtered counts **explicitly non-gating**; Q_YSZ = 1 − P_span sole primary; median-based, anti-outlier clause added (§2.1) |
| 5 | P_span definition assumed from Project 1 | **Fully frozen and written out** — component rule, connectivity, boundaries, axis, domain mapping (§1) |
| 6 | EDT estimator loosely specified | **Estimator frozen** — mask, isotropy, units, domain (§2.2) |
| 7 | A5 stated as "never re-derive" | **Clarified**: one pre-damage recomputation from the measured generator calibration curve is permitted and must be logged (§4.3) |
| 8 | Kill test K1–K5 | **K0 topology sanity check added and made mandatory**; K3 given a one-shot recalibration path (§5) |
| 9 | — | **New explicit limitation** on sub-grain fragments (§2.1) |

---

## 1. Frozen definition of P_span and Q_YSZ

**Q_YSZ ≡ 1 − P_span**, the isolated (non-spanning) volume fraction of the YSZ
phase.

**P_span is the SPANNING-CLUSTER volume fraction, not the largest-component
fraction.** Precisely, as implemented in `cmlib.percolation.percolation_summary`
(frozen in Project 1, unmodified):

| element | frozen choice |
|---|---|
| **component rule** | voxels of the phase in connected component(s) that contain at least one voxel in the **first** slice *and* at least one voxel in the **last** slice along the transport axis |
| **denominator** | total voxels of that phase (not of the domain) |
| **connectivity** | **6-neighbour (face-sharing)** for phase connectivity. Not 18, not 26, not graph-edge. This is the rule for which the Phase-0 validation gate (p_c = 0.3116077) is defined and the rule used for every percolation number in Projects 1 and 2 to date. |
| **boundary conditions** | **free (non-periodic) on all six faces.** No wrap-around, no artificial boundary pores. |
| **transport axis — synthetic** | **axis 0 (z)**, the generator's designed percolation axis (z-boundary particle layers sit exactly on the domain faces) |
| **transport axis — real data** | **axis 2 (x)**, through-thickness, as used in Phase 5 and Step 0 |
| **multiple spanning clusters** | unioned, then divided by phase voxels |

**Note, recorded once so it is never re-litigated:** on all three *pristine* real
stacks P_span and P_largest are numerically identical (0.99891 / 0.98803 /
0.92458), so this choice is empirically inconsequential for G1-i, which is a
pristine-only gate. It matters only in the degraded state (`coarse_post`:
P_span = 0.0000 vs P_largest = 0.1235), which G1-i does not use.

**Domain size / compressed-scale mapping.** Synthetic domains are the
`DESIGN_MEMO` §1.2 per-analog shapes at a uniform **20 nm** voxel:
fine 169×180×180, medium 181×190×190, coarse 181×188×188. Real stacks are at
19.53×19.53×20.0 / 24.41×24.41×25.0 / 29.14×29.14×30.0 nm. Synthetic particle
diameters are **~3× compressed** relative to real (`DESIGN_MEMO` §1.4).
**Consequence, frozen:** Q_YSZ is dimensionless and resolution-independent and
is therefore compared **directly** across synthetic and real. Any count- or
density-based quantity is **not** comparable across the two and is never gated.

---

## 2. Amended gates

### 2.1 G1-i (v2.1) — YSZ pristine fragility

**Raw cluster counts are dropped as a gate metric in every form** (raw,
per-Mvoxel, per-µm³). **Filtered cluster counts are dropped as a primary gate**
and are diagnostic only.

> **G1-i pass criteria**, on **class medians** over 5 seeds (full qualification):
> 1. **Q_fine < Q_medium < Q_coarse** (strict);
> 2. **Q_coarse / Q_fine ≥ 10**;
> 3. **the ordering is not driven by a single outlier seed** — operationalised
>    as: removing any one seed from any class leaves criteria 1 and 2 intact.

Real-data reference values (`step0_ysz_percolation.csv`): Q = 0.00109 / 0.01197 /
0.07542, ratio **69×**.

> **Explicit limitation, frozen (quoted verbatim into every downstream
> artifact):** "The synthetic YSZ generator will model grain-scale backbone
> connectivity and pristine percolation fragility. It will not reproduce the
> measured sub-grain fragment population, fragment-size distribution, or
> fragment surface area."

Justification for that limitation is measured, not assumed: real isolated YSZ
consists of ~100 nm sub-grain fragments (median isolated cluster 68 / 125 / 78
voxels) that exceed the whole-grain volume ceiling by 22–232×, an object class a
whole-grain generator cannot represent
(`AUDIT_RESPONSE_AND_PREREG_V2.md` §2.2–2.3).

### 2.2 G1-h (v2.1) — YSZ length scale

> **Estimator, frozen:** Euclidean distance transform over **all YSZ voxels** of
> the same mask and domain used for P_span; **isotropic** sampling at
> **20.0 nm/voxel** (`scipy.ndimage.distance_transform_edt`, scalar sampling);
> statistic = **arithmetic mean over YSZ voxels**; units **nm**.
>
> **Pass criterion:** strict ordering of **class medians**,
> mean-EDT(fine) < mean-EDT(medium) < mean-EDT(coarse).

EDT p50 is recorded but **not gating** (it failed Step 1 on an exact
quantization tie at 40.00 nm, not on direction). Specific surface area or any
other continuous length-scale metric may be reported as **supplementary only**
and is **not** a gate, since it is not frozen here.

### 2.3 G1-a (v2.1) — composition control

`neck_scale` is solved **per structure seed** to hit the analog's Φ_Ni target.
Gate remains **±2 %**. This is a numerical-precision fix, not a physics change.

> **Per-structure log, mandatory:** seed; solved `neck_scale`; achieved Φ_Ni;
> achieved Φ_YSZ; achieved grain/node size; number of necks (Ni) and candidate
> and sintered contacts (YSZ); solve iterations; failure reason if unconverged.

### 2.4 Unchanged gates

G1-b (Φ_YSZ ±2 %), G1-c (Ni percolates, single cluster, P_span = 1.000),
G1-d (YSZ percolates, P_span > 0, value recorded, no threshold),
G1-e (particles ≥ 30), G1-f (Ni neck p10 ≥ 4 vox), G1-g (Ni SNOW node size
orders fine < medium < coarse) stand as written.

---

## 3. A4 — particulate YSZ generator requirements

Implemented only after this document is committed. The generator must satisfy,
each verified rather than asserted:

1. **YSZ grains are explicit geometric objects** with recorded centres and size.
2. **Contacts are explicit graph edges** in a recorded candidate-contact list.
3. **Sintered contacts are topologically connected** (6-connectivity).
4. **Unsintered contacts are topologically disconnected** (6-connectivity).
5. **Changing `p_sinter` changes topology without changing grain volume beyond
   the Protocol-A tolerance.**
6. **No reliance on accidental raster overlap for connectivity** — every
   connection is an explicitly rasterized neck. Verified by K0 at `p_sinter = 0`.

`cmlib/damage.py::add_ysz_pore` is **not modified and not deleted**; Project 1
and Step 0/1 results depend on it. The new generator is additive.

---

## 4. Frozen protocols

### 4.1 Volume control — **PROTOCOL A (frozen)**

- **Φ_YSZ is a hard target** (real per-anode: 0.421 / 0.388 / 0.384).
- **Grain core radius may be adjusted** to compensate for neck-volume changes as
  `p_sinter` varies.
- **Mean grain diameter drift is reported** for every structure.
- **Drift must be ≤ 5 %** during the kill test, unless the Q change is
  demonstrably topological (i.e. explained by contact connectivity, not volume).

Protocol B (fixed diameter, drifting Φ) is **not** used.

### 4.2 Initial p_sinter values and their derivation source

| analog | target Q_YSZ (real) | **initial p_sinter** |
|---|---|---|
| fine | 0.00109 | **0.70** |
| medium | 0.01197 | **0.56** |
| coarse | 0.07542 | **0.42** |

**Derivation source:** bond percolation on the simple-cubic lattice
(p_c = 0.2488), computed directly at L = 24 with 3 replicates
(`AUDIT_RESPONSE_AND_PREREG_V2.md` §5), reading `p_sinter` off the measured
1 − P_span curve at the three real targets. Sign check: the required ordering
p_sinter(fine) > p_sinter(medium) > p_sinter(coarse) is decreasing with
coarseness, as Herring's sintering scaling law predicts (time to reach a given
x/R scales as Rⁿ, n = 3–4, so at fixed schedule larger grains reach smaller
relative neck size).

**These are initial predictions, not immutable constants.**

### 4.3 A5 — p_sinter anti-tuning rule (clarified)

- `p_sinter` **may not** be adjusted to improve a damage result.
- `p_sinter` **may** be derived/calibrated **before damage** using the **actual
  generator graph**.
- If the actual generator graph is not effectively SC bond percolation,
  `p_sinter` is **recomputed once** from the measured generator calibration
  curve, before damage.
- **Any recomputation is logged and frozen before damage operators run.**
- **No** `p_sinter` change after observing O1/O2/O3 outcomes.
- **No** per-damage-seed or per-damage-metric `p_sinter` tuning.

---

## 5. Kill-test specification

Not run until this document is committed. **3 seeds per condition.**

### K0 — topology sanity check (MANDATORY, gates everything else)

At fixed nominal grain geometry, obeying Protocol A:

| condition | expectation |
|---|---|
| `p_sinter = 0` | **no system-spanning YSZ network**; P_span ≈ 0 (or at most the largest isolated-grain fraction). **If it still percolates, the generator is invalid.** |
| `p_sinter = 1` | highly connected; **P_span ≥ 0.98** unless finite-size effects explain otherwise |

**If K0 fails, stop and report. Do not proceed to K1–K5.**

### K1 — volume fraction control

For every kill-test structure:
`|Φ_YSZ_achieved − Φ_YSZ_target| / Φ_YSZ_target ≤ 0.02`.

Report: achieved Φ_YSZ; grain diameter / equivalent grain size; number of YSZ
grains; number of candidate contacts; number of sintered contacts; realized mean
coordination among sintered contacts.

### K2 — decoupling / lever test — **TRUE KILL CRITERION**

At fixed nominal grain-diameter distribution and fixed volume-control protocol,
vary `p_sinter` and show Q_YSZ moves by **at least one order of magnitude**.

**Pass, on medians over 3 seeds:**
- `Q_low-p / Q_high-p ≥ 10`
- `Q_low-p ≥ 0.05`
- `Q_high-p ≤ 0.005`

equivalently low-connectivity `P_span ≤ 0.95` and high-connectivity
`P_span ≥ 0.995`.

Volume constraints during K2: Φ_YSZ drift obeys Protocol A; if grain size is
adjusted to preserve Φ_YSZ, **mean grain diameter drift is reported** and must
be ≤ 5 % unless the Q change is clearly topological.

> **K2 failure consequence (frozen):** if `Q_low/Q_high < 10`, **Option A is not
> viable with this generator architecture.** Report the achieved Q range and
> **fall back to Option B (Ni-only).** **Do not** attempt to rescue K2 by
> changing grain geometry, domain size, or contact model.

### K3 — coarse target check (3 seeds)

At frozen coarse `p_sinter`: median P_span ∈ **[0.90, 0.95]**, i.e.
Q_YSZ ∈ [0.05, 0.10].

Pass rule: median inside bracket; **≥ 2 of 3 seeds** inside bracket; **no seed
more than 0.03 absolute** beyond the bracket.

If **K2 passes but K3 fails only because the SC-derived mapping is wrong**:
recompute coarse `p_sinter` **once** from the measured generator calibration
curve, freeze the new value, rerun K3. **If K3 fails again, stop and report;**
fallback to Option B may be required.

### K4 — fine target check (3 seeds)

At frozen fine `p_sinter`: median **P_span ≥ 0.995**, with **≥ 2 of 3 seeds
≥ 0.995**.

If the generator only reaches ~0.990, this does **not** automatically pass: K2
must still show `Q_low/Q_high ≥ 10`. If it does not, **the kill test fails.**

### K5 — filtered cluster diagnostics — **SECONDARY ONLY, NEVER GATING**

Report: filtered cluster count; filtered cluster volume fraction; the minimum
cluster-size threshold used; threshold sensitivity if cheap.

---

## 6. Seed counts

| activity | seeds |
|---|---|
| **kill test** (K0–K5) | **3 per condition** — generator validation |
| **full Step 1 re-qualification** | **5 per analog class** |

---

## 7. Final rule

**No damage modeling, no O1/O2/O3, and no full qualification until the amended
kill test passes and the advisor reviews the result.**

---

## 8. Deviations from the advisor's instructions

**None.** All of A1–A5, the three v2.1 additions, and the K0–K5 specification
are adopted as written. Two items were determinate rather than ambiguous and are
recorded here as choices made under the instruction to freeze them, not as
deviations:

1. **P_span component rule** — the instruction offered largest-component vs
   spanning-cluster. **Spanning-cluster** is frozen (§1), because it is the
   Project-1 definition already used for every committed percolation number and
   for the Step 0 real-data measurement. Empirically inconsequential for G1-i
   (identical on all pristine stacks).
2. **Transport axis differs between synthetic (axis 0) and real (axis 2)** — each
   is its own domain's designed/physical transport direction. This is inherited
   from Step 0 and Step 1, stated here explicitly rather than changed.
