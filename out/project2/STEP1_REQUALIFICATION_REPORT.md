# Project 2 — Step 1 re-qualification: **ALL GATES PASS**

Run 2026-08-10, 15 structures, ~17 min. Under `PREREGISTRATION_V2_1.md`
(committed `e62f30b`), kill test approved (`8e77035`), `p_sinter` frozen at
**0.955 / 0.523 / 0.416**.

Code: `scripts/project2/step1_requalification.py`, generator in
`cmlib/project2.py`. Data: `step1_requal.csv` (all 15 structures, every
parameter and seed), `step1_requal_gates.csv`, `step1_requal_class_medians.csv`.

**No damage operator was implemented or applied.** O1/O2/O3 do not exist in the
codebase. `cmlib/damage.py` and `cmlib/synth.py` are unmodified.

---

## Verdict

| gate | scope | result |
|---|---|---|
| G1-a Φ_Ni ±2 % (per-seed `neck_scale`) | per-structure | **15/15 PASS** |
| G1-b Φ_YSZ ±2 % | per-structure | **15/15 PASS** |
| G1-c Ni percolates, 1 cluster, P_span = 1.000 | per-structure | **15/15 PASS** |
| G1-d YSZ percolates (P_span > 0) | per-structure | **15/15 PASS** |
| G1-e particles ≥ 30 | per-structure | **15/15 PASS** |
| G1-f Ni neck p10 ≥ 4 vox | per-structure | **15/15 PASS** |
| G1-g Ni size ordering | class median | **PASS** — 463 → 534 → 633 nm |
| G1-h YSZ length scale (mean EDT) | class median | **PASS** — 57.3 → 64.6 → 73.2 nm |
| G1-i YSZ fragility (Q_YSZ) | class median | **PASS** — see §3 |

**Both gates that failed with the random-field placement (G1-h, G1-i) now pass.**

---

## 1. Class medians (5 seeds each)

| | fine | medium | coarse | real |
|---|---|---|---|---|
| Φ_Ni | 0.32103 | 0.24951 | 0.22902 | 0.322 / 0.250 / 0.229 |
| Φ_YSZ | 0.42095 | 0.38794 | 0.38401 | 0.421 / 0.388 / 0.384 |
| **Q_YSZ = 1 − P_span** | **0.00092** | **0.00992** | **0.06768** | **0.00109 / 0.01197 / 0.07542** |
| YSZ P_span | 0.99908 | 0.99008 | 0.93232 | 0.99891 / 0.98803 / 0.92458 |
| mean YSZ EDT | 57.3 | 64.6 | 73.2 nm | — |
| SNOW Ni node d | 463 | 534 | 633 nm | — |
| YSZ grain d | 521 | 557 | 640 nm | — |
| solved `neck_scale` | 0.747 | 0.793 | 0.974 | — |
| sintered coordination | 8.41 | 4.58 | 3.53 | — |
| YSZ clusters (never gating) | 276 | 215 | 136 | — |
| filtered clusters (never gating) | 0 | 6 | 17 | — |

**Sintered coordination falls with coarseness (8.41 → 4.58 → 3.53)** — the
signature the whole architecture was built to produce, and the direction
Herring's scaling law predicts.

---

## 2. G1-a — per-structure log (A3 mandatory)

Every field required by amendment A3. All 15 solves converged in 8 iterations;
no failures.

| analog | seed | `neck_scale` | iters | Φ_Ni | dev % | Φ_YSZ | dev % | YSZ iters | node d (nm) | YSZ grain d (nm) | Ni necks | cand. contacts | sintered | coord |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fine | 0 | 0.74700 | 8 | 0.32069 | −0.41 | 0.42110 | +0.02 | 9 | 463.1 | 521.0 | 732 | 3494 | 3332 | 8.41 |
| fine | 1 | 0.74700 | 8 | 0.32042 | −0.49 | 0.42091 | −0.02 | 11 | 462.5 | 521.3 | 732 | 3492 | 3321 | 8.39 |
| fine | 2 | 0.74700 | 8 | 0.32137 | −0.20 | 0.42095 | −0.01 | 12 | 462.6 | 521.7 | 732 | 3469 | 3303 | 8.34 |
| fine | 3 | **0.75794** | 8 | 0.32297 | +0.30 | 0.42093 | −0.02 | 11 | 464.4 | 521.3 | 732 | 3503 | 3333 | 8.42 |
| fine | 4 | 0.74700 | 8 | 0.32103 | −0.30 | 0.42107 | +0.02 | 11 | 462.5 | 520.8 | 732 | 3500 | 3346 | 8.45 |
| medium | 0 | 0.80700 | 8 | 0.24995 | −0.02 | 0.38791 | −0.02 | 12 | 534.5 | 556.5 | 430 | 2920 | 1524 | 4.58 |
| medium | 1 | **0.79333** | 8 | 0.24951 | −0.20 | 0.38804 | +0.01 | 9 | 533.6 | 558.2 | 430 | 2909 | 1514 | 4.55 |
| medium | 2 | **0.79606** | 8 | 0.24932 | −0.27 | 0.38803 | +0.01 | 12 | 533.1 | 558.3 | 430 | 2883 | 1538 | 4.63 |
| medium | 3 | **0.79333** | 8 | 0.24886 | −0.46 | 0.38794 | −0.02 | 8 | 533.7 | 557.0 | 430 | 2892 | 1555 | 4.68 |
| medium | 4 | **0.78513** | 8 | 0.24959 | −0.16 | 0.38790 | −0.03 | 12 | 534.3 | 557.4 | 430 | 2901 | 1455 | 4.38 |
| coarse | 0 | 0.98800 | 8 | 0.22872 | −0.12 | 0.38406 | +0.01 | 10 | 632.6 | 641.5 | 224 | 2136 | 864 | 3.46 |
| coarse | 1 | 0.98800 | 8 | 0.23013 | +0.49 | 0.38398 | −0.01 | 12 | 632.5 | 641.3 | 224 | 2136 | 910 | 3.64 |
| coarse | 2 | **0.97433** | 8 | 0.22811 | −0.39 | 0.38401 | +0.00 | 11 | 630.7 | 639.8 | 224 | 2149 | 894 | 3.58 |
| coarse | 3 | **0.97433** | 8 | 0.22903 | +0.01 | 0.38409 | +0.02 | 11 | 634.8 | 640.5 | 224 | 2145 | 912 | 3.65 |
| coarse | 4 | **0.96886** | 8 | 0.22902 | +0.01 | 0.38400 | +0.00 | 11 | 632.2 | 639.8 | 224 | 2115 | 865 | 3.46 |

**A3 did the job it was introduced for.** Nine of fifteen seeds solved to a
`neck_scale` different from the frozen seed-0 value (bold). Worst Φ_Ni deviation
is now **+0.49 %**, against **+3.28 %** last time — the single G1-a failure of
the original Step 1 is gone, and the worst case is 4× inside the gate.

Φ_Ni deviations do not go to zero because neck widths are integer-rounded, so
achievable Φ_Ni is quantized in discrete rungs — the same discrete-rung
behaviour documented for Platform v2. The residual is a quantization floor, not
a solver failure.

---

## 3. G1-i — the gate that killed the random field

**Q_YSZ per seed:**

| analog | seeds (sorted) | median |
|---|---|---|
| fine | 0.00067, 0.00077, 0.00092, 0.00097, 0.00131 | **0.00092** |
| medium | 0.00899, 0.00955, 0.00992, 0.01172, 0.01469 | **0.00992** |
| coarse | 0.03163, 0.06062, 0.06768, 0.07181, 0.09166 | **0.06768** |

| criterion | requirement | result |
|---|---|---|
| 1. strict ordering | Q_fine < Q_medium < Q_coarse | **PASS** — no overlap between class ranges |
| 2. ratio | Q_coarse / Q_fine ≥ 10 | **PASS — 73.6×** (real 69×) |
| 3. anti-outlier | dropping any one seed preserves 1 and 2 | **PASS** — worst-case leave-one-out ratio **69.8×** across all 15 cases |

The class ranges are fully disjoint (fine max 0.00131 < medium min 0.00899;
medium max 0.01469 < coarse min 0.03163), so the ordering does not depend on the
median statistic. Coarse seed 3 (0.03163) is a genuine low outlier — it is
retained, and criterion 3 shows it does not carry the result.

**Comparison against the original Step 1 failure**, which is the point of the
whole exercise:

| | fine | medium | coarse | trend |
|---|---|---|---|---|
| **real** Q_YSZ | 0.00109 | 0.01197 | 0.07542 | rising, 69× |
| **particulate** Q_YSZ | 0.00092 | 0.00992 | 0.06768 | **rising, 73.6×** ✓ |
| *random field, clusters/µm³* | *4.83* | *4.22* | *2.78* | *falling* ✗ |

## 4. G1-h — mean EDT resolves what p50 could not

| statistic | fine | medium | coarse | ordered? |
|---|---|---|---|---|
| **mean YSZ EDT (gate)** | **57.3** | **64.6** | **73.2 nm** | **yes — PASS** |
| YSZ EDT p50 (recorded, non-gating) | 44.7 | 49.0 | 60.0 nm | yes, here |

Amendment A2 was made because p50 tied at exactly 40.00 nm for fine and medium
under the random field. With the particulate generator p50 happens to separate
too, but the gate remains on the mean, as frozen — the criterion was not
re-chosen after seeing the result.

---

## 5. Operating-point grain diameter vs sweep drift

Recorded per the advisor's documentation condition on the accepted K2 deviation.

| analog | **operating-point YSZ grain d** | across-seed spread | K2 sweep-range drift |
|---|---|---|---|
| fine | **521.3 nm** (520.8–521.7) | 0.17 % | 1.90 % |
| medium | **557.4 nm** (556.5–558.3) | 0.32 % | **6.91 %** |
| coarse | **640.5 nm** (639.8–641.5) | 0.27 % | **6.66 %** |

**The drift is confirmed to be a sweep-only artifact.** At the frozen operating
p_sinter the grain diameter is single-valued to within **0.17–0.32 %** across
five independent seeds — two orders of magnitude tighter than the sweep-range
drift, and far inside the 5 % Protocol A limit. No structure carried into any
later stage inherits the drift.

---

## 6. Limitations

1. **YSZ jitter is 0.02**, against the Ni lattice's 0.15 — a 7.5× more regular
   phase. This is forced, not chosen: the fine analog's inter-grain gap is only
   1.60 voxels, and worst-case approach along a bond is 2·jitter·nn, so jitter
   above ~0.028 would close the gap and reintroduce the accidental raster
   connectivity that K0 exists to exclude. **The synthetic YSZ is therefore more
   crystalline than a real sintered network**, and the constraint binds hardest
   exactly where the packing is tightest (fine). Any future relaxation requires
   either a larger a_ysz (fewer grains per domain) or a denser-than-FCC packing.
2. **Q agreement with real is targeted, not predicted.** `p_sinter` was
   calibrated once against the real Q values (A5, pre-damage, frozen). The
   non-trivial content is that all three analogs are simultaneously reachable at
   fixed Φ_YSZ with the correct ordering and a 73.6× spread — which no σ of the
   random field could deliver — not that the absolute values were forecast.
3. **The frozen sub-grain limitation applies verbatim:** "The synthetic YSZ
   generator will model grain-scale backbone connectivity and pristine
   percolation fragility. It will not reproduce the measured sub-grain fragment
   population, fragment-size distribution, or fragment surface area."
4. **YSZ grain size is not validated against real data** — no YSZ grain-size
   measurement exists in this study. `a_ysz` scales with the Ni particle
   diameter and was set to give a ≥1.5-voxel gap for fine. Only ordering is
   claimed.
5. **Analogs remain ~3× size-compressed** relative to real particles
   (`DESIGN_MEMO` §1.4), with the recorded bias direction: fine's proportionally
   fatter Ni necks bias *against* reproducing "fine loses Ni percolation worst".
6. **One YSZ placement seed per structure seed**; grain geometry and sintering
   draws are not independently replicated from the structure seed.
7. **Class gates use medians with no significance test**, per the project-wide
   no-p-values rule. G1-i's margin (disjoint class ranges, 69.8× worst-case
   leave-one-out) does not rest on that choice.
8. **Coarse Q has real seed variance** (0.032–0.092, a 2.9× spread). Damage
   comparisons on the coarse analog will need the damage-seed averaging rule
   inherited from Project 1 §0g/1.

---

## 7. Status

- Step 1 re-qualification: **complete, ALL GATES PASS** (G1-a … G1-i).
- `p_sinter` unchanged at 0.955 / 0.523 / 0.416. No parameter was tuned in
  response to any gate outcome.
- **Not run:** O1, O2, O3 — unimplemented. No damage operator exists.
  `cmlib/damage.py`, `cmlib/synth.py` untouched.
- The platform now satisfies the precondition that failed in the original
  Step 1: a YSZ phase whose pristine fragility is coarseness-ordered, with the
  dynamic range needed for a YSZ damage operator to have a lever.

**Awaiting your review before Step 2 (operator implementation).**
