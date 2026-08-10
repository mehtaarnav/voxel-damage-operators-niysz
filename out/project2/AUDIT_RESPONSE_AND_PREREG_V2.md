# Audit response and corrected pre-registration (v2)

Written 2026-08-10 in response to the Trap 1 / Trap 2 audit. Every claim below
is checked against measurement or derivation; where the audit is right I say so,
where it is wrong I say that too, and I found one error of my own that neither
of us had caught.

New evidence produced for this response (no generator code written):
`scripts/project2/audit_ysz_cluster_sizes.py` → `audit_ysz_cluster_sizes.csv`
(real YSZ cluster-size distributions, all three pristine stacks), plus a
bond-diluted lattice percolation calculation reported in §5.

---

## 0. Verdict in one table

| audit claim | status |
|---|---|
| **T1** — linear contact scaling ⇒ λ constant ⇒ N_iso ∝ 1/d³ ⇒ same inverted failure | **Correct in substance. Accepted.** Two corrections to the reasoning (§1). |
| **T1 fix** — break geometric similarity; coarse gets smaller relative neck x/R | **Correct, and grounded in real sintering physics (Herring scaling). Accepted.** Now quantified (§5). |
| **T2 arithmetic** — isolated volume supports ~0.088 whole particles/Mvox at real coarse | **Arithmetic method right, number wrong by 3.2×** — it assumes 20 nm voxels for all stacks; coarse_pre is 29.14 × 29.14 × 30 nm. Correct ceiling is **0.279**/Mvox (§2.1). |
| **T2 premise** — the 6.07 count is "dominated entirely by segmentation speckle" | **REFUTED by measurement.** Voxel-scale specks are 2.6–20.3 % of clusters and carry **0.003–0.09 %** of isolated volume. Median isolated cluster is 68–125 voxels (§2.2). |
| **T2 conclusion** — a particulate generator cannot reproduce real raw `n_clusters`; the metric is broken as a cross-model comparator | **UPHELD, for a stronger reason than speckle** (§2.3). |
| **T2 fix** — drop raw `n_clusters`; score on P_span divergence, or filtered count | **Accepted, with a refinement:** filtered count is still not commensurable for fine and medium; P_span (a *volume* fraction) is the defensible primary (§4). |
| **My Step 1 report** — "4–38× over-fragmented" | **My error.** I compared clusters/Mvoxel across stacks of different voxel size. Correct physical values 36 / 26 / 12× (§3). Verdict unchanged. |

---

## 1. Trap 1 — accepted, with two corrections to the reasoning

The algebra is right: at fixed Φ, N ∝ 1/d³; linear contact scaling gives
V_int ∝ d³; so λ = N·V_int is scale-invariant, P_iso depends only on λ, and
isolated-cluster *density* inherits the 1/d³ prefactor. Fine would always
out-fragment coarse. **The criticism lands squarely on what I actually wrote** —
"contact size and coordination scaling with analog coarseness" is precisely the
geometric-similarity construction that fails.

**Correction 1 — the mechanism differs for the generator I proposed to reuse.**
The argument as stated is a Boolean/random-packing argument, where coordination
*emerges* from an interaction volume. The Ni generator is a **jittered lattice
with an explicit pair list**: coordination is a design input (6 nearest
neighbours), not an emergent N·V_int. So λ is not forced to be scale-invariant
by construction. The 1/d³ *number-density* prefactor, however, is unavoidable in
either model, and that is the part that does the damage. The trap survives via
the prefactor, not via λ.

**Correction 2 — the required margin is much larger than the trap implies, so
the fix has ample room.** Quantified at both scales:

| | d_fine → d_coarse | N_fine/N_coarse | P_iso ratio needed to invert |
|---|---|---|---|
| my analogs (compressed) | 420 → 560 nm | 2.37 | > **2.37** |
| real anodes | 1148 → 1715 nm | 3.33 | > **3.33** |

The real data implies an available P_iso ratio of
(1−0.9246)/(1−0.9989) = **68×**. Against a 2.4–3.3× headwind that is a margin of
20–29×. **The sintering-yield fix does not merely work; it works comfortably.**

**Correction 3, which matters more than either — switching the metric to a
volume fraction dissolves the 1/d³ term entirely.** The prefactor afflicts
*counts per unit volume*. It does not touch the isolated *volume* fraction
(= 1 − P_span), which is set by λ alone. So the corrected metric in §4 removes
Trap 1's mathematical core. **The physics fix is still required** — with
p_sinter held constant, λ is constant, P_iso is constant, and all three analogs
land on the same P_span, failing the gate a third way. Both corrections are
needed, for different reasons.

---

## 2. Trap 2 — arithmetic corrected, premise refuted, conclusion upheld

### 2.1 The ceiling calculation, redone with per-sample voxel sizes

The audit assumes 20 nm voxels throughout ("1 Mvoxel = 8 µm³"). The three
pristine stacks were imaged at resolutions roughly proportional to their feature
size:

| stack | spacing (nm) | 1 Mvoxel = |
|---|---|---|
| fine_pre | 19.53 × 19.53 × 20.0 | 7.63 µm³ |
| medium_pre | 24.41 × 24.41 × 25.0 | 14.90 µm³ |
| coarse_pre | 29.14 × 29.14 × 30.0 | **25.47 µm³** |

Redoing the coarse ceiling with 25.47 µm³: YSZ volume = 9.78 µm³, whole particle
= 2.641 µm³, isolated 7.54 % ⇒ **0.279 whole particles/Mvoxel**, not 0.087.
Observed 6.07 exceeds it by **22×**, not 69×. Direction unchanged, magnitude 3.2×
less extreme.

Measured directly (`audit_ysz_cluster_sizes.csv`), per Mvoxel:

| | fine | medium | coarse |
|---|---|---|---|
| whole-particle ceiling | 0.0044 | 0.0438 | **0.279** |
| observed n_clusters | 1.025 | 2.458 | **6.068** |
| exceeded by | 232× | 56× | **22×** |

### 2.2 The speckle premise is refuted

If the counts were "dominated entirely by segmentation speckle", voxel-scale
objects should dominate. They do not:

| | fine | medium | coarse |
|---|---|---|---|
| **median** isolated cluster | **68 vox** | **125 vox** | **78 vox** |
| mean isolated cluster | 449 | 1,891 | 4,772 |
| clusters ≤ 8 vox | 20.3 % | **2.6 %** | 5.8 % |
| **isolated volume in ≤ 8 vox clusters** | **0.09 %** | **0.004 %** | **0.003 %** |

The median isolated YSZ cluster is 68–125 voxels — an equivalent sphere of
roughly 100–150 nm. In the medium stack only **2.6 %** of clusters are
voxel-scale. Speckle carries three parts in 100,000 of the isolated volume.
**These are genuine sub-grain fragments, not segmentation noise.**

### 2.3 The conclusion survives anyway — and the correct reason is stronger

Real isolated YSZ is ~100 nm fragments of grains whose diameter is 1148–1715 nm.
They are far too large to be speckle and far too small to be whole grains: only
**0 / 1 / 25** clusters (of 975 / 1598 / 2917) exceed one whole particle volume.

**A generator that rasterizes whole spheres cannot produce this object class at
all.** Its smallest disconnectable unit is one entire grain. So the
incommensurability is not "the real data has noise the model lacks" (fixable by
filtering) but "the real data's fragments live at a length scale the model
cannot represent" (not fixable by filtering). **Raw `n_clusters` must be dropped
as a cross-model gate metric — the audit's conclusion — but filtering alone does
not rescue it either.** See §4.

---

## 3. An error of my own, corrected

My Step 1 report compared **clusters per Mvoxel** between synthetic (20 nm
voxels) and real (19.7 / 24.6 / 29.4 nm), and between real stacks of differing
resolution. A Mvoxel is not a fixed physical volume, so those ratios were wrong.
Corrected to physical density:

| | fine | medium | coarse |
|---|---|---|---|
| real clusters/µm³ | 0.1344 | 0.1650 | 0.2382 |
| synthetic clusters/µm³ | 4.826 | 4.216 | 2.779 |
| **over-fragmentation** | **36×** | **26×** | **12×** (I reported 37.5 / 13.7 / 3.7) |

Also corrected: the **real** raw-count ordering is far weaker than I presented.
In physical density it rises only **1.77×** fine → coarse, not the 5.9× implied
by the per-Mvoxel figures. This further undermines raw `n_clusters` as a gate
target and independently supports the audit's recommendation. **The Step 1
verdict — synthetic trend inverted, G1-i failed — is unchanged.**

---

## 4. The corrected metric, derived rather than asserted

Three candidates, scored on the real data:

| metric | fine | medium | coarse | ordering | ratio | dimensionless? | reachable by particulate generator? |
|---|---|---|---|---|---|---|---|
| raw clusters/µm³ | 0.1344 | 0.1650 | 0.2382 | ✓ | 1.77× | no | **no** (22–232× above ceiling) |
| filtered > 0.1 V_p, /µm³ | 0.00069 | 0.00661 | 0.01045 | ✓ | 15.2× | no | **partly** — see below |
| **1 − P_span (isolated volume fraction)** | **0.00109** | **0.01197** | **0.07542** | ✓ | **69×** | **yes** | **yes** |

**Why not the filtered count**, despite it being a large improvement: compare it
against the whole-particle ceiling in the same units (per µm³) — fine
0.00069 vs ceiling 0.00058, medium 0.00661 vs 0.00294. **The filtered count still
exceeds what whole grains can supply, by 1.2× for fine and 2.2× for medium.**
Only coarse fits (0.01045 vs 0.01096). So filtering fixes the metric for the
coarse analog and leaves it incommensurable for the other two — exactly the
classes where the gate needs to discriminate. It is also a low-statistics
quantity for fine (**5 clusters** in the entire stack).

**1 − P_span is the right primary**: it is a volume fraction, hence dimensionless
and resolution-independent; it is insensitive to how fragments are partitioned
into clusters, which is the whole problem above; it carries the largest dynamic
range of the three (69×); it is immune to Trap 1's 1/d³ prefactor; and it is
exactly the quantity O3 must move. This is the audit's own first suggestion,
now supported by measurement rather than adopted on assertion.

---

## 5. Sintering yield — derived targets, not tuning knobs

Model: grains on the lattice, each nearest-neighbour contact forming a
conductive neck independently with probability `p_sinter`. This is bond
percolation on the simple cubic lattice (p_c = 0.2488). Direct calculation
(L = 24, 3 replicates, fraction of grains in the z-spanning cluster):

| p_sinter | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.55 | 0.60 | 0.70 | 0.80 |
|---|---|---|---|---|---|---|---|---|---|
| 1 − P_span | 0.386 | 0.180 | 0.101 | 0.050 | 0.028 | 0.014 | 0.007 | 0.0013 | 0.0002 |

Reading off the real targets:

| analog | real 1 − P_span | **required p_sinter** | margin above p_c |
|---|---|---|---|
| fine | 0.00109 | **≈ 0.70** | 2.8× |
| medium | 0.01197 | **≈ 0.56** | 2.2× |
| coarse | 0.07542 | **≈ 0.42** | **1.7×** |

Three things this establishes, all of which were assumptions before:

1. **The required P_span range is reachable**, and at values comfortably above
   the percolation threshold — the coarse analog sits 1.7× above p_c, not on it,
   so it is not in the critically-fluctuating regime where seed variance would
   swamp the signal.
2. **The required yield ordering is p_sinter(fine) > p_sinter(medium) >
   p_sinter(coarse)** — decreasing with coarseness, which is exactly what
   Herring's scaling law predicts (time to reach a given x/R scales as Rⁿ,
   n = 3–4, so at fixed sintering schedule larger grains reach smaller relative
   neck size). **The physics the audit invoked and the numbers the data demands
   agree in sign and are of sensible magnitude.** That agreement is a real
   check, not a restatement.
3. **It gives the kill test a pre-specified target** rather than a knob to turn
   until something works.

Caveats: 6-coordinated lattice (real YSZ coordination is likely higher, which
would shift all p_sinter values up without changing the ordering); grain-count
fraction is used as a proxy for the phase volume fraction, exact only for
equal-sized grains and ignoring neck volume; L = 24 with 3 replicates is a
scoping calculation, not a converged one.

---

## 6. Pre-registration v2 — formal amendments

Frozen before any YSZ generator code is written. These supersede DESIGN_MEMO
§4.2 and STEP1_REPORT §5 where they conflict.

**A1 — G1-i is re-specified.** Raw `n_clusters` is **dropped as a gate metric**
in every form (raw and per-volume), on the grounds established in §2.3 and §3:
it is not reproducible in principle by a whole-grain generator, and its real
ordering is weak (1.77×) once resolution is handled correctly.

> **G1-i (v2).** On class means over ≥3 seeds:
> **(i) primary —** pristine YSZ **P_span** orders fine > medium > coarse, and
> the isolated volume fraction (1 − P_span) spans at least **one order of
> magnitude** fine → coarse (real: 69×).
> **(ii) secondary, reported not gating —** filtered cluster density (clusters
> > 0.1 V_particle, per µm³) orders fine < medium < coarse.
> **(iii) recorded, never gating —** raw `n_clusters`, with the standing note
> that it is not comparable to the real value.

**A2 — G1-h is re-specified.** The p50 EDT statistic fails on voxel
quantization, not on direction (Step 1: 40.00 / 40.00 / 44.72 nm, tie; mean EDT
42.41 / 43.72 / 48.65 nm, correctly ordered). Same failure mode as Project 1's
degenerate D2.

> **G1-h (v2).** YSZ length scale orders fine < medium < coarse on **mean YSZ
> EDT** (a non-snapped statistic). EDT p50 is recorded but not gating.

**A3 — G1-a is re-specified.** `neck_scale` frozen per analog at the seed-0
solve caused 1 failure in 15 (coarse seed 4, +3.28 %).

> **G1-a (v2).** `neck_scale` is solved **per structure seed** to hit the analog's
> Φ_Ni target; the gate remains ±2 %. The solved value is recorded per structure.

**A4 — the YSZ generator is specified as particulate with an explicit sintering
yield**, replacing `add_ysz_pore`'s thresholded random field for Project 2. Two
parameters per analog, both frozen before any damage run: grain diameter
(scaling with the analog) and **`p_sinter`**, targeted at the §5 values.
**Geometric similarity must be broken:** relative neck size x/R decreases with
grain size, per Herring scaling. `add_ysz_pore` is **not modified or deleted** —
Project 1 results depend on it.

**A5 — an anti-tuning rule, carried from Project 1.** `p_sinter` is set once per
analog from the §5 derivation and frozen before any operator runs. It may be
re-derived only if the kill test shows the lattice-percolation mapping is
quantitatively wrong, and any such re-derivation is recorded as an amendment
with its reason. **`p_sinter` may never be adjusted to improve a damage result.**

---

## 7. The refined kill test — specification

One fine seed and one coarse seed. Pass requires **all four**:

| # | criterion | target |
|---|---|---|
| K1 | Φ_YSZ on target for both analogs | within ±2 % of real (0.421 / 0.384) |
| K2 | **Decoupling demonstrated:** at *fixed* grain diameter, varying `p_sinter` moves P_span | monotone, and spanning ≥ 0.90 → ≥ 0.999 |
| K3 | **Coarse hits the real P_span** at its derived `p_sinter` ≈ 0.42 | P_span ≈ 0.92 (accept 0.88–0.96) |
| K4 | **Fine hits the real P_span** at its derived `p_sinter` ≈ 0.70 | P_span ≥ 0.99 |
| K5 | filtered cluster density rises fine → coarse | secondary, reported |

**K2 is the actual kill test** — it is the claim that grain size and sintering
yield are independent controls. If P_span cannot be moved across the required
range at fixed grain size without breaking Φ_YSZ, Option A fails and we fall
back to **Option B (Ni-only)**.

Explicitly *not* required: matching real `n_clusters`, matching any Δ magnitude,
or any numerical fit — the project's ordinal-only constraint applies unchanged.

---

## 8. What I have and have not done

**Done:** the cluster-size audit (new measurement on the real stacks), the
ceiling recomputation, the physical-density correction to my own Step 1 report,
and the bond-percolation derivation of `p_sinter` targets.

**Not done — awaiting approval:** no YSZ generator code has been written. No
line of Option A exists. `cmlib/damage.py`, `cmlib/synth.py` and
`cmlib/project2.py` are unchanged by this response; the only new code is the
read-only audit script. O1/O2/O3 remain unimplemented.

**Recommendation:** adopt amendments A1–A5 and authorize the §7 kill test.
