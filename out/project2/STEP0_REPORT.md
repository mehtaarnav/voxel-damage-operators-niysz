# Project 2 — Step 0 result and Go/No-Go decision

Run 2026-08-10. Code: `scripts/project2/step0_percolation.py`,
`scripts/project2/step0_figures.py`. Data: `out/phase5/phase5_percolation.csv`,
`out/phase5/phase5_retention.csv`, `out/project2/step0_ysz_percolation.csv`,
`out/project2/step0_ysz_retention.csv`, `out/project2/step0_percolation.png`.

Full stacks (0.48–1.11 Gvoxel), native resolution, no sub-volume.
6-connectivity, transport axis = x = array axis 2 — all inherited from Phase 0,
none re-chosen.

**Headline: the gate returns outcome (a), and more sharply than the brief's
premise claimed. The coarse degraded anode's YSZ network does not percolate at
all — `P_span = 0.0000`, largest component 12.4 % of the phase.** The divergent
signature Project 2 set out to fingerprint is real, measured, and large.

---

## 1. Phase 5 defect resolution

### 1.1 Corrected root cause

My earlier note (manifest §2.2) inferred an "interrupted write". **That
inference was wrong and is corrected here.** The committed artifact pattern —
`phase5_percolation.csv` with exactly one row (`coarse_pre`),
`phase5_retention.csv` containing only a header-less empty frame, *and a valid
43 KB `phase5_percolation.png`* — is precisely what `phase5_percolation.py`
produces when invoked as `--samples coarse_pre`. A crash would not have written
the figure. **The most likely cause is a deliberate single-sample invocation
whose output was then committed as if it were the full run**, not a failure.
The practical consequence for a reproducer is unchanged, and no reported
conclusion was ever affected.

### 1.2 Both files regenerated, complete

All six stacks, original column schema preserved exactly:

| sample | grain | state | n_clusters | spans x | P_span | P_reach | P_largest | published *P* | P_reach/pub |
|---|---|---|---|---|---|---|---|---|---|
| fine_pre | fine | pristine | 2,526 | yes | 0.9840 | 0.9865 | 0.9840 | 0.985 | 1.001 |
| medium_pre | medium | pristine | 2,774 | yes | 0.9583 | 0.9693 | 0.9583 | 0.965 | 1.004 |
| coarse_pre | coarse | pristine | 1,225 | yes | 0.9315 | 0.9668 | 0.9315 | 0.959 | 1.008 |
| fine_post | fine | degraded | 3,545 | yes | 0.6686 | 0.8457 | 0.6686 | 0.809 | 1.045 |
| medium_post | medium | degraded | 1,594 | yes | 0.8190 | 0.9490 | 0.8190 | 0.884 | 1.074 |
| coarse_post | coarse | degraded | 1,949 | yes | 0.8821 | 0.9105 | 0.8821 | 0.886 | 1.028 |

**Cross-check against the frozen Phase 6 table — exact agreement:**

| grain | phase5 `P_span_retained` | phase6 `P_span_retained` | difference |
|---|---|---|---|
| fine | 0.679503 | 0.679503 | 0.00 × 10⁰ |
| medium | 0.854693 | 0.854693 | 0.00 × 10⁰ |
| coarse | 0.946965 | 0.946965 | 0.00 × 10⁰ |

`P_reach/published` spans 1.001–1.074, reproducing REPORT.md's stated
"0.1–7.4 %" exactly. **The defect was regenerability-only, as claimed; the
regenerated numbers confirm every value Project 1 reported.** The stale
single-sample `phase5_percolation.png` has also been replaced.

### 1.3 Durable codebase fixes

1. **`cmlib.percolation.percolation_summary_lowmem`** — the label array now
   lives in a disk-backed memmap and is reduced slab-wise. This was necessary,
   not precautionary: the largest stack (coarse_post, 1.115 Gvoxel) needs a
   4.46 GB int32 label array on a machine with 2.2 GB free at run time.
   Definitions are unchanged; **equivalence against the frozen
   `percolation_summary` is asserted on 14 randomised volumes × 3 axes
   (spanning p = 0.05 … 0.8, plus all-empty and all-full) before any stack is
   read.** The gate passed, and the real-data agreement with Phase 6 above is a
   second, independent confirmation.
2. **Incremental CSV writes** — every output file is rewritten after each
   sample, so an interruption now truncates to a valid short file. Combined
   with 1.3.1 this removes both the mechanism I originally suspected and the
   one that actually applied.

---

## 2. Step 0 — YSZ percolation, measured for the first time

| sample | grain | state | Φ_YSZ | n_clusters | **spans x** | **P_span** | P_reach | **P_largest** |
|---|---|---|---|---|---|---|---|---|
| fine_pre | fine | pristine | 0.4213 | 975 | yes | 0.9989 | 0.9991 | 0.9989 |
| medium_pre | medium | pristine | 0.3881 | 1,598 | yes | 0.9880 | 0.9884 | 0.9880 |
| coarse_pre | coarse | pristine | 0.3838 | 2,917 | yes | 0.9246 | 0.9332 | 0.9246 |
| fine_post | fine | degraded | 0.3123 | 2,671 | yes | 0.9568 | 0.9632 | 0.9568 |
| medium_post | medium | degraded | 0.3764 | 2,893 | yes | 0.8549 | 0.8908 | 0.8549 |
| **coarse_post** | **coarse** | **degraded** | **0.3243** | **3,873** | **NO** | **0.0000** | 0.5029 | **0.1235** |

### Retained YSZ percolation

| grain | P_span pre → post | **P_span retained** | P_largest retained | n_clusters ratio |
|---|---|---|---|---|
| fine | 0.9989 → 0.9568 | **0.958** | 0.958 | 2.74 |
| medium | 0.9880 → 0.8549 | **0.865** | 0.865 | 1.81 |
| **coarse** | 0.9246 → **0.0000** | **0.000** | **0.134** | 1.33 |

---

## 3. What the measurement establishes

### 3.1 The divergence is real and ordinally clean

| grain | **Ni** P_span retained | **YSZ** P_span retained |
|---|---|---|
| fine | **0.680** ← worst | 0.958 |
| medium | 0.855 | 0.865 |
| coarse | 0.947 | **0.000** ← worst |

**Fine is worst at retaining Ni percolation; coarse is worst at retaining YSZ
percolation — completely so.** Both orderings are strict and three-level on the
YSZ side. This is exactly the signature the brief described, now held as a
measurement rather than a paraphrase.

### 3.2 The YSZ collapse is *not* a volume-fraction effect

This is the strongest single finding, and it kills the obvious confound before
it can be raised:

- **fine_post has the LOWEST degraded Φ_YSZ (0.3123) and still spans at 0.957.**
- **coarse_post has MORE YSZ (0.3243) and does not span at all.**

More YSZ, no percolation; less YSZ, robust percolation. Whatever destroys the
coarse YSZ backbone, it is **morphological, not volumetric**. The random-site
threshold (0.3116) is plotted in `step0_percolation.png` for reference and is
plainly not the operative variable — fine_post sits essentially *on* it and
spans fine.

### 3.3 My §0.2 concern is resolved — and resolved against volume fraction

The design memo flagged that ΔΦ_YSZ suggested *fine* was worst (−0.109 vs
coarse −0.060). **That reading was wrong, and it was wrong in an instructive
way:** fine loses the most YSZ *volume* and retains YSZ *percolation* best;
coarse loses less volume and loses percolation entirely. **Volume fraction is
anti-predictive of percolation retention for YSZ — the same failure mode
Project 1 documented for Ni**, where Φ_Ni ranked the anodes exactly backwards.
Two phases, two datasets, same lesson: volume fraction does not predict
connectivity retention in this system.

### 3.4 Coarse YSZ starts weakest, then fails completely

Pristine YSZ P_span orders **fine 0.9989 > medium 0.9880 > coarse 0.9246**, with
coarse carrying 3× the component count of fine (2,917 vs 975) *before any
degradation*. So the coarse YSZ backbone is the most marginal to begin with and
then collapses outright. This has a direct design consequence — see §5.2.

### 3.5 It is a spanning failure, not disintegration

`coarse_post` retains `P_reach = 0.503`: half the YSZ phase still touches one
face. The network did not crumble uniformly; it **severed**. The largest
component fell to 12.4 % while the cluster count rose only 1.33× — the smallest
fragmentation-count increase of the three anodes. **A few decisive breaks, not
pervasive comminution.** This is precisely the "fragmentation at near-constant
volume" mechanism O3a was designed around, and it vindicates two design choices
in the memo: omitting island-removal for YSZ, and recording `P_largest` and
`n_clusters` alongside the spanning boolean.

---

## 4. Gate evaluation

| outcome | criterion | verdict |
|---|---|---|
| **(a) coarse worst** | coarse worst on YSZ retention | **YES — 0.000 vs 0.865 vs 0.958** |
| (b) fine worst | would kill Project 2 | NO |
| (c) saturated | all six ≥ 0.99 | NO — spans 0.0000 to 0.9989 |

**Outcome (a). The pilot proceeds as written in DESIGN_MEMO §4.** No
substituted YSZ outcome variable is required: **binary YSZ percolation is a
live, sensitive outcome on real data** — it goes to exactly zero in the one case
the mechanism predicts. `P_largest` and `n_clusters` are retained as secondary
recorded diagnostics, not promoted to primary.

### C2 stands as written, with its direction now fixed by measurement

> **C2 — YSZ ordering.** Coarse loses YSZ percolation at a **strictly lower**
> damage intensity than fine: `n*_YSZ(coarse) < n*_YSZ(fine)`, by **≥ 1.0
> damage round**.

Given the real data's three-level YSZ ordering (unlike Ni, where medium vs
coarse is unresolved), C2 may additionally be scored at three levels as a
**secondary, non-binding** check: `n*_YSZ(coarse) < n*_YSZ(medium) <
n*_YSZ(fine)`. Failing the three-level version does not fail C2.

---

## 5. Two amendments to the design memo, forced by the data

Both are consequences of §3.4 and must be frozen before Step 1 runs.

### 5.1 Gate G1-d was wrong and would have rejected a correct analog

DESIGN_MEMO §4.2 set **G1-d: YSZ percolates, P_span ≥ 0.99**. Real
`coarse_pre` YSZ has **P_span = 0.9246** — a faithful coarse analog would be
*rejected by this gate for correctly reproducing reality*. Replace with:

> **G1-d (revised):** YSZ percolates (`P_span > 0`) in every analog at t = 0,
> and the pristine YSZ `P_span` is recorded per seed. No fixed threshold.

### 5.2 New gate G1-i: reproduce the pristine YSZ ordering

The measurement shows pristine YSZ quality is itself coarseness-ordered. An
analog set that starts all three classes at the same YSZ robustness is not
faithful, and would hand O3 an easier problem than reality poses.

> **G1-i (new):** pristine YSZ `P_span` must order **fine > medium > coarse**,
> and pristine YSZ `n_clusters` must order **fine < medium < coarse**.

This strengthens the §2.2/A2.3 requirement that YSZ morphology scale with
analog coarseness (the σ-scaling addition): G1-i is now the *measured* target
that σ-scaling must hit, not a plausibility argument. If σ-scaling cannot
produce this ordering, the YSZ generator — not the operator — is the thing that
needs redesign, and we will know that before any operator is written.

---

## 6. Limitations of Step 0

Stated now so they travel with every downstream use of this signature.

1. **Pristine and degraded are different specimens** (Rx36/37/38 vs
   Rx41-1/2/3), not the same volume re-imaged. Every retention number is
   therefore confounded with specimen-to-specimen variation. **This is the
   single most important caveat.** It is mitigated, not removed, by the effect
   size: a bulk topological measure going from 0.925 to *exactly zero* is far
   outside any plausible specimen-to-specimen range for a quantity that sits at
   0.92–0.99 in all five other stacks. Ordinal use is defensible; magnitude use
   is not.
2. **n = 1 per (grain, state) cell.** No replicates, so no error bars and no
   statistics. Consistent with Project 1's standing rule: no p-values are
   computed anywhere.
3. **Segmentation is taken as given.** Label assignment was verified through
   published volume fractions (worst deviation 0.20 % over 18 values, Phase 2
   gate), including `coarse_pre`'s inverted Ni/YSZ grey convention. YSZ
   percolation is nonetheless sensitive to thin-feature segmentation in a way
   volume fraction is not; a one-voxel erosion of the YSZ label would lower
   P_span, and the coarse anode has the thinnest YSZ features. **Not tested
   here.** A segmentation-sensitivity probe is the cheapest available check on
   §3.1 and is recommended before the Project 2 result is published — it is not
   required before the pilot.
4. **Anisotropy.** Post-redox medium and coarse stacks are 40 %
   voxel-anisotropic (17.9 × 17.9 × 25 nm). Percolation is purely topological
   and so is unaffected — this caveat bites SNOW-derived metrics, not these.
5. **One axis.** Only the x (through-thickness) transport axis is scored;
   `percolates_y`/`percolates_z` are recorded in the CSVs but not analysed.

---

## 7. Decision

**GO for Step 1 (analog qualification), with G1-d revised and G1-i added.**

- Step 0: **complete**, outcome **(a)**.
- Phase 5 defect: **resolved**, with the root cause corrected and two durable
  code fixes.
- Target signature for Project 2 is now **measured**: Ni retention
  0.680 / 0.855 / 0.947 (fine worst), YSZ retention 0.958 / 0.865 / 0.000
  (coarse worst).
- **Not run, per instruction:** O1, O2, O3 remain unimplemented; no synthetic
  analog has been built; `cmlib/damage.py` and `cmlib/synth.py` are untouched.

**Awaiting authorization for Step 1.**
