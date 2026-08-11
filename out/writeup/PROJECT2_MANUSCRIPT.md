# Validation criteria that exclude the mechanism they validate: six artifact classes in voxel-based models of Ni-YSZ electrode degradation

**Draft manuscript, 2026-08-11.** Every number traces to a committed CSV and a
commit hash (§7). Claims needing literature support are marked
**[UNVERIFIED — CUT BEFORE SUBMISSION]**; no source access was available and
Amendment A5 binds: unsourceable claims are cut, not softened.

**Target: *Acta Materialia*, leading with Part B.** The artifact classes are
journal-independent and are the strongest contribution. Secondary: *Journal of
Power Sources*, leading with Part A.

---

## Corrections to the brief, before anything else

**1. The C3 values in the brief do not trace.** The brief states TPB retention
"0.799 > 0.621 > 0.574". Committed measurements:

| source | fine | medium | coarse | monotone? |
|---|---|---|---|---|
| this work (`phase6_comparison_table.csv`) | 0.7993 | 0.7460 | 0.5897 | yes |
| published, digitized | 0.7434 | 0.5862 | 0.6075 | **no** |
| brief | 0.799 | *0.621* | *0.574* | — |

The brief's medium and coarse match neither. **More importantly the published
series is not monotone** — coarse (0.6075) exceeds medium (0.5862). "Fine
retains TPB best" is robust in both sources; the full three-level C3 ordering
holds only in our measurement. **The paper claims the two-level fact only**,
mirroring the Ni side where medium-vs-coarse is likewise unresolved.

**2. "Six mechanisms tested" overstates the record.** Items 4 and 6 of the
brief's Part A are operator failures, not mechanism tests: O5 was an invalid
operator (it roughened rather than agglomerated) and O5v2 never reached a
result. **Four mechanisms were tested; two operator routes closed without
testing their mechanism.**

---

## Abstract

Voxel-based models are a standard tool for simulating microstructural
degradation in solid oxide cell electrodes. We report that two principal
geometry-based classes of such model fail to reproduce the degradation ordering
measured on real Ni-YSZ tomography, and — more consequentially — that six
methodological artifacts silently corrupt this class of study. The most general
is an impossibility: **enforcing monotone surface-area reduction, the natural
way to validate a coarsening operator, structurally excludes Rayleigh-type neck
break-up, the mechanism such operators are built to represent.** We demonstrate
it by bracketing: a curvature-ranked operator that reproduces neck thinning
violates area monotonicity at its first step, while a greedy area-decreasing
operator accepts zero moves because a sphere–neck–sphere geometry is already a
single-voxel local minimum. Finite-temperature acceptance bridges them and is
exactly what an area-monotonicity gate forbids. We further show that
connectivity metrics are pruning-dependent; that voxel-scale erosion
manufactures triple-phase boundary; that demanding perfect pristine connectivity
is over-constrained relative to real electrodes, which carry 1.2–11.2 %
disconnected Ni; that a jittered regular lattice has a minimum cut of exactly
one full cross-section with zero seed variance; and that curvature-ranked moves
do not guarantee area reduction. Real Ni networks fail at 1–2 % of throats via
small, spatially variable, non-planar cuts. All work was pre-registered, with
every amendment committed before the run it governed.

---

## Part B — Six artifact classes (lead contribution)

### B6. Monotone area reduction excludes Rayleigh-type break-up *(centrepiece)*

**Claim.** A coarsening operator validated by requiring monotone surface-area
reduction cannot exhibit Rayleigh-type neck break-up, because that instability
is collective and proceeds through a positive-area barrier.

Exact discrete identity for a one-voxel swap on a 6-connected lattice, derived
in `cmlib/damage2.py`: **ΔA = 2·(nb(a) − nb(b))**, where `nb` counts Ni
neighbours of the removed surface voxel `a` and the added pore site `b`.

| operator | ΔA test | neck (63 vox) | S_spec vs pristine 0.45052 | acceptance |
|---|---|---|---|---|
| Option A, 6-conn | none | 63 → 63 (no thinning) | 0.45195 ↑ | — |
| Option A, 26-conn | none | **63 → 57 → 15 → 0** | 0.45195 ↑, then 0.44431 ↓ | — |
| Option B, greedy | ΔA ≤ 0 | 63 → 63 | **0.45052 unchanged** | **0.000** |

Option A crosses the barrier without pricing it; Option B prices it and refuses
to cross. Greedy acceptance is exactly zero because in a sphere–cylinder–sphere
geometry the least-convex surface voxel still has more Ni neighbours than the
most-concave pore site (`nb_a,min > nb_b,max`) — **the structure is already a
single-voxel local minimum.** A finite-temperature rule bridges the two and is
what an area-monotonicity gate forbids by definition.

**Consequence.** Any voxel coarsening study that both enforces monotone area
reduction and claims Rayleigh-type break-up is either not producing that
break-up, or violating its own validity criterion. **[UNVERIFIED — CUT BEFORE
SUBMISSION: Rayleigh–Plateau instability of cylindrical interfaces; that area
monotonicity is commonly used as a validity check.]**

**What hid it:** on synthetic lattices the operator was never asked to thin a
neck against an explicit area budget.

### B1. Pruning-dependence of connectivity metrics
`P_span` = spanning-cluster voxels ÷ *total phase* voxels. Operators keeping
only the largest component delete exactly the non-spanning voxels, so
`P_span → 1.0000` identically — the operator rewrites the denominator. Measured:
pristine 0.9821 / 0.9713 / 0.8878 → **exactly 1.0000** at n = 1, for both O1 and
O6. Fixed by `R_Ni` = spanning ÷ *pristine* voxels, pruning-invariant,
sanity-checked to 0.00e+00. **Hidden by** gate G1-c requiring pristine
`P_span` = 1.0000.

### B2. TPB manufacture by voxel-scale erosion
Stochastic surface removal pits the surface and creates triple line. TPB
retention at n = 1: **15.24 / 12.72 / 7.71** (O1); 14.80 / 12.27 / 7.36 (O6).
**Hidden by** a synthetic platform whose TPB was already ~8× real.

### B3. Over-constrained pristine connectivity
Real electrodes carry **1.2–11.2 %** disconnected Ni when pristine (per-ROI
`P_span` 0.9821/0.9754/0.9877; 0.9713/0.9446/0.9554; 0.8878/0.9528/0.9157). A
gate demanding exactly 1.0000 is unrepresentative of every real electrode
measured, and it is what concealed B1.

### B4. Lattice min-cut planarity
A jittered regular lattice has a minimum S–T cut of **exactly one full
cross-section — 36 = 6², 25 = 5², 16 = 4² — with zero seed-to-seed variance**
across 15 structures. Random and low-area attacks never break spanning even at
50 % neck removal. Real networks: min-cut **1–2 %** of throats, ROI-variable,
non-planar. **A regular lattice cannot fail the way a real network fails.**

### B5. Curvature-ranked moves do not guarantee area reduction
Curvature rank is a *proxy* for ΔA, not ΔA. On a discrete lattice an added voxel
exposes new faces at its free sides; at n = 1 that transient exceeded the area
removed (0.45195 vs 0.45052) under both 6- and 26-connectivity.

---

## Part A — Systematic elimination

**Four mechanisms tested; two operator routes closed without reaching a
mechanism test.**

| # | mechanism | outcome | key number |
|---|---|---|---|
| 1 | lower-tail neck widening | **no lever** | 8.54 / 8.50 / 8.50 rounds; −0.04, ~25× below threshold |
| 2 | narrow-neck severing (O2) | **no failure** | 100 % of lower-quartile throats severed, `P_span` = 1.0000 in 15/15, 1.8 % volume lost |
| 3 | YSZ fracture (O3) | **pristine-loaded, not mechanistic** | main +2.87 rounds → matched-fragility control **−0.27** |
| 4 | surface erosion, real ROIs (O6) | **reverses the real ordering** | fine 9.72 > medium 9.06/9.28, coarse 9.39; all ρ ≤ +0.213 |
| — | volume-conserving redistribution (O5) | operator invalid (surface +3.3 %) | mechanism not tested |
| — | volume-conserving coarsening (O5v2) | route closed, see B6 | mechanism not tested |

Mechanisms 1 and 2 are mutually reinforcing: widening the narrow-neck population
gives no benefit *and* destroying it entirely causes no failure. The min-cut
audit explains both — critical edges overlap the lower quartile at
**0.250 / 0.280 / 0.188** against a chance expectation of 0.25. **The narrow
necks are not load-bearing.**

---

## Part C — Positive structural findings

1. **Real Ni networks fail at 1–2 % of throats** via small, ROI-variable,
   non-planar cuts (fine 0.0141–0.0196; medium 0.0221–0.0302; coarse
   0.0055–0.0067); a jittered lattice fails only at a full cross-section.
   Targets for any future disordered generator: mean degree **3.4–3.9**,
   sd **2.0–2.4**, connected-pair-distance CV **≈ 0.41**.
2. **Pristine disconnection is real and coarseness-ordered** (B3).
3. **The coarse fragility paradox.** Coarse is 3–5× the most topologically
   fragile yet the best real retainer (0.947). **Pristine topology does not
   predict degradation** — and neither does specific surface area (all
   |ρ| ≤ 0.213, raw and partial).
4. **The corrected TPB estimator is sound on real data**: 4.48–4.66 /
   1.87–2.37 / 1.47–1.55 µm⁻² against literature 3.62 / 2.11 / 1.47 —
   **1.0–1.3×**. The ~8× inflation seen earlier was a property of the synthetic
   YSZ placement, not the estimator.

---

## What this paper does NOT claim

- **Not** that agglomeration is not the degradation mechanism. It may be; we
  could not test it (B6).
- **Not** that no geometry-based mechanism can explain the ordering.
  Electrochemical and transport-coupled mechanisms are untouched.
- **Not** that phase-field or continuum methods would fail — they are not
  subject to B6's single-voxel argument. **Future work.**
- **Not** a three-level TPB ordering (see Corrections).

---

## Limitations

3 ROIs and 3 damage seeds per anode; two `R_Ni` thresholds that are not
independent because the spanning cluster vanishes discontinuously; real
pristine and degraded stacks are **different specimens**, so every retention
value carries specimen confounding; synthetic analogs are ~3× size-compressed;
n = 3 anodes throughout, so **no p-values are reported anywhere in this work, by
design**; a `porespy.snow2` chunked-watershed bug was active in the original
Phase 3/4 extractions (assessed; no conclusion changed).

---

## §7 Traceability, pre-registration chain, self-corrections

**Every amendment was committed before the run it governed.**
`e62f30b` (prereg v2.1) → `8e77035` (kill test) → `52aa519` (Step 1 requal) →
`b278dc5` (Step 2) → `cb2ca49`→`971f0ba` (O5) → `267efd3`→`58bc6db` (vuln
audit) → `50b0bde`→`c41e4fa` (v3) → `49b1059`→`0aa6276` (gates) →
`f980562`→`ff3e5d1` (O6) → `bfc5d89`→`10f2c51` (R_Ni) → `307a220` (C1-real) →
`304aa8a`→`c38218d` (O5v2-A) → `c414e63`→`84bcc61` (O5v2-B).

**Number → source.** A1 `p10_group_experiment.csv` (`4bd7487`); A2
`step2_O2_ceiling_check.csv`; A3 `step2_O3_main.csv` + `step2_O3_control.csv`
(`b278dc5`); A4 `c1real_results.csv` (`307a220`); B1/B3 `c1real_o6_validity.csv`,
`c1real_rni_gate.csv` (`ff3e5d1`, `10f2c51`); B2 `step2_O1_main.csv`; B4
`audit_ni_vulnerability.csv` (`58bc6db`); B5/B6 `o5v2_report.md`,
`o5v2_optionB_report.md` (`c38218d`, `84bcc61`); C1 `c1real_results.csv` +
`v3_real_ni_graph_audit.csv`; C4 `step2_r2_tpb_baseline.csv`,
`c1real_rni_gate.csv`; real outcomes `phase6_comparison_table.csv` (`47b08ee`).

**Self-corrections against our own work**, part of the record:
1. **TPB periodic wrap** — `np.roll` over-counted by exactly 4× on the analytic
   case; fixed by slicing.
2. **Reversed axis extents** — `SAMPLES` stores `(nx,ny,nz)` while the array is
   `(z,y,x)`; v3 min-cut fractions were provisional until re-run.
3. **Saturated `n_secondary`** — the frozen rule chose a point past the
   transition; reported, not re-picked.
4. **A verdict line that could not detect its own failure** — the A1 check
   compared only across damage intensities, never against pristine, and printed
   PASS on a failing run.
5. **A contradictory self-prompt** — "Metropolis" and "accept only
   area-non-increasing" cannot both hold; resolved to greedy.

**Controls that caught would-be false positives:** the R3 matched-fragility
control (turned a +2.87-round C2 "pass" into −0.27); Project 1's provenance
verification of the widened structures (67,803 voxels differing; 46/224 necks
changed) before its null was accepted.

## Literature verification — BLOCKING

No source access was available. **The manuscript cannot be submitted until every
[UNVERIFIED] marker is source-checked or cut.** Outstanding: Rayleigh–Plateau
instability; Simwonis/Tietz/Stöver 2000; Sarantaridis & Atkinson 2007;
Faes/Grew FIB-SEM coarsening; Herring scaling; NiO/Ni molar volumes (≈11.2 vs
6.59 cm³/mol); SC/FCC/BCC percolation thresholds (0.2488 / 0.1803 / 0.120) and
the SC site threshold 0.3116077 (Xu et al. 2014) used in the Phase-0 gate.
Per A5: **cut, do not soften.**
