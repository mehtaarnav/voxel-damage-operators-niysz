# When the validation criterion excludes the mechanism: six artifact classes in voxel-based models of Ni–YSZ electrode degradation

**Manuscript draft — 2026-08-11**
Target: *Acta Materialia* (Part B lead) · Alternate: *Journal of Power Sources* (Part A lead)

Every quantity traces to a committed CSV and a git commit hash (§8). References
carry an explicit verification status (§9); nothing was checked against a source
in this environment, so items marked ⚠ or ✂ must be resolved before submission.

---

## Abstract

Voxel-based models are a standard tool for simulating microstructural
degradation in solid oxide cell electrodes. Working from FIB-SEM tomography of
three Ni–YSZ anodes spanning a coarseness series, we report that two principal
geometry-based classes of such model fail to reproduce the measured degradation
ordering, and — more consequentially — that six methodological artifacts
silently corrupt this class of study.

The most general is an impossibility. **Enforcing monotone surface-area
reduction, the natural way to validate a coarsening operator, structurally
excludes Rayleigh-type neck break-up, the very mechanism such operators are
built to represent.** We establish it by bracketing rather than by a single
experiment: a curvature-ranked operator reproduces neck thinning (63 → 0 voxels)
but violates area monotonicity at its first step, while a greedy
area-decreasing operator accepts *zero* moves because a sphere–neck–sphere
geometry is already a single-voxel local minimum. Finite-temperature acceptance
bridges the two and is exactly what an area-monotonicity gate forbids.

We further show that connectivity metrics are pruning-dependent; that
voxel-scale erosion manufactures triple-phase boundary before destroying it;
that demanding perfect pristine connectivity is over-constrained relative to
real electrodes, which carry 1.2–11.2 % disconnected Ni; that a jittered regular
lattice has a minimum cut of exactly one full cross-section with zero seed
variance; and that curvature-ranked moves do not guarantee area reduction. Real
Ni networks, by contrast, fail at 1–2 % of throats via small, spatially
variable, non-planar cuts.

All work was pre-registered, with every amendment committed to version control
before the run it governed.

---

## 1. Introduction

Ni–YSZ cermet anodes lose performance through microstructural change: the Ni
phase coarsens and can lose electronic percolation, the YSZ scaffold can
fracture, and the triple-phase boundary (TPB) where all three phases meet is
consumed. Simulating this on segmented tomograms — applying a voxel-scale
"damage operator" and tracking connectivity — is an established approach ⚠[R7,R8].

The measured target signature is specific and counter-intuitive. Across the
Holzer/Pecho anode series ✓[R1,R2], the **fine** anode retains Ni percolation
*worst* yet retains TPB *best*:

| anode | Ni percolation retained | TPB retained (this work) | TPB retained (published) |
|---|---|---|---|
| fine | **0.6795** | **0.7993** | **0.7434** |
| medium | 0.8547 | 0.7460 | 0.5862 |
| coarse | 0.9470 | 0.5897 | 0.6075 |

*Source: `phase6_comparison_table.csv` (`47b08ee`).*

Two cautions travel with this table for the rest of the paper. **Medium versus
coarse is unresolved on the Ni side** — it flips with the definition of the
percolation measure. And the **published TPB series is not monotone** (coarse
0.6075 > medium 0.5862), so only the two-level statement *"fine retains TPB
best"* is defensible; the three-level ordering holds in our measurement alone.

We set out to identify which damage mechanism reproduces this signature. We did
not succeed, and the reasons why are the contribution.

## 2. Methods in brief

**Data.** Six segmented FIB-SEM stacks (three anodes × pristine/degraded),
0.48–1.11 Gvoxel each, at 19.53–29.14 nm voxel pitch ✓[R1,R2]. Pristine and
degraded are **different specimens**, which confounds every retention magnitude
and is why all claims here are ordinal.

**Percolation.** 6-connectivity, free boundaries, face-to-face spanning. The
implementation reproduces the simple-cubic site threshold to 0.31218 against a
reference 0.3116077 ✓[R3].

**Networks.** SNOW watershed partitioning ✓[R4] on the solid phase; nodes are Ni
chambers, edges are throats with measured inscribed diameter.

**Outcome metric.** After discovering artifact B1 we adopted
**R_Ni = spanning-cluster Ni voxels ÷ pristine Ni voxels**, which is invariant
to island pruning and directly comparable to measured retention. Its sanity
check — R_Ni(0) must equal pristine P_span — holds to 0.00 × 10⁰.

**Statistics.** n = 3 anodes. **No p-values are reported anywhere in this work,
by design.** Correlations use Spearman ρ with leave-one-out sign checks, and a
|ρ| ≥ 0.6 threshold fixed in advance.

---

## Part B — Six artifact classes

### B6. Monotone area reduction excludes Rayleigh-type break-up

**Claim.** A coarsening operator validated by requiring monotone surface-area
reduction cannot exhibit Rayleigh-type neck break-up, because that instability
is collective and proceeds through a positive-area barrier.

For a one-voxel swap on a 6-connected lattice the area change is exact, not
approximate. Removing surface voxel *a* changes exposed Ni faces by
2·nb(a) − 6; adding at pore site *b* by 6 − 2·nb(b), where nb counts Ni
neighbours. Hence

> **ΔA = 2·[nb(a) − nb(b)]**,  so  **ΔA ≤ 0 ⟺ nb(a) ≤ nb(b)**.

![Figure 1](figs/fig1_b6_impossibility.png)

**Figure 1.** Two operators bracket an impossibility on a sphere–neck–sphere
test body (neck 63 voxels, S_spec = 0.45052). **(a)** The curvature-ranked
operator with a 26-connectivity stencil thins the neck to zero; with
6-connectivity it does nothing; the greedy ΔA ≤ 0 operator does nothing.
**(b)** The same runs in surface area. The curvature-ranked operator enters the
region forbidden by validity gate (ii) at n = 1 (0.45195 > 0.45052) before
falling below pristine at n = 3. The greedy operator never leaves the pristine
line — it accepts **zero** moves.
*Source: `o5v2_report.md`, `o5v2_optionB_report.md` (`c38218d`, `84bcc61`).*

Greedy acceptance is exactly zero because in this geometry the least-convex
surface voxel still has more Ni neighbours than the most-concave pore site
(nb_a,min > nb_b,max): **the body is already a single-voxel local minimum.**
Spheres minimise area at fixed volume, and a straight cylinder is stable to any
*single* voxel move — Rayleigh break-up requires a correlated set of voxels to
move together, transiently raising area ⚠[R5,R6].

The two options therefore bracket the space of single-swap algorithms: Option A
crosses the barrier without pricing it; Option B prices it and refuses to cross.
The only bridge is finite-temperature acceptance, which an area-monotonicity
gate forbids by construction.

**Consequence for practice.** Any voxel coarsening study that both enforces
monotone area reduction and claims Rayleigh-type break-up is either not
producing that break-up, or is violating its own validity criterion.

### B1. Pruning-dependence of connectivity metrics

P_span divides spanning-cluster voxels by *total phase* voxels. Operators that
retain only the largest connected component delete exactly the non-spanning
voxels — so P_span → 1.0000 identically. The operator rewrites the denominator.

Measured on real ROIs: pristine P_span 0.9821 / 0.9713 / 0.8878 → **exactly
1.0000** at n = 1, for two independent operators. Resolved by R_Ni (§2).
**Hidden by** a synthetic qualification gate that required pristine
P_span = 1.0000.

### B2. TPB manufacture by voxel-scale erosion

![Figure 4](figs/fig4_tpb.png)

**Figure 4.** TPB relative to pristine under reduction-only surface erosion on
real ROIs. TPB rises 7–15× before collapsing. Stochastic single-voxel removal
pits the Ni surface, and every new Ni/pore facet touching YSZ creates triple
line. Any TPB conclusion drawn before the peak is an artifact of discretisation.
*Source: `c1real_rni_gate.csv` (`10f2c51`).*

### B3. Over-constrained pristine connectivity

Real electrodes carry **1.2–11.2 %** disconnected Ni when pristine (per-ROI
P_span: fine 0.9821/0.9754/0.9877; medium 0.9713/0.9446/0.9554; coarse
0.8878/0.9528/0.9157). Synthetic platforms that require perfect pristine
connectivity are unrepresentative of every real electrode measured here — and
that requirement is what concealed B1 for the entire synthetic phase of this
study.

### B4. Lattice min-cut planarity

![Figure 2](figs/fig2_mincut.png)

**Figure 2.** Minimum source–sink cut as a fraction of all necks. A jittered
regular lattice cuts at exactly one full cross-section — 36 = 6², 25 = 5²,
16 = 4² — with **zero** seed-to-seed variance across 15 structures. Real ROIs cut
at 1–2 % of throats, with genuine ROI-to-ROI spread. Random and low-area attacks
never break the lattice even at 50 % neck removal.
*Source: `audit_ni_vulnerability.csv` (`58bc6db`), `c1real_results.csv` (`307a220`).*

A regular lattice has no locally weak plane: every cross-section is
statistically identical. **It cannot fail the way a real network fails**, and no
local stochastic operator can make it.

### B5. Curvature-ranked moves do not guarantee area reduction

Curvature rank is a *proxy* for ΔA, not ΔA. On a discrete lattice an added voxel
exposes new faces at its own free sides; at n = 1 that transient exceeds the
area removed (0.45195 vs 0.45052) under both 6- and 26-connectivity stencils
(Figure 1b). The exact identity in B6 is the remedy.

---

## Part A — Systematic elimination

Four mechanisms were tested; two operator routes closed without reaching a
mechanism test.

| # | mechanism | outcome | key quantity |
|---|---|---|---|
| 1 | lower-tail neck widening | **no lever** | 8.54 / 8.50 / 8.50 rounds; −0.04, ≈25× below threshold |
| 2 | narrow-neck severing | **no failure** | 100 % of lower-quartile throats severed; P_span = 1.0000 in 15/15; 1.8 % volume lost |
| 3 | YSZ contact fracture | **pristine-loaded, not mechanistic** | +2.87 rounds → matched-fragility control **−0.27** |
| 4 | surface erosion, real ROIs | **reverses the real ordering** | Figure 3; all \|ρ\| ≤ 0.213 |
| — | volume-conserving redistribution | operator invalid (surface +3.3 %) | mechanism untested |
| — | volume-conserving coarsening | route closed (B6) | mechanism untested |

![Figure 3](figs/fig3_reversal.png)

**Figure 3.** **(a)** Measured Ni percolation retention: fine is worst.
**(b)** Transition intensity under reduction-only surface erosion applied to the
same real microstructures: **fine is last to fail.** The simulation inverts the
experiment. Neither specific surface area (ρ ≤ +0.089) nor pristine min-cut
fraction (ρ ≤ +0.213) correlates with transition intensity, raw or partial.
*Source: `c1real_results.csv` (`307a220`), `phase6_comparison_table.csv`.*

Mechanisms 1 and 2 are mutually reinforcing: widening the narrow-neck population
confers no benefit **and** destroying it entirely causes no failure. The min-cut
audit explains both — critical edges overlap the lower-quartile neck population
at 0.250 / 0.280 / 0.188, against a chance expectation of 0.25. **The narrow
necks are not load-bearing.**

Mechanism 3 illustrates the value of controls. The main arm showed a clean
2.87-round separation in the predicted direction — a publishable-looking result.
A matched-fragility control, in which pristine YSZ fragility was equalised
across analogs, **inverted it to −0.27 rounds**: the effect was pristine-state
loading, not mechanism.

---

## Part C — Positive findings

1. **Real Ni networks fail at 1–2 % of throats** via small, ROI-variable,
   non-planar cuts (Figure 2). Any future disordered generator must hit: mean
   degree 3.4–3.9, degree sd 2.0–2.4, connected-pair-distance CV ≈ 0.41.
2. **Pristine disconnection is real and coarseness-ordered** (B3).
3. **The coarse fragility paradox.** Coarse is 3–5× the most topologically
   fragile network yet the best retainer in reality (0.947). **Pristine topology
   does not predict degradation**, and neither does specific surface area.
4. **The corrected TPB estimator is sound on real data**: 4.48–4.66 /
   1.87–2.37 / 1.47–1.55 µm⁻² against published 3.62 / 2.11 / 1.47 — a factor
   of 1.0–1.3. The ≈8× inflation observed earlier was a property of synthetic
   phase placement, not of the estimator.

---

## 6. What this paper does not claim

- **Not** that agglomeration is not the degradation mechanism. It may well be;
  B6 says we could not test it with this class of operator.
- **Not** that no geometry-based mechanism explains the ordering.
  Electrochemical and transport-coupled mechanisms are untouched.
- **Not** that phase-field or continuum methods fail — they are not subject to
  B6's single-voxel argument, and are the natural next approach.
- **Not** a three-level TPB ordering (§1).

## 7. Limitations

Three ROIs and three damage seeds per anode; two R_Ni thresholds that are not
independent because the spanning cluster vanishes discontinuously; pristine and
degraded stacks are different specimens; synthetic analogs are ≈3× size
compressed with proportionally fatter necks, a bias that works *against* the
hypothesis under test; n = 3 anodes; a `porespy.snow2` chunked-watershed
defect was active in the earliest extractions (assessed, no conclusion changed).

---

## 8. Reproducibility, pre-registration and self-correction

**Every amendment was committed before the run it governed.** Chain:
`e62f30b` → `8e77035` → `52aa519` → `b278dc5` → `cb2ca49`/`971f0ba` →
`267efd3`/`58bc6db` → `50b0bde`/`c41e4fa` → `49b1059`/`0aa6276` →
`f980562`/`ff3e5d1` → `bfc5d89`/`10f2c51` → `307a220` → `304aa8a`/`c38218d` →
`c414e63`/`84bcc61`.

**Quantity → source.** A1 `p10_group_experiment.csv`; A2
`step2_O2_ceiling_check.csv`; A3 `step2_O3_main.csv` + `step2_O3_control.csv`;
A4 and Figures 2–3 `c1real_results.csv`; B1/B3 `c1real_o6_validity.csv`,
`c1real_rni_gate.csv`; B2 and Figure 4 `c1real_rni_gate.csv`; B4
`audit_ni_vulnerability.csv`; B5/B6 and Figure 1 `o5v2_report.md`,
`o5v2_optionB_report.md`; measured outcomes `phase6_comparison_table.csv`.
Figures regenerate via `scripts/project2/make_figures.py`.

**Self-corrections issued against our own work**, recorded as part of the
result:
1. A TPB estimator that over-counted by exactly 4× through periodic wrap.
2. Reversed axis extents that made one audit's cut fractions provisional.
3. A pre-registered secondary-outcome rule that selected a saturating point.
4. A validity check that compared only across damage intensities, never against
   pristine, and printed PASS on a failing run.
5. An internally contradictory plan specifying both Metropolis acceptance and
   strict area monotonicity.

**Controls that caught would-be false positives:** the matched-fragility control
(§Part A, mechanism 3); and a provenance verification of synthetic structures
(67,803 voxels differing, 46 of 224 necks changed) before an earlier null was
accepted.

---

## 9. References and verification status

**No reference below was checked against a source in this environment.**
✓ = citation is recorded in the repository from earlier verified work;
⚠ = from memory, plausible, **must be verified**; ✂ = cut unless verified.

| # | reference | status |
|---|---|---|
| R1 | Pecho et al., *Materials* — Ni–YSZ transport/percolation (PMC5512617, `ma8095265`) | ✓ in-repo |
| R2 | Holzer et al., *Materials* — Ni–YSZ TPB (PMC5455394, `ma8105370`) | ✓ in-repo |
| R3 | Xu, Wang, Lv & Deng, *Phys. Rev. E* **89**, 012120 (2014) — simple-cubic site threshold 0.3116077 | ✓ in-repo |
| R4 | Gostick, *Phys. Rev. E* **96**, 023307 (2017) — SNOW partitioning | ✓ in-repo |
| R5 | Rayleigh / Plateau — capillary instability of a liquid cylinder | ⚠ |
| R6 | Nichols & Mullins — surface-diffusion-driven break-up of solid cylinders, *J. Appl. Phys.* (1965) | ⚠ |
| R7 | Simwonis, Tietz & Stöver, *Solid State Ionics* **132** (2000) — Ni coarsening in Ni–YSZ | ⚠ |
| R8 | Sarantaridis & Atkinson, *Fuel Cells* **7** (2007) — redox cycling review | ⚠ |
| R9 | Herring — scaling laws in sintering, *J. Appl. Phys.* **21**, 301 (1950) | ⚠ |
| R10 | Faes / Grew et al. — coarsening simulated on FIB-SEM reconstructions | ✂ specifics unknown |
| R11 | NiO/Ni molar volumes ≈ 11.2 / 6.59 cm³ mol⁻¹ | ⚠ |
| R12 | Bond-percolation thresholds SC 0.2488, BCC 0.1803, FCC 0.120 | ⚠ |

**Submission blocker.** R5–R12 must be verified or cut. **B6's measurement does
not depend on R5–R6**: the bracketing in Figure 1 stands on its own, and only
the *framing* as Rayleigh-type break-up requires those citations. If they cannot
be verified, B6 is reported as an area-barrier impossibility without the
Rayleigh attribution, and loses nothing essential.
