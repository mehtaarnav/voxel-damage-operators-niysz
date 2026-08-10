# Paper outline — negative result / mechanistic falsification

**Working title**

> Why lower-tail neck engineering did not improve Ni percolation retention:
> a falsification study with a mechanism-specific null under surface-erosion-like
> damage

*Alternatives, if a shorter title is wanted:*
- "Connectivity margin does not beat coarseness, and neck widening has no lever:
  a two-part falsification for Ni-YSZ electrodes."
- "Mechanism-matched microstructure design: a null result for lower-tail neck
  widening under surface-mediated degradation."

**Target framing.** A negative-results / methods paper. The contribution is not
"we failed to find an effect" but **"we falsified a predictor on real data,
built a platform that could have detected the corresponding causal effect, found
none, and identified the mechanistic reason the intervention has no lever."**
Venue candidates: *Journal of Power Sources*, *Journal of the Electrochemical
Society*, or a materials-informatics/negative-results venue. The
pre-registration and the null autopsy are assets to foreground, not caveats to
bury.

---

## Claims (in the order the paper establishes them)

1. **Pristine connectivity-margin metrics do not beat coarseness on the
   Holzer/Pecho dataset.** The metrics that reproduce the retention ordering are
   monotone in coarseness and indistinguishable from mean particle size at n = 3;
   the metrics that are genuinely graph-theoretic (min-cut, effective
   conductance) match no ordering.
2. **Ni volume fraction and pristine TPB density can be anti-predictive of
   retained Ni percolation** — they rank the three anodes exactly backwards —
   and retained TPB is ordered oppositely to retained percolation, so
   "degradation resistance" must name its property.
3. **A real-data-grounded synthetic platform can generate mass-conservative
   lower-tail neck widening**, moving neck p10 by 2× (120 → 240 nm) at fixed p50
   (320 nm) and Φ_Ni conserved to ≤0.02 %.
4. **Under a frozen surface-erosion-like redox surrogate, lower-tail neck
   widening does not change the percolation-loss threshold** — 8.54 vs 8.50 vs
   8.50 rounds, a −0.04-round difference, below one standard error and ~25×
   below the pre-registered threshold.
5. **The null is mechanistically explained: collapse is surface-area-mediated,
   not lower-tail-neck-mediated.** The voxels removed at the failure step are
   one-voxel-deep surface material (mean 20.00/20.02 nm against a background
   median of 28.28 nm), indistinguishable between base and widened structures.
6. **Degradation-mitigation design rules must be mechanism-specific.** An
   intervention aimed at a bottleneck has no lever when the failure mode is
   global thinning.

*Claim boundary to state explicitly in the abstract and again in the discussion:*
claim 4 is a **model-specific null**, not a universal negative about neck
engineering, and the paper must not be readable as the latter.

---

## Section outline

### 1. Introduction
- Ni-YSZ electrode degradation; Ni percolation loss as the failure mode of
  interest; the practical need for a *pristine-state* design metric.
- The intuition under test: percolation is bottleneck-limited, so widening the
  narrow tail of the neck distribution should delay percolation loss.
- Why this intuition is hard to test on real data: coarseness co-varies with
  everything.
- Contribution statement: (i) predictor falsification on real data;
  (ii) a decoupling platform; (iii) a pre-registered causal null; (iv) a
  mechanistic autopsy explaining it.

### 2. Real-data falsification
**2.1 Dataset and outcome.** Three anodes (fine/medium/coarse), pristine and
post-redox, full stacks 0.48–1.11 Gvoxel, no sub-sampling. Retained P_span,
P_reach, published *P*. Establish that **fine-is-worst is the only robust outcome
fact** and **medium vs coarse is unresolved** (definition-dependent flip).
**2.2 Candidate predictors.** Neck-width percentiles, λ₂ (raw and normalised),
face-to-face min-cut, effective conductance; SNOW network extraction; explicit
watershed particle size; Φ_Ni; TPB density.
**2.3 Validation gates.** Percolation threshold recovery (0.31218 vs 0.3116),
volume fractions within 0.20 % of published, TPB implementation exact on
analytic cases, metrics exact against series/parallel and λ₂ = 2(1−cos π/n).
The skeleton-dimensionality gate **failure**, reported and replaced rather than
worked around.
**2.4 Result.** The coarseness-proxy / failed-graph-metric split. λ₂ ∝ N^(−0.867)
and the λ₂·N collapse (7.65× → 1.24×). Neck p10 is the only well-separated
metric and reproduces the particle-size ordering exactly. Marker-parameter
sensitivity.
**2.5 The anti-predictive result.** Φ_Ni and pristine TPB rank backwards;
retained TPB runs opposite to retained percolation.
**2.6 Statistical standing.** n = 3 is directional; no p-values, by design;
1/6 and 2/6 chance baselines stated.

### 3. Synthetic decoupling platform
**3.1 Rationale.** The confound argument: coarseness varies monotonically across
the real anodes, so a predictor and a coarseness proxy are not separable at
n = 3. Decoupling requires structures matched on loading and particle size that
differ only in the neck tail.
**3.2 Generator.** Jittered lattice of spheres with explicit cylindrical necks;
161 × 168 × 168 vox @ 20 nm; 96 particles, 224 bonds; R = 12.1 vox; bimodal neck
mixture; Φ_Ni = 0.250 and YSZ/pore split taken from the medium anode.
**3.3 Intervention.** Max-clip of the lower tail plus a radius solve for mass
conservation. Achieved-ratio rungs (1.00 / 1.33 / 2.00×) rather than nominal
targets, and why (discrete value grid).
**3.4 Qualification.** P2-A composition/topology PASS (and the finding that plain
nearest-neighbour adjacency already meets the coordination band — no topology
machinery was needed); P2-B neck distribution PASS with the p50/p10 = 2.90 note;
P2-C split into C1 PASS / C2 DEFERRED.
**3.5 Known confound and what it forbids.** The 5.17–6.79 % radius shrink,
pre-registered before any outcome, with its attribution rules.
**3.6 Damage model.** D4 definition; frozen parameters and their honest
provenance (E0 spike, not re-derived); re-validation on the new geometry;
bisection rather than a fixed intensity grid, vindicated by the transition
falling inside E0's unresolved gap.

### 4. Causal experiment
**4.1 Pre-registered design.** 5 structure seeds × 3 groups × 5 damage seeds =
75 bisections; primary outcome = transition midpoint; 1.0-round interpretability
threshold; damage-seed averaging rule frozen before any widened structure was
inspected.
**4.2 Result.** 8.54 / 8.50 / 8.50; 74 of 75 midpoints exactly 8.5; zero variance
in every widened structure; effect below one standard error and negative in sign.
**4.3 Integrity checks.** (i) Damage-seed variance did not reproduce across seed
sets (sd 0.594 → 0.200) — reported against interest, and it sharpens the null.
(ii) Provenance verification of the widened masks: 67,803 voxels differ, 46/224
necks changed, balanced add/remove counts. **The null is not a pipeline
artifact.**
**4.4 Resolution ceiling.** What the design could not have seen.

### 5. Null autopsy
**5.1 Question.** Does a lower-tail advantage survive to the collapse boundary?
**5.2 D1 remaining volume** — converged within 0.5 pp, group differences below
within-group sd, widened slightly *worse*.
**5.3 D2 — a degenerate proxy, reported as a methodological warning.** An EDT
percentile over all phase voxels measures surface fraction, not neck thickness:
p10 pins to one voxel for every group **including at n = 0** where the structures
demonstrably differ. Present this as a transferable pitfall for the field, and
state plainly that it contributes no evidence in either direction.
**5.4 D3 backbone size** — P_span = 1.000 everywhere at n = 8; nothing to detect;
base retains 0.04 at n = 9 while both widened groups sit at 0.00.
**5.5 D4loc — the decisive diagnostic.** Removed-voxel EDT p10 = p50 = p90 =
20.0 nm against background median 28.28 / p90 ≈ 73 nm; indistinguishable between
base and 2.00×; removed fraction 0.17 vs 0.18.
**5.6 Interpretation.** Bypassed, not erased. Branch A in the pre-registered
wording, with the surface-mediation refinement. **Branch C (hidden continuous
benefit) explicitly not supported** — all continuous diagnostics point marginally
the wrong way.

### 6. Discussion
- The two halves are one finding: a pristine geometric margin can only predict
  retention if the degradation operator is sensitive to that margin.
- Mechanism-matched design: what surface-mediated thinning does and does not
  expose as a lever.
- Why negative results of this shape are worth publishing — the platform, the
  pre-registration, and the autopsy are reusable independently of the null.
- Scope and non-claims, restated: this is a null for one intervention under one
  frozen operator, not a universal negative.

### 7. Limitations
Transcribe §F of the memo in full: phenomenological damage model with
non-derived frozen parameters; one-round outcome resolution; synthetic TPB
7.3–19.2× real; deferred image-based size comparability; recorded (moot)
radius-shrink confound; invalid D2; no real redox calibration; no Family C;
real-data limits (sub-REV coarse ROIs, unresolved medium-vs-coarse, digitized
TPB ground truth, anisotropic post-redox voxels, conductance-weighting choice);
the porespy bug and its impact assessment; n = 3.

### 8. Future work
Neck-selective damage models (with the D1-tautology caveat); hierarchical
backbone architectures; real-data calibration of damage and size metrics; TPB
realism; targeted bottleneck interventions; a finer-grained outcome variable.

### Appendices
- A. Pre-registration, with all amendments and their dates.
- B. Validation-gate record, including the two reported gate failures.
- C. Reproducibility manifest (scripts, data, frozen parameters).
- D. The porespy `snow2` chunked-watershed bug and impact assessment.

---

## Figures and tables

| # | item | type | source |
|---|---|---|---|
| T1 | Real-data predictor/outcome table — all pristine predictors vs all retention outcomes, three anodes | table | `out/phase6/phase6_comparison_table.csv` |
| F1 | λ₂ scaling collapse — λ₂ vs N over 21 ROIs with the N^(−0.867) fit, inset showing λ₂·N collapsing 7.65× → 1.24× | scatter + inset bar | `phase4e_lambda2_scaling.py`, `out/phase4/phase4e_lambda2_scaling.png` |
| T2 | Platform v2 generator qualification — per-seed Φ_Ni, degree, achieved p10, neck p10/p50, radius shrink, node counts, P_span | table | `out/platform_v2/qualification_run.csv` |
| F2 | p10-group transition-midpoint plot — midpoint vs achieved p10 ratio, individual damage seeds as points, group means with SE, 1.0-round threshold band drawn | strip/dot plot | `out/platform_v2/p10_group_experiment.csv` |
| T3 | Autopsy remaining-volume table — D1 and D3 combined, per group at n = 0 / 8 / 9 | table | `out/platform_v2/null_autopsy.csv` |
| F3 | D4loc removed-voxel thickness distribution — histogram of removed-voxel EDT at the transition step, base and 2.00× overlaid, with the *all-voxel* n = 8 background distribution behind | overlaid histogram | `out/platform_v2/null_autopsy_localization.csv` |
| F4 | Schematic: intervention bypassed by surface erosion — a widened neck alongside a base neck, with erosion fronts advancing uniformly from all surfaces, showing both networks reaching sub-percolation at the same cumulative stripping | schematic | new artwork |

**Supporting/supplementary figures**
- S1 Example slices, pristine and post-redox, all three real anodes
  (`out/phase2/phase2_slices_*.png`).
- S2 Percolation-threshold validation and finite-size scaling
  (`out/phase0/phase0_finite_size_scaling.png`).
- S3 SNOW marker-parameter (`r_max`) sensitivity of the metric orderings
  (`out/phase3/phase3c_rmax_sensitivity.csv`).
- S4 The skeleton-dimensionality gate failure
  (`out/phase3/phase3_GATE_FAILURE_skeleton_dimensionality.png`).
- S5 Neck-width distributions, base vs 1.33× vs 2.00×, showing the tail-selective
  shift at fixed p50.
- S6 D2 degeneracy figure — EDT percentile identical across groups at n = 0,
  with a one-voxel reference line; the methodological-warning figure.

**Figure-count discipline.** F1–F4 plus T1–T3 carry all six claims. If the venue
caps figures, F4 (schematic) and F3 (D4loc) are the two that must survive — they
are the mechanism.

---

## Framing notes for drafting

- **Lead with the mechanism, not the absence.** The abstract's last sentence
  should be claim 5, not claim 4.
- **State the claim boundary in the abstract.** One clause: "under this
  surface-erosion-like operator."
- **Report the two against-interest findings prominently** (the non-reproducing
  damage-seed variance, and the degenerate D2). They are the paper's credibility.
- **Do not compute statistics the data cannot bear.** No p-values on n = 3; the
  synthetic comparison is reported as an effect size against a pre-registered
  threshold and a standard error, not as a significance test.
