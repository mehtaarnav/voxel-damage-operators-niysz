# Verdict: does a pristine-state connectivity margin predict retained Ni percolation?

**Short answer: no — not in the sense the hypothesis requires.**

The connectivity-margin metrics split into two groups. The ones that reproduce
the outcome ordering (10th-percentile neck width, algebraic connectivity)
produce *exactly the same ordering as mean Ni particle size*, so at n = 3 they
are not distinguishable from a plain coarseness measure — which is the specific
comparison the hypothesis was supposed to win. The ones that are genuinely
graph-theoretic rather than a proxy for coarseness (minimum-cut conductance,
effective network conductance) **fail**: they match no outcome ordering at all.

Separately, and independently of the hypothesis, the exercise produced a clean
result worth keeping: **Ni volume fraction and initial TPB density rank the
three anodes exactly backwards** for retained Ni percolation.

---

## 1. What the outcome actually is

Retained Ni percolation, measured on the full stacks (0.48–1.11 gigavoxel each,
no sub-sampling), degraded ÷ pristine:

| anode | P_span retained | P_reach retained | published *P* retained |
|---|---|---|---|
| fine   | **0.680** | **0.857** | **0.821** |
| medium | 0.855 | 0.979 | 0.916 |
| coarse | 0.947 | 0.942 | 0.924 |

Our absolute values track the published *P* closely: `P_reach` matches to within
0.1–7.4 % on all six stacks (ratios 1.001, 1.004, 1.008, 1.028, 1.045, 1.074).

**The one robust fact in this table is that the fine anode retains Ni
percolation worst.** The papers say so independently: Ni *P* "drops from 0.99 to
0.80" in the fine sample "compared to 0.90 in the coarse sample"
(ma8095265, §3). Medium vs coarse is a near-tie that *flips with the
definition* — `P_span` puts coarse first (0.947 vs 0.855), `P_reach` puts medium
first (0.979 vs 0.942), and the published values are 0.924 vs 0.916, a 0.8-point
gap. **Medium vs coarse should be treated as unresolved.**

This matters a great deal for how much the "matches" below are worth: if the
only reliable outcome fact is *fine is last*, then a predictor only has to put
fine last to score a match, and there are 2 such orderings out of 6.

### A second outcome, running the opposite way

Retained TPB density is **not** ordered like retained percolation:

| anode | TPB retained (this work) | TPB retained (published) |
|---|---|---|
| fine   | **0.799** | **0.743** |
| medium | 0.746 | 0.586 |
| coarse | 0.590 | 0.608 |

Fine is *best* at retaining TPB and *worst* at retaining Ni percolation. Any
claim about "degradation resistance" therefore has to name which property it
means. The hypothesis names Ni electronic percolation, so that is what is scored
below.

---

## 2. The comparison table

Pristine-state predictors, mean ± sd over non-nested 8 µm ROIs
(6 / 6 / 9 ROIs for fine / medium / coarse):

| quantity | fine | medium | coarse |
|---|---|---|---|
| **λ₂ raw** (weighted Laplacian) | 0.458 ± 0.334 | 1.228 ± 0.837 | 3.504 ± 5.905 |
| **λ₂ normalised** | 0.00094 ± 0.00049 | 0.00210 ± 0.00187 | 0.00230 ± 0.00362 |
| **min-cut, face-to-face** | 1205 ± 630 | 1323 ± 1128 | 283 ± 172 |
| **effective conductance** | 254 ± 92 | 334 ± 241 | 121 ± 69 |
| **neck width p10** (nm) | 69.5 ± 4.6 | 156.9 ± 43.0 | 205.3 ± 45.0 |
| neck width p50 (nm) | 297.9 ± 25.1 | 516.8 ± 96.8 | 631.3 ± 59.2 |
| network nodes per ROI | 509 ± 50 | 165 ± 20 | 82 ± 24 |
| mean chamber diameter, SNOW (nm) | 671 ± 24 | 901 ± 77 | 1082 ± 46 |
| **mean Ni particle size, explicit watershed** (nm, vol-weighted) | 1148 | 1445 | 1715 |
| **Ni volume fraction** | 0.3221 | 0.2497 | 0.2293 |
| *published Φ_Ni* | *0.322* | *0.250* | *0.229* |
| **pristine TPB density** (µm⁻²) | 3.624 | 2.109 | 1.473 |
| *published TPB_total* | *2.65* | *2.03* | *1.07* |

Degraded-state outcomes:

| quantity | fine | medium | coarse |
|---|---|---|---|
| still percolates after redox | yes | yes | yes |
| P_span pristine → degraded | 0.984 → 0.669 | 0.958 → 0.819 | 0.932 → 0.882 |
| P_reach pristine → degraded | 0.987 → 0.846 | 0.969 → 0.949 | 0.967 → 0.911 |
| *published P pristine → degraded* | *0.985 → 0.809* | *0.965 → 0.884* | *0.959 → 0.886* |
| TPB density pristine → degraded (µm⁻²) | 3.624 → 2.896 | 2.109 → 1.573 | 1.473 → 0.868 |

Full machine-readable version: `out/phase6/phase6_comparison_table.csv`.

---

## 3. The three questions, answered

### Q1. Do the connectivity-margin metrics rank the anodes best-to-worst in the same order as the measured retained percolation?

**Two of the five do; three do not.**

| predictor | ordering (best first) | matches retained Ni percolation? |
|---|---|---|
| 10th-percentile neck width | coarse > medium > fine | **yes** |
| λ₂ raw | coarse > medium > fine | **yes** |
| λ₂ normalised | coarse > medium > fine | **yes** |
| min-cut face-to-face | medium > fine > coarse | **no** |
| effective conductance | medium > fine > coarse | **no** |

against the outcome ordering `coarse > medium > fine` (from `P_span` retention
and from the published *P* retention, which agree).

Three qualifications that materially weaken the two "yes" results:

1. **λ₂ raw is not measuring connectivity here — it is measuring node count.**
   This is measured, not asserted (`phase4e_lambda2_scaling.py`). Fitting all
   21 ROIs pooled gives λ₂ ∝ N^(−0.867), close to a pure 1/N law. And the
   diagnostic that matters:

   | | fine | medium | coarse | spread |
   |---|---|---|---|---|
   | λ₂ raw | 0.458 | 1.228 | 3.504 | **7.65×** |
   | λ₂ × N | 232 | 200 | 247 | **1.24×** |

   Multiplying by node count collapses a 7.65-fold between-anode difference to
   1.24-fold. Almost everything λ₂_raw appeared to say about the three anodes
   was a restatement of how many watershed chambers each contains — i.e. of
   their coarseness. The normalised Laplacian was computed precisely to remove
   this, and once it is removed the separation nearly vanishes:
   0.00094 / 0.00210 / 0.00230, with standard deviations
   (0.00049 / 0.00187 / 0.00362) that overlap all three anodes.
2. **The coarse λ₂ scatter exceeds its own mean** (3.50 ± 5.90). One ROI gave
   19.1, another 0.32.
3. **The match is to `P_span`, not `P_reach`.** No predictor at all matches the
   `P_reach` ordering (`medium > coarse > fine`). Since medium and coarse are
   tied to within noise, the surviving content of the match is only "fine is
   last".

Only **neck-width p10 separates the anodes with non-overlapping error bars**
(69.5 ± 4.6 vs 156.9 ± 43.0 vs 205.3 ± 45.0). It is the sole connectivity metric
here with a defensible signal.

#### Robustness to the watershed marker parameter

Sweeping SNOW's `r_max` over {2, 4, 6, 8} voxels
(`phase3c_rmax_sensitivity.csv`):

| metric | ordering at r_max = 2 / 4 / 6 / 8 |
|---|---|
| λ₂ raw | coarse > medium > fine — **stable** |
| min-cut | medium > fine > coarse — **stable** (and consistently wrong) |
| effective conductance | medium > fine > coarse — **stable** (consistently wrong) |
| neck width p10 | *medium > coarse* > fine at r_max=2; coarse > medium > fine at 4, 6, 8 |

So the ranking of the only well-separated connectivity metric is **not fully
stable** to the marker parameter: at the finest segmentation medium and coarse
swap. The swap is small (184 vs 179 nm — a 3 % gap) and it happens in exactly
the pair that the *outcome* also fails to resolve. **Every metric puts fine last
at every parameter setting except min-cut and g_eff, which put coarse last at
every setting.** That is the honest summary: the resolvable content on both
sides is "fine is worst", and medium vs coarse is noise in predictor and outcome
alike.

Against retained **TPB** the answer reverses: retained TPB is ordered
`fine > medium > coarse`, which none of the connectivity-margin metrics match,
and which Ni volume fraction and pristine TPB density both match exactly.

### Q2. Does that ranking differ from what volume fraction or mean particle size alone would predict?

**Versus volume fraction: yes — oppositely.** Ni volume fraction orders the
anodes `fine (0.322) > medium (0.250) > coarse (0.229)`, i.e. it predicts the
fine anode is the most robust. The fine anode is in fact the *least* robust.
Volume fraction is not merely uninformative here, it is anti-correlated. The
same is true of pristine TPB density.

**Versus mean particle size: no — and this is what sinks the hypothesis.**
Mean Ni particle size, computed the way you specified (`scipy.ndimage`
distance transform → `skimage.feature.peak_local_max` markers →
`skimage.segmentation.watershed`, volume-weighted equivalent-sphere diameter),
orders the anodes `coarse (1715 nm) > medium (1445) > fine (1148)`. The SNOW
chamber diameters agree: `coarse (1082) > medium (901) > fine (671)`.

That is **the identical ordering** produced by neck-width p10 and by both forms
of λ₂. Neck width and particle size are both monotone in coarseness, so on three
samples they are indistinguishable by construction — there is no ordering either
could produce that would separate them.

This is not an artifact of the watershed marker parameter: the ordering
`fine < medium < coarse` survives every `min_distance` in {2, 3, 4, 6, 8}
(volume-weighted diameters 1024→1232, 1318→1574, 1471→1870 nm across the sweep).
So the particle-size ranking is robust, and the connectivity metrics reproduce
it rather than beating it.

The hypothesis's claim is specifically that a connectivity margin predicts
retention **better than** mean particle size. This dataset cannot support that
claim, because the connectivity metric that works produces exactly the particle
size ranking, and the connectivity metrics that are *not* reducible to
coarseness — min-cut and effective conductance, which actually weight paths by
conductance and account for network topology — are the ones that **fail**.

### Q3. Statistical power

**n = 3 is a directional check, not a statistically powered result.** There are
3! = 6 possible orderings, so any predictor reproduces the outcome ordering with
probability 1/6 ≈ 17 % by chance; and because medium vs coarse is unresolved,
the effective bar is "put fine last", which 2 of 6 orderings clear by chance —
33 %. **No correlation coefficient and no p-value is reported anywhere in this
work**, and none should be computed from these three points. The Phase 6 script
enforces this.

---

## 4. What would actually test the hypothesis

The failure mode here is a confound, not a null result: in this dataset
coarseness varies monotonically across the three samples and drags every
candidate metric with it. To separate a connectivity margin from a coarseness
proxy you need samples where they **disagree** — e.g. two anodes matched on mean
particle size and Ni volume fraction but differing in neck-width distribution
(bimodal vs narrow), or a designed/simulated microstructure series where neck
statistics are varied at fixed particle size. Three samples that differ in
"everything at once" cannot do it, regardless of how carefully the graph metrics
are computed.

It is also worth noting that min-cut and effective conductance failing is *not*
obviously a defect of those metrics: both are strongly affected by the coarse
anode's sub-REV ROIs (see limitation 1), where an 8 µm cube frequently contains
only 2–15 network nodes touching a face, and one ROI contained no spanning Ni
cluster at all. A larger domain might rehabilitate them. That is a real open
question this run cannot settle.

---

## 5. Gate record

| phase | gate | result |
|---|---|---|
| 0 | percolation threshold within ~0.02–0.03 of 0.3116 | **PASS** — 0.31218 after finite-size extrapolation (deviation 0.0006); 18- and 26-connectivity controls within 0.0013 |
| 2 | volume fractions within ~10–15 % of published | **PASS** — worst deviation 0.20 % over 18 values |
| 3 (skeleton) | skeleton visually/structurally sound | **FAIL — stopped and reported** |
| 3 (watershed) | network sound and comparable across anodes | **PASS** |
| 4a | TPB implementation correct | **PASS** — exact on analytic cases; staircase bias measured 1.713 vs bound √3 |
| 4 | TPB within factor 2 of published, ordering preserved | **PASS** — ratios 1.04–1.47, ordering fine > medium > coarse preserved |
| — | metrics implementation | **PASS** — exact against series/parallel conductance, bottleneck min-cut, and λ₂ = 2(1−cos π/n) |

The Phase 3 skeleton failure was reported rather than worked around:
`skeletonize` returned a curve skeleton for the fine anode (83 % of skeleton
voxels degree-2) but a medial *sheet* for medium and coarse (79 % and 94 %
degree ≥ 4, spiking at degree exactly 8). Because the skeleton's dimensionality
varied monotonically with coarseness, every downstream metric would have been
confounded with the variable under test. Replaced, with your approval, by SNOW
watershed network extraction.

---

## 6. Limitations that bound these conclusions

1. **The coarse anode is sub-REV for the network metrics.** Memory caps ROI
   analysis at ~120–150 Mvoxel. At a 10 µm cube the coarse Ni volume fraction is
   0.2595 against a full-stack 0.2293 (+13 %). Handled by 9 non-nested ROIs with
   reported spread, but the coarse min-cut and λ₂ values remain unreliable — one
   ROI had no spanning cluster. Volume fraction, TPB density and percolation are
   **not** affected: those were computed on full stacks.
2. **Medium vs coarse retention is within noise** (published 0.916 vs 0.924).
   Only "fine is worst" is robust.
3. **TPB ground truth is digitized from a bar chart** (Figure 7); no table
   exists. Self-consistency verified to within 0.7 percentage points against the
   paper's own panel C.
4. **Post-redox medium and coarse stacks are 40 % voxel-anisotropic**
   (17.9 × 17.9 × 25 nm). Percolation (topological) and TPB (true per-axis
   lengths) handle this correctly; the scalar-voxel SNOW extraction was not
   applied to them.
5. **The conductance weighting is a modelling choice**, not a measurement:
   cond = throat area / throat length, with intrinsic σ dropped. `neck_nm` is
   also reported unweighted so the conclusion can be checked against it.

---

## Addendum (2026-08-10, filed during follow-on work)

A `porespy.networks.snow2` default-parameter bug (silently-enabled chunked
watershed partitioning) was found and fixed for new work
(`cmlib/pnm.py`, now defaults to serial extraction). It was silently active
throughout the SNOW extractions above (Phase 3/4). **No conclusion in this
report changes** — full impact assessment, the one ROI directly re-checked,
and the reasoning bound for the other 20, are in
[`IMPACT_NOTE_porespy_parallel_bug.md`](IMPACT_NOTE_porespy_parallel_bug.md).
Nothing above this addendum was altered.
