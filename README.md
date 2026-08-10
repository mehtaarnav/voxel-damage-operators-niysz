# Connectivity margin of the pristine Ni network as a predictor of retained percolation

Falsification exercise against the Holzer/Pecho Ni-YSZ FIB-tomography dataset
([Zenodo 4056538](https://zenodo.org/records/4056538)).

**Hypothesis under test.** A graph-theoretic "connectivity margin" computed on
the *pristine* Ni network — weighted algebraic connectivity, minimum-cut
conductance, or the lower quantile of neck widths — predicts how much Ni
electronic percolation survives redox degradation, better than Ni volume
fraction, mean particle size, or initial TPB density.

---

## How to re-run

```bash
python phase0_validate_percolation.py     # percolation-code validation
```
```bash
python phase1_download.py && python phase1_read_metadata.py && python phase1_get_papers_pmc.py && python phase1_get_supplementary2.py && python phase1_build_ground_truth.py
```
```bash
python phase2_inspect_labels.py && python phase2_volume_fractions.py
```
```bash
python phase3a_rev_study.py && python phase3_extract_network.py && python phase3_snow_sensitivity.py
```
```bash
python phase4a_validate_tpb.py && python phase4b_tpb.py && python phase4c_metrics.py && python phase4d_particles.py --limit 2 && python phase4e_lambda2_scaling.py
```
```bash
python phase5_percolation.py && python phase6_verdict.py
```

Support library in `cmlib/`; probes and diagnostics are `probe_*.py` and
`diag_*.py`. All figures and tables land in `out/<phase>/`.

`probe_metrics.py` validates the metric implementations against analytically
known networks (series/parallel conductance, bottleneck min-cut, and
λ₂ = 2(1−cos π/n) for path graphs) — run it if you touch `cmlib/metrics.py`.
`phase3_extract_graph.py` is the **failed** skeleton pipeline, kept as the
record of the Phase 3 gate failure; it is not part of the result chain.

---

## Gate results

| Phase | Gate | Result |
|---|---|---|
| 0 | detected site-percolation threshold within ~0.02–0.03 of 0.3116 | **PASS** — 0.31372 raw at *L*=100; 0.31218 after finite-size extrapolation vs 0.311608 reference (deviation 0.0006). 18- and 26-connectivity controls also within 0.0013 of literature |
| 2 | voxel-counted volume fractions within ~10–15 % of published | **PASS** — worst deviation 0.20 % across all 18 (6 samples × 3 phases) |
| 3 (v1, skeleton) | skeleton visually/structurally sound | **FAIL** — see below |
| 3 (v2, watershed) | network visually sound and comparable across anodes | **PASS** |
| 4a | TPB implementation correct | **PASS** — exact on axis-aligned analytic cases; staircase bias measured at 1.713 vs theoretical bound √3 = 1.732 |

### The Phase 3 failure, and what replaced it

`skimage.morphology.skeletonize` produced a **curve** skeleton for the fine
anode but a **medial sheet** for the medium and coarse anodes:

| | fine | medium | coarse |
|---|---|---|---|
| skeleton voxels of degree 2 (curve signature) | 83.1 % | 19.0 % | 4.9 % |
| skeleton voxels of degree ≥ 4 (sheet signature) | 9.2 % | 79.2 % | 94.5 % |
| mean skeleton-voxel degree | 2.58 | 6.69 | 7.68 |
| median branch length | 6.7 vox | 2.0 vox | 2.0 vox |

Reference-shape controls confirm the algorithm is behaving (cylinder → 96.7 %
degree-2; sphere and slab → a point), so this is a property of the
microstructures: thinning preserves topology, and the thick multiply-connected
Ni domains of the coarser anodes have a genuinely 2-dimensional medial axis.

Because the skeleton's *dimensionality* varied monotonically with grain
coarseness, any cross-anode ranking derived from it would have been an artifact
confounded with exactly the variable the hypothesis is tested against. Pruning
does not fix this — a sheet stays a sheet.

**Replacement (agreed with the user):** watershed network extraction (SNOW,
Gostick, *Phys. Rev. E* **96**, 023307 (2017)) via `porespy.networks.snow2`,
applied to the Ni phase. Nodes = Ni chambers, edges = shared interfaces, edge
weights measured directly on the interface. Verified against an analytic shape:
a 6-voxel bar at 10 nm returns `throat.inscribed_diameter` = 60.0 nm and
`throat.cross_sectional_area` = 3600 nm² — both exact.

---

## Stated conventions

Every choice below is a decision, not a default.

**Voxel adjacency (percolation).** 6-connectivity (face-sharing). Required for
the *p*<sub>c</sub> = 0.311608 reference to apply; 18- and 26-connectivity have
thresholds 0.1372 and 0.09755. Also the conservative choice — edge- and
corner-touching voxels share zero interfacial area.

**Spanning.** One connected component touching both index-0 and index-(n−1)
slices along the tested axis; free boundaries.

**Transport axis.** x (array axis 2). The papers plot depth profiles against
"Distance x-coord" over 0–20 µm, matching the x extent of the image windows,
and describe the connectivity check as made "with the inlet plane on the left
(x-direction)".

**Phase labels.** Assigned from the dataset's own stated brightness ordering
(`2_3D_Data_Info.xlsx`, column AE), never from the published numbers — so the
Phase 2 comparison stays an independent test. The three largest-count labels
are the phases; any other label (only 255, ≤ 0.25 % of a volume) is
*unassigned* and stays in the denominator. This matches the papers, whose three
volume fractions for the fine pristine sample sum to 0.997 while its 255 class
occupies exactly 0.249 %.

Note `5_Rx38` genuinely has **Ni and YSZ inverted** relative to the other five
stacks, and every stack uses a different encoding (0/100/200, 1/2/3,
36/121/194, 0/76/150).

**TPB density.** Voxel-edge counting: an interior voxel edge is TPB iff its four
surrounding voxels contain all of Ni, YSZ and pore; contribution is the physical
voxel dimension along that edge, so anisotropy is handled. Computed on the full
stacks by streaming two z-slices at a time — no sub-sampling.

*Known bias, measured not assumed:* the staircase path over-estimates length.
Exact (ratio 1.000) for axis-aligned lines; 1.713 for a worst-case (1,1,1)
line against the theoretical bound √3 = 1.732. The papers used a third
convention again — "the length of each TPB line is determined based on the
skeletonization of TPB-voxels" — so the raw value and the *ratio* to published
are both reported; a constant ratio would indicate a pure convention offset,
which leaves rankings unaffected.

**Watershed markers.** σ = 0.4 voxels Gaussian blur on the distance transform;
`peak_local_max` `min_distance` = 4 voxels (SNOW `r_max` = 4). Both are in
*voxels*, so their physical meaning differs between samples — sensitivity sweeps
are run rather than assuming the default is harmless. Border-truncated regions
are excluded from size statistics.

Sweep outcomes: particle size ordering (fine < medium < coarse) survives
`min_distance` ∈ {2,3,4,6,8}. The network metric rankings survive
`r_max` ∈ {2,4,6,8} for λ₂, min-cut and effective conductance; **neck-width p10
flips medium↔coarse at r_max = 2** (184 vs 179 nm, a 3 % gap), which is reported
in `REPORT.md` rather than smoothed over.

**Edge weight.** Conductance = throat cross-sectional area / throat length, in
nm. The intrinsic conductivity σ is dropped: it multiplies every edge
identically and cancels from every ranking.

**Percolating fraction — three definitions.** The papers' *P* is **not** a
connected-component measure; it "can be obtained from the MIP-PSD analysis",
i.e. simulated intrusion, which measures the fraction *reachable from a
boundary*. So we report `P_span` (touches both faces, our strict definition),
`P_reach` (touches at least one face — the like-for-like analogue of the
published *P*), and `P_largest`.

---

## Known limitations

1. **The coarse anode has no representative volume at this resolution.** Memory
   caps analysis at ~120–150 Mvoxel. At a 10 µm cube the coarse Ni volume
   fraction is 0.2595 against a published 0.229 (+13 %), and its full stack
   gives 0.2293 — the interior is genuinely Ni-rich. Handled by using several
   non-nested ROIs per anode and reporting the spread, but the coarse network
   metrics remain sub-REV and are flagged as such. Volume fraction, TPB density
   and percolation are *not* affected: those are computed on full stacks.

2. **Voxel anisotropy.** The three pristine stacks are 2.4–3.0 % anisotropic
   (handled: `sampling`/`spacing` carry true voxel dimensions; `snow2` takes a
   scalar so the geometric mean is passed, mis-stating lengths by ≤ 1.5 %). The
   post-redox medium and coarse stacks are **40 %** anisotropic
   (17.9 × 17.9 × 25 nm) and are therefore not processed under the scalar
   assumption.

3. **TPB ground truth is digitized from a bar chart.** No table of these values
   exists in either paper or supplement. Self-consistency checked: recomputing
   Figure 7 panel C from panels A/B reproduces the printed percentages to within
   0.7 points.

4. **n = 3.** Only 6 orderings are possible, so any predictor matches by chance
   with probability 1/6 ≈ 17 %. No correlation coefficient or p-value is
   reported.

---

## Ground-truth sources

- **`T-S4`** — Pecho *et al.*, *Materials* **8**(9), 5554–5585 (2015),
  doi:10.3390/ma8095265, **supplementary Table S4**: Φ, *P*, Φ_eff, β, τ,
  M_pred for Ni/YSZ/pore × 3 grains × 2 states. Exact printed values.
- **`P-F7`** — Pecho *et al.*, *Materials* **8**(10), 7129–7147 (2015),
  doi:10.3390/ma8105370, **Figure 7** — digitized bar chart, ±0.05 µm⁻².
- **`T-S1`/`T-S2`** — same supplement: raw-powder PSD, and image dimensions
  (which agree exactly with the Zenodo metadata).

Neither paper's main text contains per-sample numeric tables; the transport
paper's only table is qualitative and the TPB paper has none.
