# Reproducibility manifest

Everything required to regenerate the results in `PATH_B_MEMO.md` from a clean
checkout. Repository: `connectivity_margin`, branch `master`, head commit
`b21c2bf` ("null autopsy: intervention bypassed, not erased — D4 collapse is
surface-mediated").

---

## 0. Environment and inputs

**Environment check:** `check_env.py` — reports presence and version of every
required package: `numpy`, `scipy`, `skimage`, `networkx`, `tifffile`, `skan`,
`matplotlib`, `pandas`, `openpyxl`, `requests`. Also uses `porespy` (SNOW
network extraction).

**Not committed, regenerable (see `.gitignore`):**

| path | what | how to regenerate |
|---|---|---|
| `data/` | Holzer/Pecho segmented tomography, ~2.2 GB (6 zipped stacks + metadata) | `phase1_download.py`, then unzip |
| `refs/` | papers and supplementary | `phase1_get_papers_pmc.py`, `phase1_get_supplementary2.py` |
| `out/networks/`, `out/graphs/` | large binary intermediates | fully determined by code + seeds |

**Real dataset.** Six segmented stacks (Rx36, Rx37, Rx38, Rx41-1/2/3) =
fine/medium/coarse × pristine/post-redox, 0.48–1.11 Gvoxel each. Papers:
`ma8095265` (transport) and `ma8105370` (TPB).

**Known environment caveat.** A `porespy.networks.snow2` default-parameter bug
(silently-enabled chunked watershed partitioning) was active during the Phase 3/4
SNOW extractions. Fixed for new work in `cmlib/pnm.py` (now defaults to serial
extraction). Impact assessed in `IMPACT_NOTE_porespy_parallel_bug.md`; **no
conclusion changes**, and nothing in Phase 3/4 was retroactively altered.

---

## 1. Scripts required to regenerate the main results

### 1.1 Library (`cmlib/`)

| module | role |
|---|---|
| `io.py` | stack loading |
| `phases.py` | phase labelling, volume fractions |
| `roi.py` | non-nested ROI selection |
| `percolation.py` | `P_span`, `P_reach`, spanning-cluster extraction |
| `graph.py` | Laplacian / λ₂, min-cut, effective conductance |
| `pnm.py` | SNOW network extraction (**serial by default — bug fix**) |
| `particles.py` | explicit watershed particle sizing |
| `tpb.py` | TPB density |
| `metrics.py` | metric assembly |
| `ground_truth.py` | published-value tables |
| `synth.py` | Platform v2 generator: lattice geometry, neck mixtures, max-clip widening, mass-conservative radius solve, rasterization |
| `synthvol.py` | ternary volume container and I/O |
| `damage.py` | YSZ/pore placement, **D4 operator**, ternary rebuild |
| `api.py` | shared entry points |

### 1.2 Part A — real-data falsification (run in order)

| # | script | produces |
|---|---|---|
| 1 | `phase1_download.py` → `phase1_read_metadata.py` → `phase1_extract_tables.py` → `phase1_get_papers_pmc.py` → `phase1_get_supplementary2.py` → `phase1_find_definitions.py` → `phase1_build_ground_truth.py` | `out/phase1/*` |
| 2 | `phase0_validate_percolation.py` | `out/phase0/*` — threshold-recovery gate |
| 3 | `phase2_inspect_labels.py`, `phase2_volume_fractions.py` | `out/phase2/*` |
| 4 | `phase3_extract_network.py`, `phase3_snow_sensitivity.py`, `phase3a_rev_study.py` | `out/phase3/*` |
| 5 | `phase4a_validate_tpb.py`, `phase4b_tpb.py`, `phase4c_metrics.py`, `phase4d_particles.py`, `phase4e_lambda2_scaling.py` | `out/phase4/*` |
| 6 | `phase5_percolation.py` | `out/phase5/*` — the outcome variable |
| 7 | `phase6_verdict.py` | `out/phase6/*` — comparison table and rankings |

*Superseded / diagnostic, retained for the record:* `phase3_extract_graph.py`
(skeleton route, **gate failure — replaced by SNOW**), `diag_skeleton.py`,
`diag_skeleton_figure.py`, and the `probe_*.py` family (`probe_memory`,
`probe_skan`, `probe_skan2`, `probe_porespy`, `probe_metrics`,
`probe_local_thickness`, `probe_snow2_parallel_bug`).
*Alternate paper-fetch paths not on the main route:* `phase1_get_papers.py`,
`phase1_get_supplementary.py`.

### 1.3 Part B — synthetic platform (run in order)

| # | script | produces |
|---|---|---|
| 1 | `scripts/next/phase0_validate_synthetic_pipeline.py` | synthetic-pipeline validation |
| 2 | `scripts/next/t5_coupling_experiment.py`, `t5b_coupling_experiment.py`, `t5b_reanalyze_strict_gating.py` | `out/next/t5*` — coupling checks |
| 3 | `scripts/next/familyB_pilot.py` | `out/next/familyB_pilot*` |
| 4 | `scripts/spike/e0_vertical_slice.py`, `e0b_saturation_bridge.py` | `out/spike/*` — E0 spike; **origin of `p_erode = 0.35`, `expand_vox = 1`** |
| 5 | `scripts/platform_v2/design_probe.py` | `out/platform_v2/design_probe_log.txt` |
| 6 | `scripts/platform_v2/qualification_run.py` | `out/platform_v2/qualification_*` — gates P2-A/B/C |
| 7 | `scripts/platform_v2/cpsd_candidate_precheck.py`, `cpsd_convergence_and_recompute.py`, `real_cpsd_variability.py` | `out/platform_v2/cpsd_*`, `real_cpsd_variability*` — the size-metric work that led to the P2-C2 deferral |
| 8 | `scripts/platform_v2/ternary_and_d4_pilot.py` | `ternary_placement.csv`, `d4_revalidation.csv`, `d4_bisection_base.csv` |
| 9 | `scripts/platform_v2/p10_group_experiment.py` | **`p10_group_experiment.csv` — the primary result** (~657 s) |
| 10 | `scripts/platform_v2/null_autopsy.py` | `null_autopsy.csv`, `null_autopsy_localization.csv` (~441 s) |

---

## 2. Committed data and reports

### 2.1 Reports (prose, with numbers)

| path | content |
|---|---|
| `README.md` | project overview |
| `REPORT.md` | **Part A verdict** — predictor falsification, gate record, limitations, porespy addendum |
| `IMPACT_NOTE_porespy_parallel_bug.md` | bug description and impact assessment |
| `out/next/preregistration.md` | **pre-registration + all amendments** |
| `out/next/EXECUTION_SPEC.md` | execution spec for the synthetic phase |
| `out/next/phase_minus1_prior_art.md` | prior-art review |
| `out/next/familyB_pilot_report.md` | Family B pilot |
| `out/next/t5_coupling_decision_report.md`, `t5b_coupling_decision_report.md` | coupling decisions |
| `out/spike/e0_vertical_slice_report.md` | E0 spike |
| `out/platform_v2/design_memo.md` | Platform v2 design rationale |
| `out/platform_v2/qualification_report.md` | gates P2-A/B/C |
| `out/platform_v2/ternary_d4_bisection_report.md` | YSZ/pore placement, D4 re-validation, base-only bisection |
| `out/platform_v2/p10_group_report.md` | **primary result + provenance verification** |
| `out/platform_v2/null_autopsy_report.md` | **null autopsy** |
| `out/writeup/PATH_B_MEMO.md` | final memo |
| `out/writeup/PAPER_OUTLINE.md` | paper outline |
| `out/writeup/REPRODUCIBILITY_MANIFEST.md` | this file |

### 2.2 CSVs backing each claim

**Part A**

| claim | file |
|---|---|
| percolation-threshold gate | `out/phase0/phase0_sweep_{coarse,fine,conn}.csv` |
| published ground truth, TPB digitization check | `out/phase1/phase1_ground_truth.csv`, `phase1_tpb_digitization_check.csv` |
| volume fractions, phase profiles | `out/phase2/phase2_volume_fractions.csv`, `phase2_profile_*.csv`, `phase2_label_inventory.csv` |
| SNOW networks, marker sensitivity, REV | `out/phase3/phase3_graphs_8.0um.csv`, `phase3_snow_8.0um_rmax4.csv`, `phase3c_rmax_sensitivity.csv`, `phase3a_rev.csv`, `diag_skeleton.csv` |
| TPB, metrics, particle size, λ₂ scaling | `out/phase4/phase4b_tpb_full_stacks.csv`, `phase4c_metrics_per_{roi,anode}_8.0um.csv`, `phase4d_particles.csv` |
| **outcome — retained percolation** | `out/phase5/phase5_percolation.csv`, **`phase5_retention.csv`** |
| **comparison table and rankings** | **`out/phase6/phase6_comparison_table.csv`**, `phase6_rankings.csv` |

**Part B**

| claim | file |
|---|---|
| coupling experiments | `out/next/t5_coupling_experiment{,_agg}.csv`, `t5b_coupling_experiment.csv`, `t5b_deviations{,_agg}.csv`, `t5b_strict_gating.csv` |
| Family B pilot | `out/next/familyB_pilot.csv`, `familyB_pilot_{gating,deviations,base_validity_log}.csv` |
| E0 spike, D4 parameter origin | `out/spike/e0_combined_all_intensities.csv`, `e0_vertical_slice.csv`, `e0b_saturation_bridge.csv` |
| **generator qualification** | **`out/platform_v2/qualification_run.csv`**, `qualification_deviations.csv`, `qualification_base_validity_log.csv`, `qualification_gating_log.txt`, `qualification_analysis_log.txt` |
| size-metric deferral evidence | `cpsd_candidate_precheck.csv`, `cpsd_convergence.csv`, `cpsd_convergence_log.txt`, `real_cpsd_variability.csv` |
| YSZ/pore placement, D4 re-validation, base bisection | `ternary_placement.csv`, `d4_revalidation.csv`, `d4_bisection_base.csv` |
| **primary result** | **`out/platform_v2/p10_group_experiment.csv`** |
| **autopsy** | **`out/platform_v2/null_autopsy.csv`**, **`null_autopsy_localization.csv`**, `null_autopsy_log.txt` |

---

## 3. Pre-registration sections and amendments

All in `out/next/preregistration.md`. Ordering below is chronological by freeze
date, not by section number — the `§0x` sections were appended as amendments and
some appear after later-numbered sections in the file.

| § | line | title | status |
|---|---|---|---|
| 1 | 345 | Claim boundaries — what this study can and cannot conclude | original |
| 2 | 374 | Primary and secondary outcome metrics | original |
| 3 | 386 | **Damage-model hierarchy D4/D2/D1, frozen before Phase 2 was built** | original |
| 4 | 407 | How the real Holzer/Pecho data is used (qualitative anchor, **not a fit target**) | original |
| 5 | 420 | Damage-intensity selection, frozen before Phase 3 results were seen | original |
| 6 | 442 | Noise floor | original |
| 7 | 448 | Positive-effect criteria for gate G3 | original |
| 8 | 473 | What would make the design untestable (gate G1) | original |
| 9 | 490 | Sign-off | original |
| 0 | 19 | Refined primary causal question and measurement hierarchy | frozen after T5 |
| 0a | 44 | T5b — corrected experiment required before Family B | amendment |
| 0b | 63 | Acceptance criteria T5b → Family B (supersedes §7–8 at that checkpoint) | amendment |
| 0c | 96 | Amendments after T5b review | frozen 2026-08-10 |
| 0d | 157 | Family B disordered pilot | amendment |
| 0e | 201 | Amendments after E0 spike review | frozen 2026-08-10 |
| 0f | 302 | Platform v2 Ni generator qualification — RESULT | 2026-08-10 |
| 0f | 498 | **Amendments after qualification review**: P2-C split (C1 PASS / C2 DEFERRED); `generator_radius_deviation` qualified; radius guardrail; **causal-interpretation obligation for the radius-shrink confound**; "fixed particle size" language banned | frozen 2026-08-10 |
| 0g/1 | 564 | **Damage-seed averaging requirement** — ≥3 seeds, 5 if cheap; no single-seed comparison; sub-1.0-round differences must be completed to 5 seeds before interpretation | frozen 2026-08-10, **before any widened structure was inspected under damage** |
| 0g/2 | 588 | **Status of `p_erode = 0.35`, `expand_vox = 1`** — E0 origin, not re-derived, deliberately frozen; changing either defines a new model and needs an amendment; sensitivity check obligated only on a positive branch | frozen 2026-08-10 |
| 0g/3 | 612 | **TPB magnitude caveat, corrected** — 7.3–19.2× real, superseding the earlier understated "8–10×" | frozen 2026-08-10 |

**Rules that actually bound the reported result:**
- §3's D1 rule: an effect under D1 but not D4 and not D2 scores **NEGATIVE**.
- §0g/1: the primary comparison ran at the maximum 5 damage seeds, so the
  seed-completion remedy is exhausted — no further seeds could rescue resolution.
- §0g/2: the sensitivity obligation was **not triggered** (it applies to a
  positive/weak-positive branch).
- §0f/4: the matched-shrink control was **not triggered** (no effect to attribute).
- §0f/1: P2-C2 remains an open, named obligation — deferred, not dropped.

---

## 4. Exact frozen parameters

### 4.1 Generator (Platform v2) — `scripts/platform_v2/p10_group_experiment.py:41-49`

| parameter | value |
|---|---|
| `R_VOX` (sphere radius) | 12.1 vox = 242 nm |
| `PITCH` | 32 vox |
| `NLAT_Z`, `NLAT_XY` | 6, 4 → 96 particles, 224 nearest-neighbour bonds |
| `MARGIN`, `JITTER` | 8 vox, 0.15 · pitch |
| domain shape | 161 × 168 × 168 vox = 3.22 × 3.36 × 3.36 µm |
| `VOXEL_NM` | 20.0 |
| `GEOM_SEED` | 999 |
| `FRAC_WEAK`, `WEAK_RANGE`, `NORMAL_RANGE` | 0.20, (4, 6) vox, (12, 20) vox |
| `MIN_RATIO` (p50/p10 validity floor) | 2.5 |
| `YSZ_FRAC_OF_REST` | 0.388 / (0.388 + 0.362) = 0.5173 |
| Φ_Ni target / achieved | 0.250 / 0.2502 mean |
| `STRUCT_SEEDS` | 0, 1, 2, 3, 4 |

### 4.2 Intervention

| group | `intended_T_vox` | achieved p10 ratio | neck p10 | neck p50 | radius shrink |
|---|---|---|---|---|---|
| base | — (NaN) | 1.000 | 120 nm | 320 nm | — |
| lower-tail | 8.5 | 1.333 | 160 nm | 320 nm | 1.70–2.13 % |
| lower-tail | 11.0 | 2.000 | 240 nm | 320 nm | 5.17–6.79 % |

### 4.3 D4 damage operator — `cmlib/damage.py:70`

| parameter | value | status |
|---|---|---|
| `p_erode` | **0.35** | **FROZEN** — E0 origin, not re-derived (§0g/2) |
| `expand_vox` | **1** | **FROZEN** — E0 origin, not re-derived (§0g/2) |
| dilation structure | `STRUCT6` (6-connectivity), YSZ never overwritten | definitional |
| intensity variable | `n_rounds` (integer) | the only variable swept |
| post-step | removal of Ni disconnected from the spanning backbone | definitional |
| `DAMAGE_SEEDS` | 200, 201, 202, 203, 204 — **independent of structure seeds** | frozen |

### 4.4 Bisection

| parameter | value |
|---|---|
| initial bracket | [1, 20], expand-only, cap 64 |
| termination | bracket width ≤ 1 |
| outcome | bracket midpoint (integer-bracketed → 1.0-round resolution floor) |
| interpretability threshold | **1.0 round** (§0g/1) |
| total runs | 5 structure seeds × 3 groups × 5 damage seeds = **75 bisections** |

### 4.5 Reference values for regression checks

- Base-only bisection (damage seeds 100–102, n = 15): mean **8.767**, sd 0.594,
  distribution 7.5×1 / 8.5×8 / 9.5×6, no bracket expansion on any run.
- p10-group (damage seeds 200–204): group means **8.54 / 8.50 / 8.50**;
  base-structure sd **0.200**; 74 of 75 midpoints exactly 8.5.
- Autopsy D1 at n = 8: retained **0.4186 / 0.4142 / 0.4139**.
- D4loc removed-voxel mean EDT: **20.00 nm** (base) / **20.02 nm** (2.00×);
  all-voxel n = 8 median 28.28 nm.
- Provenance check, structure seed 0: **67,803** voxels differ between base and
  2.00×; **46 of 224** necks changed; base-only 33,909 / high-only 33,894.

---

## 5. Known invalid or non-informative diagnostics

**Do not use these as evidence in either direction. None will be repaired or
re-run as part of this project.**

1. **D2 — EDT percentile over all Ni voxels (`null_autopsy.csv`,
   `edt_p10` / `edt_p25` columns): DEGENERATE.** Returns identical values for
   every group at every damage state **including n = 0**, where the structures
   are verified to differ in 67,803 voxels with neck p10 of 120 vs 240 nm. The
   proxy measures surface fraction, not neck thickness: 20.0 nm is exactly one
   voxel, and >10 % of voxels in any solid lie within one voxel of a surface. It
   contributes **no information**. A valid lower-tail diagnostic would need to be
   restricted to the neck population (e.g. SNOW throat inscribed diameters on the
   damaged network) — a new metric, not run.
2. **`cpsd_r50max` size metric: FAILED CONVERGENCE.** Non-monotone; no resolution
   with both seeds under 0.5 pp (`cpsd_convergence.csv`).
3. **Raw EDT as a stand-in for particle size: INVALID.** ~3× compressed relative
   to generator radius; not equivalent to local thickness either
   (`cpsd_candidate_precheck.csv`).
4. **Local thickness: UNCONVERGED** at the ladder lengths tested.
5. **Skeleton-based graph metrics (Phase 3, real data): GATE FAILURE.**
   `skeletonize` returned a curve skeleton for the fine anode (83 % degree-2) but
   a medial *sheet* for medium and coarse (79 % / 94 % degree ≥ 4, spiking at
   exactly 8). Skeleton dimensionality varied monotonically with the variable
   under test. Reported and replaced with SNOW, not worked around.
6. **λ₂ raw as a connectivity measure: CONFOUNDED with node count** (λ₂ ∝
   N^(−0.867); λ₂·N collapses a 7.65× spread to 1.24×). Usable only as a
   coarseness proxy, and reported as such.
7. **Synthetic TPB density: NOT REAL-DATA COMPARABLE** — 19.4–20.5 µm⁻² vs a
   real range of 1.07–2.65 µm⁻² (7.3–19.2×). Fit for internal percolation/damage
   work only.
8. **Any p-value or correlation coefficient on the n = 3 real-data comparison:
   FORBIDDEN by design**, and enforced in `phase6_verdict.py`.

---

## 6. Known limitations

Full statements in `PATH_B_MEMO.md` §F. Summary:

1. **D4 is phenomenological**; `p_erode` and `expand_vox` were frozen, not
   derived. The null is specific to this operator at this parameterization.
2. **Outcome resolution is 1.0 damage round** — sub-round effects are invisible
   by construction.
3. **Synthetic TPB is 7.3–19.2× real.**
4. **Image-based size comparability (P2-C2) is deferred** — no converged
   image-based size metric backs the size-matching claim; generator radius is
   exact ground truth for the *input parameter* only.
5. **Radius-shrink confound (5.17–6.79 % at 2.00×) is moot but recorded**; any
   future positive on this platform must address it first.
6. **D2 is invalid** — the conclusion rests on D1, D3, D4loc.
7. **No real redox calibration** — real data was a qualitative anchor only.
8. **No Family C / hierarchical architecture test.**
9. **Real-data limits:** coarse anode sub-REV for network metrics (one ROI had no
   spanning cluster); medium vs coarse retention within noise; TPB ground truth
   digitized from a bar chart (self-consistent to 0.7 pp); post-redox medium and
   coarse stacks 40 % voxel-anisotropic (SNOW not applied to them); conductance
   weighting (throat area / throat length, intrinsic σ dropped) is a modelling
   choice, with unweighted `neck_nm` also reported.
10. **The porespy `snow2` chunked-watershed bug** was active in Phase 3/4;
    assessed, no conclusion changes.
11. **n = 3** on the real data — directional only.

---

## 7. Explicitly not run

Recorded so the closure is a scoping decision on the record, not an omission:
further generator variants; alternative or neck-selective damage models; Family
C / hierarchical architectures; real-data calibration of damage or size metrics;
opening granulometry; any third image-based size metric; instruction item 6; any
repair or re-run of D2; any change to D4 or the generator.
