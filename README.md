# Voxel Damage Operators on Ni-YSZ Tomograms: Validation and Failure Modes

This repository contains the source code, pipelines, and manuscript write-up for a validation study of voxel-scale damage operators applied to segmented solid oxide cell electrode reconstructions.

The study is run against the Holzer/Pecho Ni-YSZ FIB-tomography dataset ([Zenodo 4056538](https://zenodo.org/records/4056538)) spanning three pristine and post-redox anode structures (fine, medium, and coarse).

---

## 1. Scientific Overview & The Six Failure Modes

Voxel-scale damage operators are widely used to simulate microstructural degradation directly on segmented tomograms. This project systematically tests these operators against a sharp measured signature: **the finest anode retains nickel percolation worst, yet retains triple-phase boundary (TPB) best.**

Our validation study shows that no standard operator class reproduces this signature due to six fundamental methodological and topological failure modes:

| # | Failure Mode | Impact & Description |
|---|---|---|
| **1** | **Planar Lattice Cuts** | Jittered regular lattices cut at exactly one flat cross-section with zero variance; they cannot model the local, disordered throat failures of real networks. |
| **2** | **Over-connected Gates** | Requiring perfect pristine connectivity in synthetic models conceals the fact that real electrodes carry $1.2\text{--}11.2\%$ disconnected nickel before service. |
| **3** | **Pruning Circularity** | Restricting outcomes to the largest connected component makes the spanning fraction ($P_{\text{span}}$) identically unity at the first step, hiding degradation. |
| **4** | **TPB Voxel Roughening** | Local voxel moves inherently pit the surface at the grid scale, inflating TPB density $3.7\text{--}15\times$ regardless of conservation budgets. |
| **5** | **Curvature-Rank vs. Area** | Stencil-based curvature proxies do not guarantee the specific surface area reduction assumed of them, rising at $n=1$ under standard stencils. |
| **6** | **Area Monotonicity Barrier** | Requiring monotonic area reduction strictly forbids the energy-raising or neutral steps ($80\%$ of accepted moves) needed to cross topological barriers and thin necks. |

For a complete pedagogical breakdown of the physics, mathematics, and operator behaviors, refer to the compiled [Standalone Primer](file:///c:/Users/ARNAV%20MEHTA/Downloads/soec/connectivity_margin/out/writeup/primer.html) (`out/writeup/primer.html`).

---

## 2. Repository File Map

To navigate the files in this repository, they are organized here by module and execution role:

### Core Support Library (`cmlib/`)
Located in [cmlib/](file:///c:/Users/ARNAV%20MEHTA/Downloads/soec/connectivity_margin/cmlib/):
* `damage.py` & `damage2.py`: Voxel damage operators (D4, O1, O2, O3, O5, etc.).
* `seqgreedy.py`: Incremental sequential greedy KMC and area-decreasing swap engines.
* `metrics.py`: Graph-theoretic metrics (algebraic connectivity, minimum-cut, conductance).
* `percolation.py`: 3D site/bond percolation solvers.
* `tpb.py`: Interfacial triple-phase boundary (TPB) estimators.
* `particles.py`: Particle size and volume-weighted diameter calculators.
* `synth.py` & `synthvol.py`: Synthetic FCC/jittered sphere and sintering-yield generators.
* `roi.py` & `project2.py`: Region of Interest boundary handling and coordinate helpers.
* `api.py` & `io.py`: File storage and stack indexing utilities.

### Phase 1: Predictive Metrics Pipeline
Root scripts running the original connectivity margin pipeline:
* `phase0_validate_percolation.py`: Code validation comparing site percolation thresholds against literature values.
* `phase1_download.py` to `phase1_build_ground_truth.py`: Automated fetching, extracting, and formatting of Zenodo and supplementary tables.
* `phase2_inspect_labels.py` & `phase2_volume_fractions.py`: Label verification and volume fraction calculations.
* `phase3_extract_network.py` & `phase3_snow_sensitivity.py`: Watershed network extraction (SNOW) and parameter sweeps.
* `phase4a_validate_tpb.py` to `phase4e_lambda2_scaling.py`: TPB validation, particle size ordering, and Laplacian eigenvalue scaling audits.
* `phase5_percolation.py`: Spanning and reachability analysis.
* `phase6_verdict.py`: Final report tables builder.

### Phase 2: Operator Validation & Primer (`scripts/project2/`)
Located in [scripts/project2/](file:///c:/Users/ARNAV%20MEHTA/Downloads/soec/connectivity_margin/scripts/project2/):
* `test_operators.py`: Verification suite checking operator volume conservation and monotonicity (runnable via `pytest`).
* `o7_gate_a1v2_real.py` & `o7_strict_inequality.py`: Evaluates area-decreasing swaps and strict-inequality controls on real tomograms.
* `o7_tiebreak_sensitivity.py`: Demonstrates the effect of random vs. LIFO tiebreaking.
* `o7_o5v2b_rerun.py`: Runs zero-temperature greedy KMC.
* `build_primer_html.py`: Builds the standalone HTML primer with embedded SVG plots.
* `primer_figures.py` & `primer_figures_change.py`: Generates the figures used in the primer.
* `audit_ysz_cluster_sizes.py` & `ni_vulnerability_audit.py`: Structural vulnerability diagnostics.

### Diagnostic Probes & Notes
Root files for memory, thickness, and Skan/Porespy inspection:
* `probe_local_thickness.py` & `probe_skan.py`: Local thickness and skeleton line checks.
* `probe_snow2_parallel_bug.py`: Demonstrates a parallel processing bug in SNOW2.
* `IMPACT_NOTE_porespy_parallel_bug.md`: Summary of the SNOW2 parallel dispatch issue.

---

## 3. Running the Code

### Environment Verification
To verify your Python environment carries the required scientific stack, run:
```bash
python check_env.py
```

### Running Unit Tests
Validate the damage operators, TPB counters, and percolation estimators using `pytest`:
```bash
pytest
```

### Building the Primer Report
To compile the standalone report and generate all primer figures:
```bash
python scripts/project2/build_primer_html.py
```
This generates:
* Markdown: `out/writeup/PRIMER_voxel_operators.md`
* Inlined HTML: `out/writeup/primer.html` (fully self-contained, no external requests).
