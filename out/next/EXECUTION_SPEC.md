# Execution spec — neck-size decoupling study (Phases −1 → 5)

Status (2026-08-10): **G0, G−1, G0c (T5b), G0d (Family B pilot) PASSED. E0
vertical-slice spike run (`out/spike/`) — exploratory/non-confirmatory,
returned an honest INCONCLUSIVE result. Platform-v2 targets corrected against
real data (§0e): coordination target moved from ungrounded "~6-8" to measured
"~3.5-4.5". Platform v2 Ni generator qualification RUN (§0f,
`out/platform_v2/`): Gate P2-A (composition/topology) PASSED cleanly — plain
lattice adjacency already meets the coordination target, no topology
modification built. Gate P2-C (lower-tail decoupling) NOT YET FEASIBLE at
either tested ratio — traced to a c-PSD violation caused by reduced
mass-conservation headroom at the lower Φ_Ni=0.250 target (not a bug, not
retuned). YSZ/pore placement, damage models, Family C, and real-dataset
calibration remain OUT OF SCOPE until Gate P2-C's failure is reviewed.**

Repository root for all paths below:
`C:\Users\ARNAV MEHTA\Downloads\soec\voxel-damage-operators-niysz\`

---

## 0. Scope restatement

At **fixed Ni volume fraction** and **fixed mean Ni particle size**, does
engineering the neck-size distribution — especially reducing the vulnerable
narrow-neck tail, or adding a coarse percolating Ni backbone — *independently*
improve retained Ni electronic percolation under redox-like damage, without
unacceptable loss of initial TPB density?

This is a **causal/design** question, not a predictive one. The prior study
failed because in n=3 real anodes every candidate metric was collinear with
coarseness. Synthetic structures exist solely to break that collinearity.

Two acceptable outcomes: **Path A** (neck engineering gives an independent,
robust gain) or **Path B** (once particle size and loading are fixed, no
independent effect, or decoupling is not achievable in this framework).

---

## 1. Risks that shape the design

### R1 — Widening necks raises *measured* particle size (constraint may be self-defeating)
Watershed cuts particles at necks. Widening necks merges neighbouring regions,
so watershed mean particle size rises as a mechanical consequence. The Family B
gate (neck p10 ×1.5–2, particle size ±5%) may be unsatisfiable **by
construction**.

**Mitigation / primary design axis:** vary the *lower tail only*. Selectively
widen the narrowest necks (bottom ~20% by throat radius) and leave the rest
untouched. This moves neck p10 by ~2× while altering a small fraction of Ni
volume and barely perturbing neck p50 or watershed particle size.
Report **three** size measures for every structure so the constraint is
auditable and not hostage to one definition:
  - `d_watershed_volwt` — volume-weighted equiv. sphere diameter (primary, as
    used in phase4d on real data)
  - `d_cPSD_r50max` — opening/continuous-PSD based size (neck-insensitive)
  - `d_nominal` — generator sphere diameter (ground truth)
Pre-register: primary constraint on `d_watershed_volwt`; if it cannot be held
within ±5% while `d_cPSD_r50max` and `d_nominal` *are* held, that is reported as
a finding about the measure, not a silent switch.

### R2 — Damage model D1 is near-circular
"Destroy narrow necks → structures with fewer narrow necks survive" is close to
a tautology. **Recommend promoting D4 (morphological redox surrogate) to
REQUIRED and demoting D1 to a sanity/upper-bound model.** In D4 neck
sensitivity is emergent, not asserted.
Pre-register: an effect appearing under D1 but not under D4 **and** D2 is
scored **negative**.

### R3 — Resolution vs REV squeeze; 128³ is too small for the main sweep
From the real data: neck p50 / particle diameter ≈ 0.26–0.37, but neck **p10** /
particle diameter ≈ 0.06–0.12. Resolving a p10 neck at ≥3 voxels needs particle
diameters ~25–50 voxels; ≥8 particles per side then needs ~200–400 voxels/side.
At 128³ with 8 particles across, p10 necks are ~1 voxel — the quantity under
test becomes quantization noise.

**Recommendation:** validation at 64³/100³; **main sweeps at 192³**, particle
diameter 24–28 voxels, ~7–8 particles per side. Note the binding constraint has
inverted vs. the original study: 192³ = 7.1 Mvox is trivial for memory
(~16 MB bool); **wall-clock SNOW time is the cost driver** (~10 s/extraction).

> **Amendment (post-E0 review, see `preregistration.md` §0e):** re-computed
> precisely from `phase4c_metrics_per_anode_8.0um.csv` / `REPORT.md`:
> p10/diameter = 0.0606 (fine) / 0.1086 (medium) / 0.1197 (coarse) — a ~2x
> range, not resolved by a single particle-size choice. 24–28 voxels (used by
> the Family B pilot, up to 32 in some framings) resolves medium/coarse-like
> necks at ≥3 voxels but under-resolves fine-like necks by ~2x (fine needs
> ~50–66 voxels for 3–4 voxel neck resolution). **Platform v2 must state
> explicitly which regime it targets** (medium/coarse-like at ~24–32 voxels,
> or fine-like at ~50–66 voxels, pushing well past 192³) — this was
> previously implicit, not a deliberate scope choice.

### R4 — Damage intensity can manufacture or erase the effect
Mild damage → everything retains ~1.0. Severe → everything → 0. The comparison
point must be pre-registered.
**Pre-register:** primary comparison at the damage intensity where the *Family B
mean* retained `P_span` falls in **[0.5, 0.8]**, selected from Family A
calibration and then **frozen**; the full intensity sweep is reported regardless.

> **Amendment (post-E0 review):** this risk was fully realised — E0's D4 on
> the pilot geometry showed a step transition (P_span exactly 1.0 or exactly
> 0.0, no intermediate value at any of 6 tested intensities) with zero
> [0.5, 0.8] regime ever observed. **Superseded for the platform-v2
> calibration pilot** by the per-structure bisection procedure in
> `preregistration.md` §0e ("Damage-intensity bisection procedure"), which
> brackets each structure's own transition individually rather than
> assuming one shared intensity lands every structure in [0.5, 0.8]
> simultaneously.

### R5 — Statistics: synthetic replication is legitimate, the n=3 anchor is not
Seeds give genuine replication, so effect sizes and CIs **across seeds** are
valid for synthetic structures. The n=3 real anodes remain a **qualitative
anchor only** — no p-values, no correlation coefficients, and the two must never
be pooled. (Flagged as an open question to confirm.)

---

## 2. What already exists vs. what must be written

`cmlib/` currently provides:

| module | provides | reusable as-is? |
|---|---|---|
| `percolation.py` | `structure_for`, `label_phase`, `spanning_labels`, `percolates`, `percolation_report`, `percolating_mask` | yes — but **no `P_reach`** |
| `pnm.py` | `extract_ni_network(ni_mask, spacing_nm, axis, connectivity, r_max, sigma, parallel_kw)` → `(G, diag, extras)` | **yes, directly** — now defaults to serial SNOW extraction (`parallel_kw=None`) after a porespy chunked-mode bug was found; see `IMPACT_NOTE_porespy_parallel_bug.md` |
| `tpb.py` | `tpb_density_volume(vol, lab, spacing_nm)` — already ternary | **yes, directly** |
| `metrics.py` | `summarise_network`, `algebraic_connectivity`, `mincut_between_faces`, `effective_conductance`, `neck_quantiles` | **yes, directly** |
| `phases.py` | `assign_labels`, `volume_fractions` | real-data metadata path only |
| `roi.py`, `io.py`, `ground_truth.py` | tiling, TIFF streaming, published values | not needed for synthetic |
| `graph.py` | skeleton pipeline | **DEPRECATED — do not use for cross-sample ranking** |

**Gaps to close (this is the whole Phase 0 refactor):**

1. `P_reach` / `P_largest` logic currently lives inside
   `phase5_percolation.py::analyse()` — a script, not the library.
2. Watershed particle sizing lives inside
   `phase4d_particles.py::watershed_particles()` / `size_stats()` — also a script.
3. No ternary-volume volume-fraction helper (existing one is tied to dataset
   metadata).
4. No single-call API surface.

### New/modified library files

```
cmlib/percolation.py   MODIFY  + percolation_summary(mask, axis, conn)     DONE
                                 -> P_span, P_reach, P_largest, n_clusters, percolates
                                 (moved from phase5_percolation.analyse; that
                                 script now delegates and was regression-
                                 tested bit-identical against its prior output)
cmlib/particles.py     NEW     watershed_particles(), size_stats(),         DONE
                                 cpsd_r50max() (moved from phase4d_particles,
                                 regression-tested bit-identical; c-PSD added,
                                 see the porespy local_thickness unit-handling
                                 trap documented in its module docstring)
cmlib/synthvol.py      NEW     ternary volume container +                   DONE
                                 volume_fractions_from_volume(),
                                 save_ternary()/load_ternary() (T6 round trip)
cmlib/api.py           NEW     the five wrappers (below)                    DONE
cmlib/synth.py         NEW     T5b generator scaffolding: mass-             DONE 2026-08-10
                                 conservative + percentile-targeted neck
                                 widening on a cubic sphere lattice.
                                 NOT the full Family A/B/C generator (no
                                 config, no YSZ/pore placement, no damage
                                 model) -- scoped to T5b only, per the
                                 post-T5 review decision (preregistration.md
                                 amendment). Family B proper still AWAITS
                                 APPROVAL, gated on T5b's outcome.
cmlib/damage.py        NEW     Phase 2 damage operators      [AWAIT APPROVAL, gated on T5b]
```

### The five wrapper functions (`cmlib/api.py`)

```python
extract_network(binary_ni, spacing_nm, *, axis=2, r_max=4, sigma=0.4)
    -> (G, diag)                       # thin wrapper over cmlib.pnm

compute_percolation(ternary_volume, labels, *, axis=2, connectivity=6)
    -> {P_span, P_reach, P_largest, percolates, n_clusters}

compute_tpb(ternary_volume, labels, spacing_nm)
    -> {tpb_density_um-2, tpb_length_um, volume_um3, edges_z/y/x}

compute_particle_stats(ni_volume, spacing_nm, *, min_distance=4, sigma=0.4)
    -> {d_watershed_volwt, d_watershed_mean, d_watershed_median,
        d_cPSD_r50max, n_regions_used, n_regions_border_excluded}

compute_network_metrics(G, face_lo, face_hi)
    -> {lambda2_raw, lambda2_norm, mincut, g_eff,
        neck_p10/p25/p50/p90, n_nodes, n_edges, mean_degree}
```

All five take/return plain dicts so every phase writes one tidy CSV row per
structure.

### Script layout

```
scripts/next/phase_minus1_prior_art.py          # search-log helper (optional)
scripts/next/phase0_validate_synthetic_pipeline.py
scripts/next/phase1_generate_synthetic.py       [AWAIT APPROVAL]
scripts/next/phase2_damage_models.py            [AWAIT APPROVAL]
scripts/next/phase3_decoupling_experiment.py    [AWAIT APPROVAL]
scripts/next/phase4_robustness.py               [AWAIT APPROVAL]
configs/next/synthetic_generator.yaml
configs/next/damage_models.yaml
out/next/...
```

Every script bootstraps the repo root from `__file__` (two levels up) so it runs
from the repository root. `pyyaml` is the one new dependency (check first).

---

## 3. Output tables (canonical schemas)

**`out/next/structures.csv`** — one row per generated structure:
`struct_id, family, seed, config_hash, nx,ny,nz, voxel_nm,
phi_Ni, phi_YSZ, phi_pore,
d_watershed_volwt, d_cPSD_r50max, d_nominal,
neck_p10, neck_p25, neck_p50, mean_degree, n_nodes, n_edges,
P_span, P_reach, P_largest, tpb_density, gen_params_json`

**`out/next/damaged.csv`** — one row per (structure × damage model × intensity × seed):
`struct_id, damage_model, intensity, damage_seed,
ni_loss_frac, phi_Ni_post,
P_span_post, P_reach_post, tpb_post,
P_span_retained, P_reach_retained, tpb_retained,
neck_p10_post, n_nodes_post`

**`out/next/effects.csv`** — Phase 3 contrasts:
`family, damage_model, intensity, neck_p10_low, neck_p10_high,
d_retained_P_span, seed_noise_floor, passes_criterion, monotonic`

---

## 4. Go / no-go gates

| gate | criterion | fail action |
|---|---|---|
| **G−1** | No direct prior-art equivalent, or claim narrowed | narrow scope, re-register |
| **G0** | Ternary pipeline exact on analytic cases (§6) | fix before any generation |
| **G0b (T5, informal)** | naive neck-widening scaffold decouples neck p10 from particle size at fixed Φ_Ni | **FAILED 2026-08-10** for two specific, fixable reasons (Ni mass not conserved; "widen bottom 20%" cannot move measured p10 by construction) — see `out/next/t5_coupling_decision_report.md`. NOT scored as Path B; superseded by G0c |
| **G0c (T5b)** | corrected (mass-conservative, percentile-targeted) scaffold meets `out/next/preregistration.md` §0b acceptance criteria | **PASSED 2026-08-10** — under the strict per-seed gating in §0c amendment B, lower-tail is feasible (≥4/5 seeds) at BOTH 1.5× and 2.0×; 2.5× and all of uniform mode fail 0/5. Accepted as a real partial success; proceed toward Family B WITH the mandatory amendments in §0c. See `out/next/t5b_coupling_decision_report.md` |
| **G0d (Family B disordered pilot)** | disordered/perturbed-lattice pilot (§0d) meets the SAME §0c amendment B per-seed criteria at {base, 1.5×, 2.0×}, 5 seeds | if it fails, diagnose and iterate on the pilot BEFORE any larger Family B run — do not proceed to YSZ/pore placement or damage models regardless |
| **G0e (E0 vertical-slice spike)** | exploratory pipeline spike (YSZ/pore placement, D4, retained-metric recomputation) runs cleanly and non-confirmatory rules are honoured | **RAN 2026-08-10, INCONCLUSIVE** (no differentiation by achieved p10 ratio; correctly barred from declaring Path A/B). Review of its platform-v2 follow-on found the coordination target ungrounded/backwards and other targets loosely grounded — corrected in `preregistration.md` §0e |
| **G0f (Platform v2 Ni qualification)** | Gate P2-A (composition/topology), P2-B (base distribution), P2-C (lower-tail decoupling) per `preregistration.md` §0f | **P2-A PASSED cleanly 2026-08-10** (Φ_Ni=0.2502, mean_degree=4.193 in-band, plain topology sufficient — no modification built). **P2-B validity floor PASSED, target-range mean missed** (2.90 vs 3.0-4.3). **P2-C NOT FEASIBLE at either ratio** — traced to c-PSD failing (−9% to −12%) from reduced mass-conservation headroom at Φ_Ni=0.250; not a bug, not retuned. See `out/platform_v2/qualification_report.md`. YSZ/pore placement and damage models remain blocked on this gate |
| **G1** | ≥5–10 Family B realizations; neck p10 spans ≥1.5× (target 2×); Φ_Ni within ±0.005 abs (hard ceiling ±5% relative, target ±2%, per preregistration §0b — **c-PSD or generator-known size**, not watershed alone); pristine `P_span` ≥ 0.95 recorded | **do not claim physical impossibility** — report "not testable in this framework", start Path B |
| **G2** | ≥1 damage model reproduces fine-worst percolation loss; produces valid ternary voxel output; not a single-parameter artifact | stop Path A, Path B memo |
| **G3** | Δ retained `P_span` ≥ `max(0.05, 3×noise_floor)`; monotonic/near-monotonic; survives ≥1 alternative damage model (must include D4 or D2, not D1 alone); initial TPB loss ≤20% vs best same-φ_Ni baseline | score negative, Path B |
| **G4** | Effect survives seeds, intensity sweep, alt damage model, domain size, `r_max` sweep | downgrade to "suggestive, not robust" |

**Noise floor definition (pre-registered):** pooled within-structure standard
deviation of retained `P_span` across ≥5 damage seeds at fixed generator
parameters and fixed damage intensity.

---

## 5. Phase −1 — prior-art search queries

Run against Google Scholar / Web of Science / Scopus, plus forward-citation
sweeps of the two Pecho 2015 papers and Holzer 2013 JPS.

**Cluster 1 — neck-size engineering / decoupling (closest to our claim)**
1. `Ni-YSZ anode "neck" size distribution redox tolerance microstructure`
2. `sintering neck radius independent particle size percolation SOFC electrode simulation`
3. `"neck size" OR "contact radius" effective conductivity Ni cermet degradation`
4. `decouple neck size particle size percolation porous electrode`

**Cluster 2 — bimodal / backbone architectures**
5. `bimodal Ni particle size distribution SOFC anode redox cycling`
6. `"coarse backbone" OR "scaffold" fine decoration Ni-YSZ TPB percolation`
7. `infiltrated Ni backbone SOFC anode redox stability microstructure`

**Cluster 3 — TPB vs percolation tradeoff under redox**
8. `triple phase boundary percolation tradeoff Ni-YSZ redox cycling quantitative`
9. `active TPB loss Ni percolation loss redox cycles tomography`

**Cluster 4 — synthetic microstructure generation for Ni redox damage**
10. `synthetic microstructure generation Ni-YSZ anode percolation TPB sphere packing`
11. `stochastic geometry model SOFC electrode redox degradation simulation`
12. `phase-field OR DEM sintering Ni-YSZ electrode neck growth effective transport`
13. `NiO reduction volume change 41% microstructure damage model voxel simulation`

**Cluster 5 — named groups likely to have priority**
14. `Shikazono OR Kishimoto Ni-YSZ synthetic microstructure percolation TPB`
15. `Holzer OR Pecho constrictivity bottleneck Ni-YSZ redox` (+ citing works)
16. `Cronin OR Wilson OR Barnett Ni-YSZ redox cycling 3D microstructure`
17. `Neumann OR Schmidt stochastic 3D microstructure model SOFC anode`

**Deliverables:** `out/next/phase_minus1_prior_art.md` (memo + hit table +
explicit "closest prior work and how we differ"),
`out/next/preregistration.md` (criteria frozen **before** Phase 3).

---

## 6. Phase 0 — synthetic validation tests

All on hand-built tiny ternary volumes (≤64³) with analytically known answers.
Extends, and does not replace, the already-passing `probe_metrics.py`,
`phase4a_validate_tpb.py`, `probe_porespy.py`.

**T1 — volume fractions exact.** Ternary block with prescribed voxel counts.
Assert φ recovered exactly (integer counts / total).

**T2 — TPB exact, ternary, anisotropic.** Reuse the Phase-4a analytic cases:
three phases meeting along one axis-aligned line → TPB length exactly `n·d`,
zero on the other two axes. Add a **new** case with all three phases in a
2×2 checkerboard column to confirm multiple parallel TPB lines add linearly.
Also re-confirm the measured staircase bias (1.713 on the (1,1,1) case).

**T3 — percolation, ternary input.** (a) Solid Ni slab spanning x → `P_span`=1,
`P_reach`=1. (b) Ni slab touching only x=0 → `P_span`=0, `P_reach`=1,
`percolates`=False. (c) Isolated central Ni cube → both 0. (d) Two parallel
bars, one spanning one not → `P_span` = spanning volume fraction exactly.
**(b) is the test that distinguishes `P_span` from `P_reach` and must be in.**

**T4 — SNOW throats on known shapes.** Two cubes joined by an `w`-voxel square
bar → `throat.inscribed_diameter` = `w·d`, `cross_sectional_area` = `(w·d)²`,
2 nodes / 1 throat. Sweep `w ∈ {4,6,8}`. (The `w`=6 case already passes in
`probe_porespy.py`; generalize and assert.)

**T5 — particle stats on known packing.** N non-touching spheres of radius R →
`n_regions_used` = N (minus border-truncated), `d_watershed_volwt` ≈ 2R within
one voxel. Then add a thin neck and confirm the **documented** merge behaviour
(this is the R1 quantification: measure how much `d_watershed_volwt` moves as
neck width grows at fixed R — a required input to the Family B design).

**T6 — round-trip.** A generated structure saved and reloaded reproduces every
metric bit-identically; config hash + seeds recorded.

**Deliverables:** `scripts/next/phase0_validate_synthetic_pipeline.py`,
`out/next/phase0_validation_report.md`.
**Gate G0:** all of T1–T6 pass; T5's neck-vs-size coupling curve is quantified
and carried into the Phase 1 design.

---

## 7. Minimal implementation order (this session's proposal)

1. `cmlib/percolation.py` += `percolation_summary` (move from phase5 script) — **done**
2. `cmlib/particles.py` (move from phase4d script; add `cpsd_r50max`) — **done**
3. `cmlib/synthvol.py` (tiny; ternary container + φ helper) — **done**
4. `cmlib/api.py` (five wrappers) — **done**
5. `scripts/next/phase0_validate_synthetic_pipeline.py` (T1–T6) — **done, gate G0 PASS 31/31**
6. Phase −1 searches → prior-art memo → pre-registration — **done, gate G−1 PASS, no direct equivalent found**
7. T5 coupling-decision experiment — **done, FAILED for two specific fixable reasons** (see
   `out/next/t5_coupling_decision_report.md`); Q1/Q2 below are resolved by this, not left open.
8. **T5b (this amendment): `cmlib/synth.py`** (mass-conservative + percentile-targeted
   widening) + `scripts/next/t5b_coupling_experiment.py`, gated by
   `out/next/preregistration.md` §0b before any Family B / full generator work.

Nothing in steps 1–7 was a large run: all tests ≤64³ (T1-T6) or a single small
lattice (T5/T5b, ~2 Mvoxel). T5b is bounded (~30-40 structures, background run).

---

## 8. Open questions

**Q1, Q2 — RESOLVED 2026-08-10 by the T5 experiment and the review decision**
(see `out/next/preregistration.md` §0 for the frozen replacement text):
the primary axis is percentile-targeted lower-tail widening (not naive
"bottom 20%"), made mass-conservative (not naive material addition), and the
particle-size gate is c-PSD / generator-known geometry (not watershed alone).
T5b tests this corrected design. Superseded, not deleted, for the audit
trail: the original Q1/Q2 text below shows what was asked before T5 ran.

<details><summary>original Q1/Q2 (superseded)</summary>

**Q1 (blocking Family B design).** Confirm the primary Family B axis:
selectively widen only the **narrow tail** (moves p10 ~2×, leaves p50 and
particle size ~fixed), rather than shifting the whole neck distribution? The
latter will very likely fail the ±5% particle-size gate for the reason in R1.

**Q2 (blocking, methodological).** If watershed particle size cannot be held to
±5% while neck p10 moves 2× — even under the tail-only design — do you want:
(a) report as "not decoupleable under the watershed size definition" and go to
Path B; (b) relax the primary size constraint to the neck-insensitive c-PSD
measure and report both; or (c) relax the gate to ±10%? My recommendation is
(a)+(b) reported together: it is the honest answer and it is itself a result.

</details>

**Q3.** Promote **D4 (morphological redox surrogate) to required** and demote
D1 to sanity/upper-bound, per R2? I recommend yes.

**Q4.** Main-sweep domain **192³** with 24–28-voxel particles (not 128³),
per R3? Confirm this is acceptable given no-HPC.

**Q5.** Confirm R5: seed-level statistics on synthetic structures are permitted
(effect sizes + CIs across seeds), while the n=3 real anodes stay qualitative
and the two are never pooled.

**Q6.** YSZ/pore placement rule — this materially sets TPB and is not yet
specified. Proposed default: second interpenetrating sphere packing for YSZ,
remainder pore, targeting the **medium** real anode
(φ_Ni 0.250 / φ_YSZ 0.388 / φ_pore 0.362). Confirm or specify otherwise.

**Q7.** Should Family A be tuned to reproduce the real anodes' *measured* neck
p10 (69.5 / 157 / 205 nm) and particle sizes (1148 / 1445 / 1715 nm) in
physical units, or run dimensionless with only the *ratios* matched? Given R3,
dimensionless-with-matched-ratios is far more tractable; physical-unit matching
would force ≥400³ domains.

**Q8.** Non-blocking: the folder is not under version control
(`git init` not run). Given "keep scripts reproducible", do you want it
initialized so configs/seeds are tracked?
