# Platform v2 — Ni generator design memo

Ni-only generator qualification. Per instruction: no YSZ/pore placement, no
D4/damage modelling, no Family C, no calibration against the real
Holzer/Pecho dataset in this phase. Code: `cmlib/synth.py`
(`platform_v2_lattice_geometry`, new; all mass-conservation and
percentile-targeting functions reused unmodified from Family B/T5b),
`scripts/platform_v2/design_probe.py`, `scripts/platform_v2/qualification_run.py`.

## Geometry choice

Perturbed (jittered) simple-cubic sphere lattice, same construction family as
the Family B pilot, generalised to an **asymmetric particle count per axis**
(`platform_v2_lattice_geometry`, new function):

| parameter | value |
|---|---|
| sphere radius R | 12.1 voxels (diameter 24.2 voxels = 484 nm at 20 nm/voxel) |
| pitch | 32 voxels |
| particles along z (percolation axis) | 6 |
| particles along y, x | 4 |
| margin (y/x only) | 8 voxels |
| positional jitter | ±15% of pitch (interior z, all y/x) |
| domain shape | (161, 168, 168) voxels — every axis inside the 160–192 target |
| n_sites / n_pairs | 96 / 224 |

**Why the particle count is asymmetric (6 vs 4), not a mistake.** The
domain-spanning guarantee (validated since T5b/Family B) requires the
z-boundary sphere layers to sit exactly on the domain faces:
`nz = (nlat_z−1)·pitch + 1`, no margin term. The y/x extent instead carries a
margin on both sides: `n_xy = nlat_xy·pitch + 2R + 2·margin`. At a pitch set
by the Φ_Ni target, these two formulas cannot both land in [160,192] with the
*same* particle count — solved by decoupling `nlat_z` from `nlat_xy` rather
than distorting pitch/margin to force a false symmetry. `nlat_z=6, nlat_xy=4`
was the first combination tried that satisfied the domain-range constraint
(no further search was needed).

## Topology strategy — the open question, resolved empirically

**Plain nearest-neighbour (6-connectivity) lattice adjacency was tested
first, per instruction, before any topology modification was designed.**
Result (`scripts/platform_v2/design_probe.py`, full log in
`design_probe_log.txt`):

- Raw topological mean degree (pairs-list, before any SNOW/watershed effect):
  2×224/96 = **4.667**.
- SNOW-measured mean degree on the actual rasterized base structure:
  **4.204** — inside the target band [3.5, 4.5].
- Degree distribution: min 1, max 7, median 5; 14/98 nodes (14.3%) at degree
  1 (pendant chambers), 0 at degree 0 (no isolated fragments).

**Finding, stated plainly: plain lattice adjacency already satisfies the
coordination target. No topology modification (face-diagonal bonds,
proximity edges, or any other mechanism) was built.** This confirms the
premise in the review: the earlier "coordination too low, needs topology
work" framing was an artifact of the ungrounded 6–8 target, not a real
platform limitation.

## Φ_Ni control method

Φ_Ni is controlled by the sphere radius **R relative to pitch**, independent
of topology (which pairs get a neck) — the two were conflated in the Family
B pilot (R/pitch=0.5, spheres exactly tangent) but are genuinely separate
knobs in this generator. R was tuned by direct numerical search (rasterize
and measure, not analytic approximation — the same discipline used to tune
the base neck-width mixture in Family B):

| R (voxels) | measured Φ_Ni (spheres + base necks) | deviation from 0.250 |
|---|---|---|
| 12.0 | 0.2491 | −0.38% |
| **12.1** | **0.2502** | **+0.10%** |
| 12.2 | 0.2515 | +0.59% |
| 12.3 | 0.2529 | +1.16% |

R=12.1 was selected: essentially exact, and diameter (24.2 voxels) still
inside the 24–32 voxel target range. Note the necks themselves contribute
substantially to Φ_Ni: spheres alone give only Φ_Ni=0.1323 at this R — the
224 base necks add roughly 0.118 (about 47% of the total Ni volume). This
means the base neck-width distribution and Φ_Ni control are **not
independent** and were tuned together, not sequentially.

## Neck-distribution design

Reused unmodified from Family B: a mixture (not a plain uniform draw — a
uniform draw provably cannot reach p50/p10 ≥2.5 at a 4-voxel resolution floor
within a physically sensible neck-width range). Ranges scaled down from the
Family B pilot's (R=14) design in proportion to the new R=12.1
(scale factor ≈0.86, rounded to clean voxel counts):

| | Family B pilot (R=14) | Platform v2 (R=12.1) |
|---|---|---|
| weak population | 4–6 vox, 20% | 4–6 vox, 20% (floor unchanged — cannot shrink further) |
| normal population | 13–22 vox, 80% | 12–20 vox, 80% |

Measured on the representative base draw (n_pairs=224, larger than the
pilot's 144 — less sampling noise expected): p10=5.0 vox, p50=15.0 vox,
**p50/p10=3.00** — inside the target 3.0–4.3 range, comfortably above the
2.5 validity floor, and base p10 resolved at 5 voxels (exceeds the ≥3,
preferably ≥4 requirement).

## Domain size

161×168×168 voxels (4.54 Mvoxel) — inside the 160³–192³ target range on
every axis. Roughly 2.1× the Family B pilot's domain (2.18 Mvoxel); SNOW
extraction cost scales accordingly (~14–40 s per structure measured here vs
~5–15 s for the pilot).

## Scope decision carried forward, restated

Per the explicit instruction: **medium/coarse-like neck geometry, medium
real anode as primary anchor. Fine-like geometry is out of scope for
platform v2.** This was already the implication of using R≈24–32 voxel
particles (adequate resolution for medium/coarse-like p10/diameter ratios
~0.11–0.12; fine-like ratios ~0.06 would need ~50–66 voxel particles, pushing
well past the 192³ ceiling — see `preregistration.md` §0e). No fine-anode
comparison is attempted or implied by this platform.

## Expected achievable p10 rungs

Confirmed by the qualification run (`scripts/platform_v2/qualification_run.py`,
5 seeds, n_pairs=224): as in Family B, achieved values cluster at discrete
rungs set by the base distribution's own value grid, not at the nominal
targets.

| nominal target | achieved p10 ratio (5 seeds) |
|---|---|
| 1.45× | 1.33, 1.33, 1.33, 1.33, 1.33 (all five identical) |
| 2.0× | 2.00, 2.00, 2.00, 2.00, 2.00 (all five identical) |

Every seed landed on exactly the same achieved rung within its nominal bin —
tighter clustering than Family B (which had scattered rungs, e.g. 1.33/1.61
within its "1.5×" bin). This particular base mixture and n_pairs=224 produces
an unusually clean, reproducible two-rung ladder (1.33×, 2.0×) rather than a
continuum. Whether that generalises or is a property of this specific
mixture/population size is untested.

## Target comparison table against the real medium anode

| quantity | platform v2 (achieved) | real medium anode (measured) | source |
|---|---|---|---|
| Φ_Ni | 0.2502 (mean of 5 base seeds, range 0.2474–0.2561) | 0.250 (design target) | `phase4c`/`REPORT.md` |
| SNOW mean_degree | 4.193 (range 4.184–4.206), stable 4.10–4.23 across all widening levels | 3.41 ± 0.28 | `phase4c_metrics_per_anode_8.0um.csv` |
| p50/p10 | 2.90 (range 2.50–3.00) | 3.29 | same |
| particle diameter | 24.2 voxels = 484 nm | 1445 nm (watershed, volume-weighted) | `REPORT.md` |

**Two honest gaps, not glossed over.** (1) Platform v2's coordination (4.19)
sits above the real medium anode's own measured value (3.41) — inside the
approved 3.5–4.5 band, but on the high side rather than centred on the real
figure; the band was set from the *range* across all three real anodes
(3.41–3.65), not tuned to medium specifically. (2) p50/p10 mean (2.90) falls
slightly *below* the 3.0–4.3 target range, though every individual seed
clears the 2.5 validity floor (one seed sits exactly at 2.50). Neither gap
blocked Gate P2-A or P2-B (see qualification report), but both are worth
noting for anyone reading this table as a claim of a tight match — it isn't
one, by design (this platform targets medium/coarse-*like* geometry, not a
calibrated fit to the medium anode specifically, which is explicitly out of
scope until real-dataset calibration is approved).

**The load-bearing finding is in Gate P2-C, not this table** — see
`qualification_report.md`: lower-tail widening at both tested ratios drives
a c-PSD deviation outside ±5% for most seeds, traced to the mass-conservation
radius-shrink mechanism having less headroom at Φ_Ni=0.250 than it had at the
Family B pilot's Φ_Ni≈0.32–0.33. Composition and topology (Gate P2-A) pass
cleanly; the decoupling mechanism itself (Gate P2-C) does not yet, at this
Φ_Ni.
