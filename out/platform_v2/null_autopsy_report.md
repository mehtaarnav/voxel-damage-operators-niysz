# Null autopsy — p10-group experiment

Run 2026-08-10, 441 s. Code: `scripts/platform_v2/null_autopsy.py`. Data:
`null_autopsy.csv`, `null_autopsy_localization.csv`, `null_autopsy_log.txt`.

**Secondary and explanatory only.** Does not change the frozen primary
outcome (no resolvable effect / inconclusive / not Path B / not a universal
negative). No D4 parameter changed: `p_erode=0.35`, `expand_vox=1` as frozen.
No new primary outcome introduced.

Primary diagnostic question: **at n_rounds=8 (the last spanning state), do
the widened structures still retain a measurable lower-tail geometric
advantage over base?**

## D1 — Remaining Ni volume

| group | achieved p10 | vox @ n=0 | vox @ n=8 (sd) | retained @ n=8 (sd) | vox @ n=9 | retained @ n=9 |
|---|---|---|---|---|---|---|
| base | 1.000× | 1,137,007 | 476,038 (14,674) | 0.4186 (0.0082) | 392,747 | 0.3454 |
| lower-tail | 1.333× | 1,137,025 | 471,078 (14,431) | 0.4142 (0.0079) | 386,737 | 0.3401 |
| lower-tail | 2.000× | 1,137,013 | 470,650 (11,539) | 0.4139 (0.0049) | 383,640 | 0.3373 |

Starting volumes are matched to ~18 voxels in 1.14 M (the mass-conservation
design working). By n_rounds=8 all three groups have lost ~58% of their Ni
and sit within **0.5 percentage points** of each other in retained fraction —
with group differences (0.0044, 0.0047) **smaller than the within-group sd
(0.005–0.008)**.

**The widened groups retain very slightly *less* volume, not more.** The
intervention does not produce a volume advantage at the collapse boundary.

## D2 — Lower-tail thickness proxy — **the diagnostic is degenerate, not null**

| group | EDT p10 @ n=0 | EDT p10 @ n=8 | EDT p25 @ n=0 | EDT p25 @ n=8 |
|---|---|---|---|---|
| base | 20.0 nm | 20.0 nm | 40.0 nm | 20.0 nm |
| lower-tail 1.33× | 20.0 nm | 20.0 nm | 40.0 nm | 20.0 nm |
| lower-tail 2.00× | 20.0 nm | 20.0 nm | 40.0 nm | 20.0 nm |

**Every value is identical across all three groups, including at n_rounds=0
where the structures are known to differ.** That is the tell: this proxy
cannot distinguish structures that we have *already verified* differ in
67,803 voxels with neck p10 of 120 vs 240 nm.

**Why:** 20.0 nm is exactly **one voxel**. The EDT over *all* Ni voxels is
dominated by the surface shell — in any solid, well over 10% of voxels are
within one voxel of a surface, so the 10th percentile pins to 1 voxel
regardless of internal geometry. The p25 at n=0 (40 nm = 2 voxels) shows the
same quantization one step up, and collapses to 1 voxel after damage because
erosion increases the surface-to-volume ratio.

**This must not be read as "the intervention was erased."** It is a
measurement failure of the proxy as specified: an EDT percentile over all
phase voxels measures surface fraction, not neck thickness. The instruction
specified this definition, and I ran it as specified — but reporting its
null as evidence would be wrong. **D2 contributes no information here.**

## D3 — Spanning-cluster size at n_rounds=8

| group | P_span @ n=8 | spanning vox @ n=8 (sd) | fraction of original (sd) | P_span @ n=9 |
|---|---|---|---|---|
| base | 1.000 | 476,038 (14,674) | 0.4186 (0.0082) | **0.04** |
| lower-tail 1.33× | 1.000 | 471,078 (14,431) | 0.4142 (0.0079) | 0.00 |
| lower-tail 2.00× | 1.000 | 470,650 (11,539) | 0.4139 (0.0049) | 0.00 |

At the last spanning state, **P_span = 1.000 for every group** — the entire
remaining Ni network is one spanning cluster in all cases, so there is no
backbone-size advantage to detect. Widened structures have marginally
*smaller* backbones, again within noise.

The `P_span @ n=9` column is informative: base retains a residual 0.04
(one seed still spanning at 9), both widened groups are at exactly 0.00. This
is the same −0.04-round difference reported in the primary result, seen from
the other side — and it is in the *unfavourable* direction for the
hypothesis.

## D4loc — Failure-step localization (the decisive diagnostic)

Voxels removed between n=8 and n=9, with their EDT measured in the n=8 mask
(structure seed 0, damage seed 200):

| group | n removed | frac of n=8 | removed p10 | removed p50 | removed p90 | removed mean | *all* n=8 p50 | *all* n=8 p90 | *all* n=8 mean |
|---|---|---|---|---|---|---|---|---|---|
| base | 83,020 | 0.17 | 20.0 | 20.0 | 20.0 | **20.00** | 28.28 | 74.83 | 38.36 |
| high 2.00× | 86,174 | 0.18 | 20.0 | 20.0 | 20.0 | **20.02** | 28.28 | 72.11 | 36.95 |

**This answers the question the autopsy was posed to answer.** The voxels
whose removal breaks spanning are **not drawn from a broad size range and
are not the lower-tail neck population in any structure-specific sense** —
they are, essentially exclusively, **surface voxels at exactly one voxel
depth** (p10 = p50 = p90 = 20.0 nm; mean 20.00/20.02 against a background
distribution whose median is 28.28 and p90 is ~72–75 nm).

That is D4's erosion operator by construction: each round strips surface
voxels with probability `p_erode`. **Collapse occurs when cumulative uniform
surface stripping has thinned the network globally — not when a specific
narrow neck is severed.** The removed set is statistically indistinguishable
between base and 2.00× (20.00 vs 20.02 nm), and the *fraction* removed is
nearly identical (0.17 vs 0.18).

## Recommended interpretation

Per the pre-registered rules, this is **branch A — "intervention erased
before transition"** — but with an important refinement that the rules did
not anticipate, and that I think matters more than the branch label:

> **The intervention is not so much *erased* as *bypassed*.** D1 shows both
> groups converge to within 0.5 pp of the same retained volume by n=8, and
> D4loc shows why: the failure mode is uniform surface erosion, which
> operates on total surface area and is indifferent to *where* the material
> sits in the neck-width distribution. A lower-tail intervention has nothing
> to act on if collapse is not mediated by the lower tail.

Recorded conclusion, in the pre-registered wording:

> **"The intervention does not survive to the collapse boundary under D4."**

With this addition, which the diagnostics support directly:

> Under this D4 parameterization the percolation-loss transition is driven by
> global surface-area-mediated thinning, not by failure of the narrowest
> necks. The removed-voxel distribution at the transition step is confined to
> one-voxel-depth surface material (mean 20.00 nm base / 20.02 nm high)
> against a background median of 28.28 nm.

**This explains the null without concluding that lower-tail necks are
irrelevant in all damage models** — and it identifies the specific property
of *this* damage model responsible.

**Branch C (hidden continuous benefit) is NOT supported.** Every continuous
diagnostic (retained volume, spanning-cluster size, spanning fraction of
original) points marginally the *wrong* way, all within noise. There is no
sub-round benefit hiding under the integer outcome.

## Caveat on D2, stated plainly

One of the three requested diagnostics returned no usable information because
the specified proxy measures surface fraction rather than neck thickness at
this voxel resolution. The conclusion above rests on **D1, D3 and D4loc**,
which are mutually consistent. Were a genuine lower-tail thickness diagnostic
wanted, it would need to be restricted to the *neck* population (e.g. SNOW
throat inscribed diameters on the damaged network) rather than an EDT
percentile over all Ni voxels — that is a new metric, not run here.

## Status

Stopping before any Path B write-up or further amendment, as instructed. Not
run: Family C, real-data calibration, opening granulometry, third size
metric, item 6. No D4 parameter altered.
