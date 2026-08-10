# Platform v2 — p10-group damage-bisection experiment

Run 2026-08-10, 657 s. Code: `scripts/platform_v2/p10_group_experiment.py`.
Data: `p10_group_experiment.csv`. Design frozen in `preregistration.md`
§0c–§0g; nothing was chosen after seeing results.

15 structures (5 structure seeds × {base, 1.33×, 2.00× achieved p10}) × **5
damage seeds** (200–204, independent of structure seeds) = 75 bisections,
each narrowed to a bracket of width 1. `p_erode=0.35`, `expand_vox=1` frozen
throughout; only `n_rounds` varied.

## Result: **no resolvable effect**

| group | achieved p10 | group mean transition | spread of structure means | n structures |
|---|---|---|---|---|
| base | 1.00× | **8.54** | 0.089 | 5 |
| lower-tail | 1.33× | **8.50** | 0.000 | 5 |
| lower-tail | 2.00× | **8.50** | 0.000 | 5 |

Differences vs base: **−0.040 rounds at 1.33×, −0.040 rounds at 2.00×.**
Both are far below the pre-registered 1.0-round interpretability threshold
(§0g/1), and both are *negative* — i.e. not even a sub-threshold hint of the
hypothesised benefit.

**Per-structure detail:** 74 of 75 bisections returned a midpoint of exactly
8.5. The single exception is base structure seed 0, where one damage seed
returned 9.5 (giving that structure a mean of 8.7 ± 0.447). Every widened
structure — all 10, across both p10 groups and all 5 damage seeds — returned
8.5 with **zero** variance.

### Pre-registered interpretation

This is the **"no effect"** branch. Per the E0 rules and §0f/4, that means
**inconclusive on the scientific hypothesis** — *not* Path B, *not* a
negative result, *not* a claim that neck engineering does not work. The
experiment was already run at the maximum 5 damage seeds, so the §0g/1
"complete to 5 seeds before interpreting" remedy is exhausted; there is no
further seed-completion available to rescue resolution.

### §0f/4 confound note

Moot in the event, but recorded as required: the 2.00× structures carry a
5.17–6.79% sphere-radius shrink (1.70–2.13% at 1.33×). Since there is **no
effect to attribute**, the radius-shrink confound does not arise here — and
the matched-shrink control is correspondingly **not** triggered. Had an
effect appeared, it could not have been attributed to neck widening alone.

## An inconsistency worth flagging: damage-seed variance did not reproduce

§0g/1 was written because the base-only bisection showed within-structure
damage-seed variance comparable to across-structure variance (midpoints
8.5/9.5/8.5 *within* one structure). **That variance did not reproduce with a
different, larger damage-seed set:**

| run | damage seeds | n | midpoint distribution | mean | sd |
|---|---|---|---|---|---|
| base-only (§C) | 100–102 | 15 | 7.5×1, 8.5×9, 9.5×5 | 8.767 | **0.594** |
| p10-group, base structures | 200–204 | 25 | 8.5×24, 9.5×1 | 8.540 | **0.200** |

Same structures, same operator, same frozen parameters — only the damage
seeds differ. The earlier sd of 0.594 was **driven by the particular seeds
100–102**, and a larger independent set gives 0.200.

Two honest consequences:

1. **The §0g/1 requirement was still correct to impose** — it was cheap, it
   is now empirically checked rather than assumed, and requiring ≥3 seeds is
   exactly what exposed this. But the *motivation* for it (large
   within-structure variance) is weaker than the base-only run suggested.
2. **More importantly, it sharpens the null.** With sd ≈ 0.2 rounds and 5
   seeds per structure, the standard error on a group mean is ≈0.09 rounds.
   The observed group difference is 0.04 rounds — **smaller than one standard
   error, and ~25× below the 1.0-round threshold.** This is not a
   marginal miss that more seeds would resolve; the effect is absent at the
   resolution this design can reach.

## Why the transition is so degenerate, and what that implies

The bisection returns integer-bracketed midpoints, so the finest resolution
attainable is 1.0 round. Essentially every structure fails between
`n_rounds` 8 and 9. Neck widening at 1.33× and 2.00× **does not move the
structure across that integer boundary at all**.

That is a resolution ceiling of the outcome variable, not obviously a
statement about the physics: an effect smaller than one erosion round is
invisible to this design by construction. A finer-grained outcome (e.g.
fractional damage intensity via `p_erode`, or retained-conductance rather
than binary percolation loss) would be needed to detect sub-round effects —
but changing `p_erode` defines a different damage model and requires an
amendment (§0g/2), and is **not** taken here.

## Status

- p10-group experiment: **complete, null.**
- Damage-parameter sensitivity check (§0g/2): **not triggered** — that
  obligation applies to a positive or weak-positive branch. This is neither.
- Not run: Family C, real-data calibration, opening granulometry, any third
  image-based size metric, item 6.

---

## Provenance verification of the widened structures (2026-08-10, post-review)

Requested before accepting the null: confirm the "widened" structures fed
into the damage run are the previously-validated ones, not an accidental
reuse of base geometry. This was a live risk — `p10_group_experiment.py`
rebuilds each structure from `intended_T_vox` looked up in
`qualification_run.csv`, so a lookup or seed-indexing slip would silently
produce base geometry while still labelling the row `lower_tail`.

### Check 1 — thresholds and recorded metrics (cross-referenced, not derived)

| mode | nominal | `intended_T_vox` | achieved p10 | neck p10 (nm) | neck p50 (nm) |
|---|---|---|---|---|---|
| base | 1.00× | NaN (no widening) | 1.000 | 120.0 | 320.0 |
| lower_tail | 1.45× | **8.5** | 1.333 | **160.0** | 320.0 |
| lower_tail | 2.00× | **11.0** | 2.000 | **240.0** | 320.0 |

Thresholds are non-null and distinct (8.5 vs 11.0). Neck p10 genuinely moves
120 → 160 → 240 nm while p50 stays pinned at 320 nm — the tail-selective
signature, and exactly the values recorded at qualification. The damage run's
`achieved_p10_ratio` column carries `nunique = 1` per group at 1.000 / 1.333 /
2.000, so no base row leaked into a widened group.

### Check 2 — raw arrays, structure seed 0 (not derived metrics)

**Neck-width arrays:**

| | min | p10 | p50 | max |
|---|---|---|---|---|
| base | 4 | 5.0 | 15.0 | 20 |
| 2.00× | **11** | **11.0** | 15.0 | 20 |

`np.array_equal(base, high)` → **False**. **46 of 224 necks changed** — i.e.
the max-clip raised every neck below 11 voxels, leaving the rest untouched,
which is the intended lower-tail rule. Note p50 is identical (15.0) while the
minimum jumps 4 → 11: the intervention is confined to the tail.

**Ni mask arrays:**

| quantity | value |
|---|---|
| base voxels | 1,137,116 |
| 2.00× voxels | 1,137,101 |
| `np.array_equal` | **False** |
| voxels differing | **67,803** (1.49% of domain) |
| in base only / in 2.00× only | 33,909 / 33,894 |
| `r_final` applied | 11.3412 (vs base R = 12.1) |

The two masks differ in 67,803 voxels, with the base-only and high-only
counts nearly balanced (33,909 vs 33,894) — the mass-conservation signature:
~34k voxels removed from particle bodies, ~34k added at the necks, net −15.

### Verdict

**Both checks pass. The widened structures are genuine and distinct from
base.** There is no path-, seed-, or threshold-indexing slip. The
zero-variance pattern is therefore a real property of this design's damage
response, not a pipeline artifact, and **the null stands as reported**.

Worth stating plainly what this makes the result: two structures that differ
in 67,803 voxels, with a 2× difference in measured neck p10 and a 6.25%
difference in particle radius, fail percolation at *identical* damage
intensity (8.5 rounds) under all five damage seeds. That is a substantive
finding about the damage response, not an absence of signal in the pipeline.
