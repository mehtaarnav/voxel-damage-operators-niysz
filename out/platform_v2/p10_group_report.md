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
