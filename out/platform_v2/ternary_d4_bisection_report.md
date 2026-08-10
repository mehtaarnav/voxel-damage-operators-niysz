# Platform v2 — YSZ/pore placement, D4 re-validation, base-only bisection

Run 2026-08-10, 91 s. Code: `cmlib/damage.py` (operators promoted from the E0
spike, definitions unchanged), `scripts/platform_v2/ternary_and_d4_pilot.py`.
Data: `ternary_placement.csv`, `d4_revalidation.csv`, `d4_bisection_base.csv`.

## A. YSZ/pore placement — **PASS**

Ni phase untouched (asserted by array equality, not assumed). YSZ share of
the non-Ni remainder targeted at the medium anode's own
0.388/(0.388+0.362) = 0.5173.

| seed | Φ_Ni | Δ vs qualification | Φ_YSZ | Φ_pore | YSZ/rest | P_span before→after | TPB (µm⁻²) |
|---|---|---|---|---|---|---|---|
| 0 | 0.2502 | +3.2e−05 | 0.3879 | 0.3619 | 0.5173 | 1.000 → 1.000 | 20.17 |
| 1 | 0.2493 | +3.4e−05 | 0.3884 | 0.3623 | 0.5173 | 1.000 → 1.000 | 19.98 |
| 2 | 0.2475 | +3.4e−05 | 0.3893 | 0.3632 | 0.5173 | 1.000 → 1.000 | 20.49 |
| 3 | 0.2480 | +3.0e−05 | 0.3891 | 0.3630 | 0.5173 | 1.000 → 1.000 | 20.09 |
| 4 | 0.2562 | +3.2e−05 | 0.3848 | 0.3590 | 0.5173 | 1.000 → 1.000 | 19.43 |

- **Φ_Ni holds at the Platform-v2 target**: mean 0.2502, max |Δ| vs the
  qualification value **3.4e−05**. Holds by construction (placement only
  labels non-Ni voxels) — verified, not assumed, as requested.
- **YSZ share of remainder hits 0.5173 exactly on every seed.** Absolute
  Φ_YSZ ≈ 0.385–0.389 vs the medium anode's 0.388, Φ_pore ≈ 0.359–0.363 vs
  0.362 — close because Φ_Ni is now itself on-target.
- **P_span unchanged to 1e−12 on every seed.**
- **TPB nonzero and structurally plausible: 19.4–20.5 µm⁻².**

**One honest caveat on TPB magnitude (corrected 2026-08-10 — the earlier
"8–10×" understated the low end of the real range).** These values are
roughly **7.3–19.2×** the real-anode range of 1.07–2.65 µm⁻² (`REPORT.md`),
ratio bounds being 19.4/2.65 and 20.5/1.07, depending on which real anode and
which TPB convention is used. That is expected and not a placement
failure: the YSZ/pore field is a smoothed random field with a 3-voxel
correlation length, which produces far more Ni/YSZ/pore triple contacts than
a real sintered microstructure with micron-scale domains. Placement was
explicitly specified as *minimal, not optimised for TPB*. **It is fit for
percolation/damage work now, but TPB density is not yet a quantity to compare
against real data** — that belongs with the deferred P2-C2 size-comparability
obligation in the calibration phase.

## B. D4 re-validation on Platform-v2 geometry — **ALL CHECKS PASS**

E0's intensity grid was **not** assumed to carry over (see §C). Verification
set re-run on the new geometry at n_rounds=5:

| check | result |
|---|---|
| reconstruction integrity (rebuild → bit-identical Ni mask) | 5/5 True |
| YSZ untouched by damage | 5/5 True |
| Ni removed by D4 becomes pore | 5/5 True |
| all metrics finite (no NaN propagation) | 5/5 True |
| largest-component behaviour sensible | yes — see below |

Largest-component behaviour: D4 fragments Ni into ~17,300–17,800 components
before island removal, then keeps the single largest. Post-damage `P_span`
is 1.000 at n_rounds=5 with TPB dropping 20.2 → 16.4 µm⁻². So at this
intensity the network is heavily fragmented in absolute component count yet
still spans — the spanning backbone survives while loose material is
stripped. That is coherent for an erosion-plus-island-removal operator, not
an artifact.

## C. Damage-calibration bisection — BASE structures only

Widened structures were **not inspected**, per instruction. Frozen procedure:
integer `n_rounds`, initial bracket [1,20], expand-only, narrow to width ≤1.
`p_erode=0.35` and `expand_vox=1` unchanged throughout. 5 seeds × 3 damage
seeds = 15 brackets, 6–7 evaluations each.

**Transition-intensity distribution: mean 8.77, sd 0.59, range [7.5, 9.5].**

| bracket midpoint | count |
|---|---|
| 7.5 | 1 |
| 8.5 | 8 |
| 9.5 | 6 |

Every one of the 15 runs bracketed cleanly to width 1 — no bracket expansion
was needed, and no run hit the floor or the [1,20] ceiling.

### Why this matters, and what it says about E0

E0's fixed grid was {2, 5, 6, 7, 8, 10}. The transition on Platform v2 sits
at **8.77 ± 0.59**, i.e. **inside E0's 8→10 gap** — exactly the interval E0
could not resolve, which is why it saw only saturated 1.0/0.0 values. The
decision to replace a shared fixed grid with per-structure bisection is
vindicated on the data: a grid would again have straddled the transition
rather than resolving it.

The spread is also informative for the p10-group experiment to come: **damage
seed matters more than structure seed here** (midpoints vary 8.5/9.5/8.5
across damage seeds *within* a single structure, comparable to the variation
across structures). Any future retention comparison must therefore average
over damage seeds, not treat one as representative.

## Stop condition

Stopping here per instruction. **Not run:** the p10-group damage experiment,
Family C, real-data calibration, opening granulometry, any third image-based
size metric, item 6. Widened structures remain uninspected under damage.
