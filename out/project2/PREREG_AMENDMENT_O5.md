# Pre-registration amendment — O5 (volume-conserving Ni agglomeration)

**Frozen 2026-08-11, before any O5 code is executed.** Amends
`PREREGISTRATION_V2_1.md` (`e62f30b`). All of its anti-tuning rules carry over
unchanged. Motivated by the Step 2 failure (`b278dc5`).

---

## 1. Why

Step 2's decisive mismatch is a magnitude, not an ordering: at the point Ni loses
percolation, O1 has destroyed **99.5 %** of TPB, while the real fine anode
**retains 79.9 %** of its TPB at its worst Ni-percolation retention (0.680) — a
~164× discrepancy. Separately, real Φ_Ni **rises** in the coarse anode (+0.0146),
which no removal-only operator can produce.

**Hypothesis.** Ni percolation loss is driven by **volume-conserving Ni
redistribution** (dewetting / agglomeration by surface diffusion), not by Ni
removal. Redistribution breaks long-range connectivity while preserving the
Ni–YSZ contact perimeter where TPB lives.

**Physical basis.** In the dewetting/agglomeration regime matter flows from high
mean curvature to low: a neck of radius r has mean curvature 1/(2r) against a
particle of radius R with 1/R, so for r < R/2 material leaves the neck and joins
the particle — the Rayleigh-type instability that thins necks and grows bodies.
O5 is a phenomenological voxel implementation of that flow direction, not a
solver.

## 2. Frozen operator definition

**O5, per round:**

1. Ni **surface** = Ni voxels 6-adjacent to non-Ni.
2. Ni **growth front** = pore voxels 6-adjacent to Ni (never YSZ — YSZ is never
   overwritten).
3. Local Ni density is a uniform filter of the Ni mask over a cube of edge
   `W_CURV`; low density at a surface voxel ⇒ locally convex ⇒ high curvature.
4. `K = round(MOVE_FRAC × n_surface)` voxels are **removed** from the surface
   sites of *lowest* local density, and `K` voxels are **added** at growth-front
   sites of *highest* local density. Ties are broken by a seeded random jitter.
5. Total Ni voxel count is conserved exactly except where the growth front is
   smaller than K.

**Frozen parameters — set now, never tuned:**

| parameter | value |
|---|---|
| `MOVE_FRAC` (fraction of Ni surface relocated per round) | **0.05** |
| `W_CURV` (curvature-proxy window, voxels) | **5** |
| connectivity | 6 |
| intensity variable | `n_rounds` (integer) |

**No largest-component pruning**, deliberately and unlike O1/D4. O1's pruning
modelled Ni *loss*; O5 models Ni *redistribution*, and deleting disconnected Ni
would destroy the volume conservation that is the whole point. Disconnected Ni
remains present and is simply not counted by `P_span`.

## 3. Frozen success criteria

**PRIMARY (new, and the real test):**

> **TPB retention at the Ni-percolation transition ≥ 0.50 in at least one analog
> class**, measured against the frozen R2 pristine baseline
> (`step2_r2_tpb_baseline.csv`: 27.708 / 17.228 / 11.294 µm⁻²).

Real is 0.590–0.799. Step 2's O1 gave 0.005. **≥ 0.50 is a qualitative change of
mechanism**, and is deliberately a loose, order-of-magnitude bar — the
ordinal-only constraint forbids fitting to the real value.

**SECONDARY:** C1 ordering — fine loses Ni percolation at strictly lower
intensity than medium and coarse, with ≥ 1.0 round separation.

**Volume-conservation gate (must pass, else the run is void):** |ΔΦ_Ni| / Φ_Ni
≤ 0.5 % at the transition state for every structure.

## 4. Frozen decision rules

| outcome | conclusion |
|---|---|
| TPB retention ≥ 0.50 **and** C1 passes | the agglomeration mechanism explains the divergence — headline result |
| TPB retention ≥ 0.50, C1 fails | mechanism right, ordering driver missing; investigate the size-compression fix (`DESIGN_MEMO` §1.4–1.5, ~400³) |
| **TPB retention < 0.50** | **the agglomeration hypothesis is falsified as implemented.** Report as a second mechanism-specific null. **Do not add operators. Do not tune.** Fall back to reporting Project 2 as a two-null study. |

## 5. n_secondary rule — corrected

Step 2's rule picked the midpoint of the *shallowest transition bracket*, which
for O1 landed **past** the transition (n = 10 vs transitions at ≈ 8.6), so the
secondary saturated at zero and was uninformative.

**Corrected rule for O5:** from a 1-seed scoping sweep, take the bracket
**below** the shallowest observed transition and use its **lower** endpoint,
clipped to [2, 15]. If no transition is observed, use 5. This guarantees a
non-saturating measurement point.

## 6. Required metrics

Per structure × damage seed: transition midpoint and bracket flag; Φ_Ni
conservation error; TPB pristine and at transition; **TPB retention**; Ni–YSZ
interfacial area pristine and at transition; Ni `P_span` / `P_reach` /
`P_largest` / `n_clusters`; retained `P_span` at `n_secondary`. Plus one
before/after visualization per analog class showing neck thinning and body
growth.

## 7. Scope

Same 15 cached Step-1 structures, 3 damage seeds, integer bisection [1, 20],
expand-only, final width ≤ 1, damage-seed averaging mandatory, group differences
< 1.0 round unresolved. `cmlib/damage.py` and `cmlib/synth.py` remain unmodified.
No exact numerical fitting; no multi-physics solver.
