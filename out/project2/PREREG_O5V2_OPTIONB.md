# Amendment — O5v2 Option B: greedy (zero-temperature) KMC

**Frozen 2026-08-11 before implementation.** Amends `PREREG_O5V2.md`
(`304aa8a`). Last pre-registered implementation path. **If the gate fails, the
agglomeration route CLOSES — no third implementation.**

## Correction adopted
My self-prompt said both "Metropolis" and "accept only area-non-increasing",
which is contradictory: Metropolis accepts ΔA>0 with nonzero probability and so
cannot guarantee gate (ii). **Greedy / zero-temperature KMC is used.**

## Algorithm (frozen)
Propose a swap of one surface Ni voxel `a` → one pore front site `b`.
**Accept iff ΔA ≤ 0 and |ΔV| ≤ 1.** Reject otherwise.

Exact discrete ΔA, derived not assumed: removing `a` changes exposed faces by
`2·nb(a) − 6`; adding at `b` by `6 − 2·nb(b)`; so **ΔA = 2·(nb(a) − nb(b))**,
where `nb` is the count of Ni 6-neighbours. **ΔA ≤ 0 ⟺ nb(a) ≤ nb(b)** — remove
where Ni neighbours are few (convex), add where they are many (concave). ΔV = 0
by construction (one out, one in).

## Frozen constants
| constant | value |
|---|---|
| γ (surface energy weight) | **1.0** |
| λ (volume penalty) | **∞ — enforced structurally**, ΔV = 0 by construction |
| kT | **0** (pure greedy; no ΔA>0 acceptance) |
| `p_coarsen` | **0.03** (unchanged) |
| swap partners | required **non-adjacent** so the ΔA algebra holds exactly |

**Recorded per run:** acceptance rate at n = 1, 3, 5; per-round volume error.

## Gate A1v2 — unchanged, all five conditions
1 ROI/anode, n = 1,3,5, seeds 300/301/302. (i) |ΔΦ_Ni| ≤ 0.005; (ii)
`S_spec(1) < S_spec(0)` strict and monotonic; (iii) `TPB(n) ≤ TPB(0)`;
(iv) `R_Ni` non-increasing; (v) YSZ untouched.

**Fail ⇒ STOP and CLOSE the agglomeration route.** Pass ⇒ frozen bisection
(fine/medium 8 µm ×3, coarse 12 µm ×3, 3 damage seeds, [1,20], thresholds 0.50
and 0.10, raw+partial Spearman). Decision rules unchanged from `PREREG_O5V2.md`.

## Fifth artifact class (confirmed, for the methods paper)
**Curvature-ranked voxel moves do not guarantee surface-area reduction on a
discrete lattice.** Distinct from the four already documented (pruning-dependence
of connectivity metrics; TPB manufacture by voxel-scale erosion;
over-constrained pristine connectivity; lattice min-cut planarity).

## Limitation recorded
Literature verification cannot be performed without source access. Amendment A5
binds: unsourceable claims are **cut, not softened**. O5v2 proceeds independently.
