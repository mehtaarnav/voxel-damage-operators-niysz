# O5v2 Option B (greedy KMC) — **gate FAIL. Agglomeration route CLOSED.**

Under `PREREG_O5V2_OPTIONB.md` (frozen `c414e63` **before** implementation).
γ = 1.0, λ structural, kT = 0, `p_coarsen` = 0.03 — none adjusted.
`cmlib/damage.py`, `cmlib/synth.py` untouched.

## Result

| n | volume error | ΔV (vox) | S_spec | acceptance rate | neck |
|---|---|---|---|---|---|
| 1 | 0.000000 | 0 | 0.45052 | **0.000** | 63 |
| 3 | 0.000000 | 0 | 0.45052 | **0.000** | 63 |
| 5 | 0.000000 | 0 | 0.45052 | **0.000** | 63 |
| 8 | 0.000000 | 0 | 0.45052 | **0.000** | 63 |

Pristine S_spec = 0.45052. **No move was ever accepted.** The structure is
returned unchanged, so `S_spec(1) = S_spec(0)` exactly and gate (ii), which
requires a **strict** decrease, fails.

**Per the frozen rule: the agglomeration route is CLOSED. No third
implementation.**

## Why acceptance is exactly zero — and it is not a bug

The exact identity **ΔA = 2·(nb(a) − nb(b))** means a swap is accepted only if
some surface voxel has **no more** Ni 6-neighbours than some available pore
front site has. In this geometry the *least* convex surface voxel still has more
Ni neighbours than the *most* concave pore site — so `nba_min > nbb_max` and
every candidate pair is rejected at the first comparison.

Physically: **spheres joined by a straight cylinder are already at a local
minimum with respect to single-voxel swaps.** Every single-voxel move increases
area. That is correct for the energy functional, and it is why the greedy rule
does nothing.

## The deeper reason, and the sixth artifact class

**Rayleigh-type neck thinning is a collective, long-wavelength instability.** A
cylinder is unstable to perturbations of wavelength greater than its
circumference, but *any single-voxel perturbation raises the area* — the
instability only lowers energy once a correlated set of voxels has moved
together. **A greedy single-voxel algorithm cannot reach it: it is blocked by a
barrier it is forbidden to cross.**

This is the exact trade-off the two options bracket, and neither side is viable
as specified:

- **Option A** (curvature-ranked, no ΔA test) *moves* material and thins the
  neck 63 → 0, but increases area at n = 1 — it crosses the barrier without
  checking the cost.
- **Option B** (greedy ΔA ≤ 0) *never* increases area, but is frozen at the
  barrier and never moves at all.

**A finite-temperature method is what bridges them — and the amendment
explicitly forbids it** (kT = 0, ΔA > 0 acceptance probability < 1e-6). That
prohibition was adopted to guarantee gate (ii); it also guarantees the operator
cannot exhibit the instability it was built to model.

**Sixth artifact class, for the methods paper:**
> *Greedy surface-area minimisation cannot reproduce Rayleigh-type neck
> break-up, because the instability requires a collective move through a
> positive-area barrier. Enforcing monotone area reduction — the natural way to
> validate a coarsening operator — structurally excludes the mechanism it is
> meant to validate.*

This joins: (1) pruning-dependence of connectivity metrics; (2) TPB manufacture
by voxel-scale erosion; (3) over-constrained pristine connectivity; (4) lattice
min-cut planarity; (5) curvature-ranked moves do not guarantee area reduction.

## Status, stated exactly

- Volume conservation: **exact** (ΔV = 0, error 0.000000) — both options met the
  hard constraint.
- Gate (ii): **FAIL** — Option A rises at n = 1; Option B never decreases.
- Conditions (iii)–(v) not evaluated; (ii) is a stop.
- **Bisection not run. C1-real not tested for agglomeration.**
- No parameter adjusted at any point.

**The agglomeration hypothesis is CLOSED as untestable by this operator class on
a voxel lattice — not falsified.** The distinction matters: nothing here says
agglomeration does not drive real degradation. It says a discrete-lattice
voxel-swap operator, under a validity gate that forbids transient area increase,
cannot represent it.

## Limitation recorded

Literature verification could not be performed — no source access in this
environment. Amendment A5 binds: unsourceable claims are **cut, not softened**.
