# O5v2 — validity gate **FAIL at condition (ii). STOP before bisection.**

Under `PREREG_O5V2.md` (frozen `304aa8a` **before** implementation).
`p_coarsen = 0.03` frozen and untouched. No erosion parameter exists in O5v2.
`cmlib/damage.py`, `cmlib/synth.py` untouched; operator in `cmlib/damage2.py`.

## Result on the analytic two-sphere-plus-neck case

Pristine `S_spec` = 0.45052, neck = 63 voxels.

| stencil | n | volume error | `S_spec` | neck voxels |
|---|---|---|---|---|
| 6-conn | 1 | **0.000000** | 0.45195 ✗ | 63 |
| 6-conn | 3 | 0.000000 | 0.45219 ✗ | 63 |
| 6-conn | 5 | 0.000000 | 0.45362 ✗ | 63 |
| **26-conn** | 1 | **0.000000** | **0.45195 ✗** | 57 |
| 26-conn | 3 | 0.000000 | **0.44431** ✓ | 15 |
| 26-conn | 5 | 0.000000 | **0.44216** ✓ | 0 |

**Gate (ii) requires `S_spec(1) < S_spec(0)` strictly. Both stencils fail at
n = 1** (0.45195 > 0.45052). **STOP per the frozen rule; the bisection was not
run and no parameter was adjusted.**

## What the 26-connectivity upgrade did and did not fix

The upgrade was the pre-registered implementation fix for a (ii) failure, and it
**changed the operator's behaviour qualitatively**: with 6-connectivity the neck
never thins (63 → 63 → 63) — the operator does nothing useful. With
26-connectivity the neck thins monotonically **63 → 57 → 15 → 0** and `S_spec`
falls **below** pristine from n = 3 onward.

So **the mechanism direction is correct** — matter does leave convex surface and
fill concave necks, at exactly conserved volume (error 0.000000 at every
intensity). **What fails is the first round only.**

## Why n = 1 increases surface area — diagnosis, not excuse

The move set is chosen by curvature rank, not by the surface-area change it
produces. A single removal from a convex site and a single addition at a concave
site are each *individually* favourable in curvature terms, but on a discrete
lattice the added voxel creates new exposed faces at its own free sides, and at
n = 1 that transient exceeds the area removed. By n = 3 the added material has
coalesced into the neck region and the net area falls. **Curvature rank is a
proxy for the area change; it is not the area change.**

The remaining pre-registered implementation path, **Option B (KMC with a
Metropolis criterion on ΔE = γ(ΔA + λΔV))**, does not have this failure mode by
construction: it accepts a move only if it reduces area, so `S_spec` cannot rise
at any n, including n = 1. **Option B has not been run.**

## Status, stated exactly

- Gate condition (i) volume conservation: **PASS**, error 0.000000 throughout.
- Gate condition (ii) surface-area reduction: **FAIL at n = 1**, both stencils.
- Conditions (iii)–(v) not evaluated on real ROIs, because (ii) is a stop.
- **Bisection not run. C1-real not tested for agglomeration.**
- `p_coarsen` = 0.03 unchanged; `N_move` unchanged.

**The agglomeration hypothesis remains untested, not falsified** — the same
status as after O5. The difference is that the failure is now localised to a
single, understood, and fixable property of the move-selection rule, with the
remedy already specified in the pre-registration.

## Recommendation

Implement Option B (KMC, area-decreasing acceptance) and re-run the identical
gate. It is the last pre-registered implementation path; if it also fails (ii),
the operator class is physically mis-specified and the agglomeration route
closes.
