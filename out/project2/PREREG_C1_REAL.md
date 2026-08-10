# Pre-registration — C1-real: rate vs topology on real Ni-YSZ ROIs

**Frozen 2026-08-11, before any run.** Amends `PREREGISTRATION_V2_1.md`; all
anti-tuning rules carry over. `cmlib/damage.py`, `cmlib/synth.py` frozen. O1 is
used **unchanged** at `p_erode = 0.35`, `expand_vox = 1`.

## Hypothesis

Ni percolation retention is governed by **how fast the damage mechanism reaches
the throats that matter**, not by how few throats must be cut. Fine anodes lose
percolation first because their specific surface area is highest, so
surface-mediated thinning consumes their critical throats soonest.

Motivation (`c41e4fa`): real Ni graphs are disordered and fragile (min-cut
1–2 % of throats), but fine's pristine min-cut fraction (0.0210) is **not**
smaller than medium's (0.0233), while fine measurably loses Ni percolation worst
(retention 0.680 vs 0.855 / 0.947). Pristine topology therefore does not
currently explain C1.

---

## A1 — O1 validity check on real voxels (gates everything)

O1 was designed and validated on the synthetic lattice. Before any bisection,
run O1 at `n_rounds = 1, 3, 5` on **one real ROI per anode class** and verify:

- (b) Ni volume loss is **monotonic** and in a sensible range — not 0 %, not
  ≥ 90 % at n = 5;
- (c) `P_span` decreases **monotonically**;
- (d) TPB is computable on real voxels at every step;
- (e) no pathology at n = 1 (removes nothing / removes everything / instant
  fragmentation).

**If O1 behaves pathologically, STOP and report. `p_erode` and `expand_vox` are
NOT adjusted.** This is a validity check, not a tuning step.

## A3 — Coarse REV sizing (gates coarse inclusion)

Extract coarse ROIs at **12, 14, 16 µm**, record the Ni-graph node count for
each, and **use the smallest size giving ≥ 150 nodes per ROI**. If no size
reaches 150 nodes inside the memory ceiling (~120–150 Mvoxel; coarse voxel
29.14 nm ⇒ 12 µm ≈ 70 Mvox, 14 µm ≈ 111 Mvox, 16 µm ≈ 165 Mvox — the last is
likely infeasible), **report fine-vs-medium only and state the limitation.**

**Coarse results from ROIs with < 150 nodes are not interpreted**, per the v3
finding that 48–72-node ROIs produce 1–2-edge min-cuts and 0.000/1.000 overlap
statistics.

## A4 — Statistical power

**3 ROIs per anode minimum; 5 if memory and time permit.** Report **both**
per-ROI and class-level results. With n = 3 a ≥ 1.0-round separation is a
single-ROI-level claim; with n = 5 it is a class-level claim. Which was achieved
must be stated explicitly.

## A6 — Pre-registered correlation tests (frozen now)

- **Test:** Spearman rank correlation between transition intensity and
  (a) specific Ni surface area, (b) pristine min-cut fraction, **pooled across
  all ROIs**.
- **"Correlates" threshold:** |ρ| ≥ 0.6 **with consistent sign under
  leave-one-out**.
- **If BOTH correlate, report both and do not force a choice.**

## Decision rules

- **PRIMARY C1-real:** fine loses Ni percolation at strictly lower damage
  intensity than medium, **≥ 1.0 round**, on ≥ 3 ROIs per anode.
- **MECHANISM TEST:** A6 above.
- **SECONDARY C3-real (A2):** TPB retention at the Ni-percolation transition,
  **reported for every ROI regardless of whether C1-real passes**.
  - ≥ 0.50 in at least one anode class ⇒ positive C3-real signal.
  - < 0.10 everywhere (as in Step 2's O1) ⇒ confirms the Step 2 TPB failure was
    **operator-specific, not platform-specific**.
- **If C1-real passes and surface area correlates while min-cut does not** ⇒
  rate hypothesis supported; headline.
- **If C1-real passes and min-cut correlates** ⇒ topology governs; revise.
- **If C1-real fails on real structures** ⇒ definitive scope boundary for this
  class of microstructural operator. **Do not build more operators without
  approval.**

## A5 — Publication-claim scoping (binding on all future text)

The claim "first test of a damage operator on real reconstructed Ni-YSZ
microstructures" is **withdrawn**. Prior work (Faes, Grew, and others) has
simulated coarsening on real FIB-SEM data. The defensible claim is:

> "First systematic comparison of surface-erosion-driven percolation loss across
> fine/medium/coarse real Ni-YSZ analogs, testing whether pristine topology or
> damage kinetics controls the ordering."

**This must be verified against the literature before any manuscript draft**;
my prior citations are from memory and unverified in this environment.

## Protocol

Integer bisection [1, 20], expand-only, width ≤ 1, on the **real ROI voxel
mask**, scoring loss of voxel `P_span` under the frozen Step-1 definition
(6-connectivity, free boundaries, spanning-cluster volume fraction, transport
axis = 2). **3 damage seeds per ROI**, damage-seed averaging mandatory. Group
differences < 1.0 round are unresolved, never passes.

**Metrics per ROI × damage seed:** transition midpoint and bracket flag; Ni
volume loss at transition; pristine specific Ni surface area; pristine min-cut
fraction; fraction of pristine min-cut throats destroyed at transition; TPB
pristine and at transition with retention; Ni `P_span` / `P_reach` /
`P_largest` / `n_clusters`.

## Constraints

No tuning of any O1 parameter. No new operators. No multi-physics. Ordinal
comparisons only. `cmlib/damage.py`, `cmlib/synth.py` unmodified.
