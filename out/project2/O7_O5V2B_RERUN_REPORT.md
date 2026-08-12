# O5v2 Option B re-run with corrected stencil — the closure does not reproduce

**Status: finding requires adjudication before any downstream work proceeds.**
Scripts: `scripts/project2/o7_o5v2b_rerun.py`,
`scripts/project2/o7_struct6_offset_check.py`,
`scripts/project2/o7_derivation_check.py`.

## What was checked

`o5v2_optionB_report.md` records acceptance rate **exactly 0.000** at n = 1, 3,
5, 8, concludes the structure is returned unchanged so gate (ii) fails, and
**CLOSES the agglomeration route** ("no third implementation"). The stated
physical reason is that in this geometry `nba_min > nbb_max`, so every candidate
pair is rejected at the first comparison.

## Finding 1 — the implementation contradicts its own frozen specification

`PREREG_O5V2_OPTIONB.md` freezes `nb` as "the count of Ni **6-neighbours**" and
the predicate as `dA <= 0 <=> nb(a) <= nb(b)`.

`cmlib/damage2.py:31` sets `STRUCT6 = ndi.generate_binary_structure(3, 1)`, which
has **sum 7 — the centre voxel is included**. Convolving the Ni mask with it
gives

    nb(a) = nN(a) + 1   for a Ni site      (centre contributes 1)
    nb(b) = nN(b)       for a pore site    (centre contributes 0)

so the implemented predicate is `nN(a) + 1 <= nN(b)`, i.e. **`nN(a) < nN(b)`** —
strictly stronger than the frozen `nN(a) <= nN(b)`. Every **area-neutral**
(dA = 0) move is silently rejected. On a reconstructed dumbbell this discards
850-1022 admissible pairs.

The exact identity `dA = 2*(nN(a) - nN(b))` is itself correct and was
re-verified independently (`o7_derivation_check.py`, max error 0 against
brute-force bond enumeration). The defect is the neighbour count fed into it,
not the algebra.

## Finding 2 — the reported 0.000 does not reproduce on any structure

Committed operator, unmodified, `ysz_mask` empty:

| structure | n | acceptance (committed) | acceptance (corrected) |
|---|---|---|---|
| dumbbell R13 neck4 | 1 / 3 / 5 | 0.983 / 0.953 / 0.942 | 1.000 / 0.959 / 0.957 |
| dumbbell R13 neck3 | 1 / 3 / 5 | 0.966 / 0.941 / 0.913 | 0.990 / 0.960 / 0.964 |
| dumbbell R16 neck5 | 1 / 3 / 5 | 0.900 / 0.850 / 0.850 | 1.000 / 0.973 / 0.969 |
| project lattice s300 | 1 / 3 / 5 | 0.997 / 0.990 / 0.986 | 1.000 / 0.999 / 0.999 |
| project lattice s301 | 1 / 3 / 5 | 0.997 / 0.992 / 0.990 | 1.000 / 0.999 / 0.998 |
| project lattice s302 | 1 / 3 / 5 | 0.997 / 0.993 / 0.970 | 1.000 / 0.999 / 0.999 |

Volume conservation is exact (dV = 0) throughout, as originally reported.

The stated physical explanation — "spheres joined by a straight cylinder are
already at a local minimum with respect to single-voxel swaps", `nba_min >
nbb_max` — **does not hold on a straight-cylinder dumbbell**, where the
committed operator accepts 85-98% of proposals.

## Finding 3 — 0.000 is a degenerate return value, not a measurement

Sweeping the YSZ fraction of pore space on the project lattice:

| YSZ frac of pore | pore-front sites | acceptance (committed) |
|---|---|---|
| 0.00 - 0.99 | 42070 - 443 | 0.988 - 0.997 |
| **1.00** | **0** | **0.0000** |

`apply_o5v2b` breaks out of its loop when `fi.size == 0` (no pore front) or
`k_move <= 0`, leaving `proposed = 0`, and returns
`accepted / max(proposed, 1) = 0/1 = 0.0`.

**Acceptance 0.000 is therefore returned both when every move is rejected and
when no move is ever proposed.** The two cases are indistinguishable in the
recorded output, and they have opposite meanings. The report interprets the
value as the former; the latter is also consistent with everything recorded
(unchanged structure, `S_spec(1) = S_spec(0)` exactly, dV = 0).

## What this does NOT establish

- The original O5v2 driver **is not in the repository** — only a transcription of
  its results (`scripts/project2/o5v2_transcribe.py` copies numbers out of
  markdown). Exact re-execution is impossible, so the true cause of the recorded
  0.000 on the original structures is **not** determined here.
- Pristine `S_spec` was 0.45052 in the original run; the structures used here
  give 0.276-0.617. None is the original geometry.
- Nothing here shows the agglomeration hypothesis is true, or that a greedy
  operator can reach Rayleigh break-up. Finding 2 shows moves are *accepted*,
  not that the right physics results.

## Consequence for Failure Mode 6 / 7

The manuscript's area-barrier argument holds that a monotone area-reduction gate
forbids the energy-raising moves the instability requires. That argument is
unaffected as *algebra*. What is now in question is the *empirical* claim resting
on it — "the greedy operator accepts no move" — which is the measurement cited
for closing the route.

The candidate Failure Mode 7 (physical kT is ~4e4 times too small to cross a
voxel-face barrier) remains arithmetically correct and independently verified,
but its stated consequence — that the operator is therefore frozen — is
contradicted: at T\* = 0 the corrected operator accepts 5.5% of moves on a
dumbbell and ~100% on the project lattice, because area-lowering and
area-neutral moves are abundant and require no barrier crossing.

## RESOLUTION — structure recovered, closure OVERTURNED

The original structure was recovered (built inline in a heredoc, never saved):
60^3, two R=10 spheres at z=15 and z=45, 3x3 neck, YSZ slab at x>=55.
`scripts/project2/o7_sequential_greedy.py` reproduces it and confirms
**S_spec = 0.45052 exactly**, with 2326 pore-front sites.

### Finding 3 does not apply to the original run

The original had 2326 pore-front sites and gave `proposed = 1, accepted = 0`.
The recorded 0.000 was a **genuine rejection at t = 0**, not the degenerate
`proposed = 0` path. The degenerate return is a real latent defect and should
still be fixed, but it did not produce the recorded number.

### Finding 1 is the cause, confirmed on the original structure

| stencil | nba_min | nbb_max | `nba_min > nbb_max`? |
|---|---|---|---|
| committed (7-element) | 4 | 3 | **True** -> reject at t = 0 |
| true 6-neighbour (frozen spec) | 3 | 3 | **False** -> first pair admissible (dA = 0) |

The recorded physical explanation is true *as the committed code computes it*
and false *under the frozen specification*. The +1 centre bias is the whole
difference.

### Fourth defect — batch composition

`dA = 2*(nN(a) - nN(b))` is exact for **one isolated move** against a current
neighbour field. `apply_o5v2b` applies ~59 moves per round, ranked once against
a **stale** field. The moves interact and the per-move `dA <= 0` guarantees do
not sum. With the stencil corrected but still batched, acceptance is ~96% and
S_spec **rises** at every intensity (0.46365 / 0.46461 / 0.46604 against
pristine 0.45052).

### Sequential greedy — never previously run — PASSES gate (ii)

One move, field recomputed after every accept. This is the operator the dA
identity actually describes.

| budget | proposed | accepted | rate | S_spec | dS | neck vol | dV |
|---|---|---|---|---|---|---|---|
| 61 | 61 | 61 | 1.000 | 0.44861 | **-0.00191** | 99 | 0 |
| 305 | 305 | 305 | 1.000 | 0.43930 | **-0.01122** | 99 | 0 |
| 1220 | 1081 | 1081 | 1.000 | 0.42235 | **-0.02817** | 99 | 0 |

Strictly decreasing and monotonic, at exact volume conservation. It terminates
at 1081 moves because no admissible pair remains — a true local minimum, after a
6.3% area reduction.

**Gate (ii) requires `S_spec(1) < S_spec(0)` strict. It passes.** Under the
frozen rule "Pass => frozen bisection", the agglomeration route should have
proceeded. **The closure is overturned**, caused by two independent departures
from the frozen specification: the centre-included stencil and batched
application against a stale field.

### What is still NOT established

- **The neck never thins** (volume 99, min cross-section 9, unchanged at every
  budget). Area is reduced by smoothing the spheres, not by neck evolution. No
  Rayleigh break-up is demonstrated. Passing gate (ii) is not evidence that the
  operator reproduces the target physics.
- Gate A1v2 conditions (i), (iii), (iv), (v) were not evaluated here; only (ii),
  the condition that closed the route.
- The neck metric used here (free span z = 25..35, volume 99) does not match the
  originally reported 63, so the neck definitions differ. S_spec, which carries
  the gate, matches exactly.

## Manuscript consequence

The B6 empirical claim **"the greedy operator accepts no move"** is wrong and
must be cut. What actually holds:

> A batched greedy operator accepts most proposals and still fails to reduce
> surface area, because per-move dA guarantees do not compose across a batch
> ranked on a stale neighbour field. A sequential greedy operator reduces area
> monotonically and does not reproduce the recorded failure.

The area-barrier **algebra** is untouched. The sentence resting on it is not.
`o5v2_area_barrier.csv` rows `greedy_area, acceptance_rate = 0.0` are transcribed
from the affected run and must be regenerated, not re-transcribed.

## Gate A1v2, all five conditions, sequential operator

`scripts/project2/o7_gate_a1v2.py`, recovered structure, seeds 300/301/302,
n = 1/3/5, k = 61 moves per round (0.03 x 2026 surface voxels).

| seed | n | prop | acc | rate | S_spec | dPhi | TPB | R_Ni | YSZ ok |
|---|---|---|---|---|---|---|---|---|---|
| 300 | 1 / 3 / 5 | 61 / 183 / 305 | all | 1.000 | 0.45028 / 0.44909 / 0.44837 | 0 | 0 | 0 | True |
| 301 | 1 / 3 / 5 | 61 / 183 / 305 | all | 1.000 | **0.45052** / 0.44933 / 0.44909 | 0 | 0 | 0 | True |
| 302 | 1 / 3 / 5 | 61 / 183 / 305 | all | 1.000 | 0.45004 / 0.44980 / 0.44837 | 0 | 0 | 0 | True |

**Verdict:** (i) PASS, (ii) **FAIL**, (iii) PASS, (iv) PASS, (v) PASS.

(ii) fails on **strictness only**: seed 301 at n = 1 returns S_spec exactly equal
to pristine, because all 61 of its moves were area-NEUTRAL (dA = 0). Seeds 300
and 302 decrease. Monotonicity holds for all three seeds. So the gate outcome is
**seed-dependent at n = 1**, and the failure is equality, not increase.

This corrects the earlier statement in this report that gate (ii) passes: that
came from a deterministic variant without tie-breaking, which is not what the
gate specifies.

### Degeneracy audit — four of five conditions are vacuous here

| condition | value on this structure | informative? |
|---|---|---|
| (i) \|dPhi_Ni\| <= 0.005 | dV = 0 exactly by construction, dPhi = 0 | **No** |
| (ii) S_spec strict decrease | 0.45052 -> 0.45004..0.45052 | Yes, but marginal + seed-dependent |
| (iii) TPB(n) <= TPB(0) | **TPB(0) = 0.000000**; 0 Ni voxels adjacent to YSZ | **No** |
| (iv) R_Ni non-increasing | **R_Ni(0) = 0.0000**; Ni touches neither z face | **No** |
| (v) YSZ untouched | pore front excludes YSZ by construction | **No** |

The YSZ slab is at x >= 55; Ni reaches x ~ 40. The two phases never touch, so the
structure has no triple-phase boundary at all. Ni is also interior to the box, so
face-spanning percolation is identically zero. `PREREG_RNI_METRIC.md`'s mandatory
sanity check (`R_Ni(0)` must equal pristine `P_span`) is satisfied only as
0 == 0.

## Standing conclusion

- The recorded reason for the closure — "the greedy operator accepts no move" —
  **is wrong**. Acceptance is 96-100% once the stencil matches the frozen spec.
- The closure is **not thereby overturned**: gate (ii) still fails as specified,
  on strictness, seed-dependently.
- But the gate that produced the closure **cannot settle the question on this
  structure**, because four of its five conditions are vacuous and the fifth is
  marginal. A gate with one informative, seed-sensitive condition is not a basis
  for closing a route.
- The neck never thins at any budget or seed. No Rayleigh break-up is
  demonstrated by any variant.

## Gate A1v2 on REAL ROIs — the run the pre-registration actually specifies

`scripts/project2/o7_gate_a1v2_real.py`, 1 ROI/anode, n = 1/3/5, AXIS = 2,
CONN = 6, k = 0.03 x surface voxels per round. Operator: `cmlib/seqgreedy.py`,
sequential greedy with incremental neighbour-field updates (validated against
brute force in `o7_seqgreedy_validate.py`: nN field exact, dV = 0, area
non-increasing, YSZ intact).

All conditions are now **informative**: TPB(0) > 0 and R_Ni(0) > 0 for all three
anodes, and the mandatory `R_Ni(0) == P_span` sanity check passes in every case.

| anode | ROI | Mvox | Phi_Ni | S_spec(0) | TPB(0) | R_Ni(0) | k/round |
|---|---|---|---|---|---|---|---|
| fine | 400x410x410 | 67.2 | 0.319302 | 0.15696 | 4.4774 | 0.9821 | 51774 |
| medium | 320x328x328 | 34.4 | — | 0.11670 | — | 0.9710 | 17589 |
| coarse | 400x412x412 | 67.9 | 0.248401 | 0.13991 | 1.5365 | 0.8878 | 35011 |

Outcome at n = 1 / 3 / 5, after the seed defect was fixed and the tie-breaking
policy frozen (seed spread is now real but small; representative seed 300):

| anode | S_spec(0) -> n=5 | dS at n=5 | TPB(0) -> n=5 | R_Ni(0) -> n=5 |
|---|---|---|---|---|
| fine | 0.15696 -> 0.15137 | **-0.00559** | **4.4774 -> 22.51 (5.0x)** | 0.9821 -> **0.9825** |
| medium | 0.11670 -> 0.11419 | **-0.00251** | **1.8660 -> 9.23 (4.9x)** | 0.9713 -> **0.9714** |
| coarse | 0.13991 -> 0.13524 | **-0.00468** | **1.5365 -> 5.69 (3.7x)** | 0.8878 -> **0.8881** |

| condition | fine | medium | coarse |
|---|---|---|---|
| (i) \|dPhi_Ni\| <= 0.005 | PASS | PASS | PASS |
| **(ii) S_spec strict + monotonic** | **PASS** | **PASS** | **PASS** |
| (iii) TPB(n) <= TPB(0) | **FAIL** | **FAIL** | **FAIL** |
| (iv) R_Ni non-increasing | **FAIL** | **FAIL** | **FAIL** |
| (v) YSZ untouched | PASS | PASS | PASS |

Under the earlier LIFO tie-breaking the same gate gave TPB inflation of only
1.5x and (iv) passed on medium. Randomising among equal-dA moves — the frozen,
correct policy — makes the failure far more severe and uniform. The area-neutral
moves are what manufacture TPB: they scatter Ni voxels across the Ni/YSZ/pore
junction at no energy cost.

### What this settles

**Condition (ii) — the condition that closed the agglomeration route — PASSES on
all three real ROIs**, strictly and monotonically, at exact volume conservation.
The recorded closure ("gate (ii) fails because the operator accepts no move") is
wrong on every count: the operator accepts, and it reduces area.

**But the gate still fails, on conditions that are informative here and were
vacuous before:**

- **(iii) TPB is manufactured, not destroyed.** +52% on fine, +34% on coarse.
  This is the manuscript's **second artifact class** (TPB manufacture by
  voxel-scale operations) recurring in a volume-conserving swap operator, not
  only in erosion. Real Ni coarsening reduces TPB; this operator increases it.
- **(iv) R_Ni can rise.** Fine 0.9821 -> 0.9823, coarse 0.8878 -> 0.8880. A swap
  can bridge previously disconnected Ni, so percolation is not monotone under
  the operator. Medium happens not to move.

So the route closure is **correct in outcome and wrong in every stated reason**.
The operator is invalid for coarsening — but because it manufactures TPB and can
increase connectivity, not because it is frozen.

### Defect found in this work, fixed, and its consequence

The first version of `cmlib/seqgreedy.py` accepted a `seed`, constructed
`np.random.default_rng(seed)`, and **never used it** — the same latent defect
flagged in `apply_o5v2b` above. All three seeds returned byte-identical results,
so the gate's "seeds 300/301/302" requirement was not met.

**Fixed.** The tie-breaking policy is now explicit and **frozen** in
`cmlib/seqgreedy.py`:

> Among the valid candidates in the extremal occupancy bucket — i.e. among moves
> with **identical dA** — one is chosen **uniformly at random** from the seeded
> generator. Ties are broken by chance, never by insertion order.

This is part of the operator specification, not an implementation detail,
because area-neutral moves dominate: 208199 of 258870 accepted moves on the fine
ROI at n = 5. On the recovered synthetic structure LIFO reached dS = -0.00024
where randomised reaches -0.00788, a 33x difference from tie-breaking alone.

## RESTATEMENT

Replacing the B6 empirical claim. What is now established, on the three real
ROIs the pre-registration specifies, with a validated operator whose stochastic
policy is frozen:

1. **"The greedy operator accepts no move" is false and must be cut.** Acceptance
   is 96-100%. The recorded 0.000 was produced by a centre-included stencil
   (`damage2.py:31`) that biases Ni sites by +1 and so rejects every area-neutral
   move, contradicting the frozen spec's "count of Ni 6-neighbours".

2. **Gate (ii) — the condition that closed the agglomeration route — PASSES on
   all three real ROIs**, strictly and monotonically, at exact volume
   conservation. The route was closed on a measurement that does not hold.

3. **The route nonetheless does not reopen.** The gate fails on two conditions
   that were vacuous on the synthetic structure and are informative on real
   ROIs:
   - **(iii) TPB is manufactured at 3.7-5.0x**, not destroyed. Real Ni coarsening
     reduces TPB.
   - **(iv) R_Ni rises** on all three anodes: a volume-conserving swap can bridge
     previously disconnected Ni, so percolation is not monotone under the
     operator.

4. **The mechanism of failure is the area-neutral plateau.** Roughly 80% of
   accepted moves have dA = 0. They cost nothing in surface energy and are free
   to scatter Ni voxels across the Ni/YSZ/pore junction, inflating TPB and
   occasionally reconnecting the network. A monotone-area validity gate cannot
   detect this, because area genuinely does decrease.

5. **This is the second artifact class recurring in a new setting.** TPB
   manufacture was previously attributed to voxel-scale *erosion*. It also
   occurs in a strictly volume-conserving *swap* operator that passes its
   surface-area gate. The artifact is not a property of erosion; it is a
   property of voxel-scale moves at a three-phase junction.

6. **The area-barrier algebra is untouched.** `dA = 2*(nN(a) - nN(b))` was
   re-verified exactly, and its three-phase generalisation
   `dE = J_NP*[2*(nN(a)-nN(b)) + (1+cos_t)*(nY(a)-nY(b))]` was derived and
   verified independently (`o7_derivation_check.py`).

7. **No Rayleigh break-up is demonstrated by any variant, seed, or budget.** The
   neck never thins. The mechanism question the milestone set out to answer
   remains open.

### Consequent edits required outside this report

- Cut the B6 sentence asserting zero acceptance.
- `o5v2_area_barrier.csv`: the five `greedy_area, acceptance_rate = 0.0` rows are
  transcribed from the affected run and feed a manuscript figure via
  `scripts/project2/make_figures.py:250`. **Regenerate, do not re-transcribe.**
- `cmlib/damage2.py:31`: the centre-included `STRUCT6` remains **unfixed**.
  Fixing it changes prior committed results and is a scope decision, not a
  correction to be made silently.
