# O6 validity gate — **FAIL. The expansion step was not the cause.**

Run 2026-08-11 under `PREREG_O6.md` (frozen `f980562` before implementation).
O6 = O1 with `expand_vox = 0`. Data: `c1real_o6_validity.csv`. `p_erode`
NOT adjusted. No bisection run. O1 untouched in `cmlib/damage.py`.

## Result

| anode | pristine P_span | P_span n=1 | pristine TPB | TPB n=1 | retention |
|---|---|---|---|---|---|
| fine | 0.9821 | **1.0000** | 4.477 | **66.264** | **14.80** |
| medium | 0.9713 | **1.0000** | 1.866 | **22.902** | 12.27 |
| coarse | 0.8878 | **1.0000** | 1.536 | **11.313** | 7.36 |

| criterion | fine | medium | coarse |
|---|---|---|---|
| (i) `P_span(n) ≤ pristine` | **FAIL** | **FAIL** | **FAIL** |
| (ii) volume loss monotone, sane | PASS | PASS | PASS |
| (iii) `TPB(n) ≤ pristine` | **FAIL** | **FAIL** | **FAIL** |
| (iv) YSZ untouched | PASS | PASS | PASS |

**Removing the expansion step changed almost nothing** (fine TPB ratio 15.24 →
14.80; P_span still exactly 1.0000). **The physics ruling was correct about the
expansion step being inappropriate, but the expansion step was not what caused
the A1 failure.** My earlier attribution was wrong, and so was the diagnosis it
rested on.

## The two actual causes

**1. Largest-component pruning makes criterion (i) unsatisfiable by
construction.** `P_span` is *spanning-cluster voxels / total phase voxels*. Both
O1 and O6 keep only the largest component. Once every non-spanning cluster has
been deleted, the surviving component is the spanning one, so `P_span = 1.0000`
identically — **the pruning step rewrites the denominator.** On the synthetic
platform this was invisible because pristine `P_span` was already exactly
1.0000 (gate G1-c). On any structure with pristine disconnection — i.e. every
real electrode — it is guaranteed.

**This is a metric/operator incompatibility, not an operator defect.** The fix
is not to the operator: either drop island pruning, or score on a
pruning-invariant quantity (e.g. spanning-cluster voxels as a fraction of
*pristine* Ni, which is the physically meaningful "how much conducting Ni is
left" and does not move when isolated material is deleted).

**2. Voxel-scale stochastic erosion manufactures TPB.** At `p_erode = 0.35` the
Ni surface becomes pitted at the single-voxel scale, and every new Ni/pore facet
that touches YSZ creates triple line. TPB rises 7–15× before falling. **This is
the same roughening artifact O5 exhibited** (surface area +3.3 %, TPB ×2).

So criterion (iii) will fail for *any* voxel-scale stochastic erosion operator at
this `p_erode`, independent of expansion or pruning.

## What this means

Three operators (O1, O5, O6) have now failed on real voxels for two recurring,
now-identified reasons. **Neither is about the degradation physics; both are
about how a voxel operator interacts with a real, imperfect microstructure.**

The synthetic platform hid both: G1-c guaranteed pristine `P_span = 1.0000`
(hiding cause 1), and TPB was already ~8× inflated (hiding cause 2).

## Not done

`p_erode` not adjusted (Amendment A, explicit). Bisection, Spearman tests, and
C3-real not run. Amendment D (v3 re-run with corrected extents) and Amendment E
(TPB verification) not reached — **though E is partly answered**: the corrected
estimator gives pristine fine TPB = **4.477 µm⁻²** against a literature ~3.62,
i.e. **1.24×**, not 8×. The estimator is sound on real data; the 8× inflation
was a property of the synthetic YSZ placement, not of the estimator.

## Recommendation

**The next decision is about the outcome metric, not the operator.** Score Ni
percolation loss on a **pruning-invariant** measure — spanning-cluster voxels
divided by *pristine* Ni voxels — which is monotone under erosion, unaffected by
island deletion, and directly comparable to the real retention values already in
hand (0.680 / 0.855 / 0.947). Gate C3-real on TPB **retention ratio** only, per
Amendment E, and treat any operator whose retention exceeds 1.0 as roughening
rather than eroding.
