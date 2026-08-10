# C1-real gates — **A3 PASS, A1 FAIL. STOP before bisection.**

Run 2026-08-11 under `PREREG_C1_REAL.md` (frozen `49b1059` before any run).
Code: `scripts/project2/c1real_gates.py`. Data: `c1real_a3_roi_sizing.csv`,
`c1real_a1_o1_validity.csv`. O1 unchanged at `p_erode = 0.35`, `expand_vox = 1`.
**No parameter was adjusted. No bisection was run.**

---

## Verdict

| gate | result |
|---|---|
| **A3** coarse ROI sizing | **PASS** — 12 µm gives **298 nodes** (≥ 150). Frozen at 12 µm. |
| **A1** O1 validity on real voxels | **FAIL** — O1's expansion step *heals* disconnection and multiplies TPB 8–15×. |

**Per A1(e): STOP and report. `p_erode` and `expand_vox` are NOT adjusted.**

---

## 1. A geometry bug in my own v3 audit, found while implementing A3

`SAMPLES` stores `(nx, ny, nz)` and `(vx, vy, vz)`, but the loaded array is
`(z, y, x) = (s[6], s[5], s[4])` with spacing `(s[9], s[8], s[7])`. **I passed
them reversed in the v3 audit.** Corrected extents:

| anode | I used (z,y,x) | **actual** | I used spacing | **actual** |
|---|---|---|---|---|
| fine | 995, 1304, 733 | **733, 1304, 995** | 19.53, 19.53, 20.0 | **20.0, 19.53, 19.53** |
| medium | 960, 1110, 610 | **610, 1110, 960** | 24.41, 24.41, 25.0 | **25.0, 24.41, 24.41** |
| coarse | 744, 1417, 456 | **456, 1417, 744** | 29.14, 29.14, 30.0 | **30.0, 29.14, 29.14** |

**Impact assessment.** Phase 5 and Step 0 are **unaffected** — they read shapes
from the files rather than from `SAMPLES`. The v3 audit's ROIs were still valid
cubes of real data of approximately the intended size (the spacings differ by
≤ 2.4 %, and SNOW uses their geometric mean, which is exactly volume-preserving),
but they were **placed using wrong extents**, so the "3 non-nested ROIs" tiling
was computed against the wrong axis lengths and the ROIs are not where the
pre-registration intended. The intrinsic disorder metrics (CV of pair distance,
coordination sd) are properties of the local network and are unaffected; the
min-cut *fractions* depend on ROI placement and should be regarded as
provisional. **The v3 conclusion — real networks are disordered with small,
variable, non-planar cuts, unlike the lattice — does not depend on placement and
stands.**

## 2. A3 — coarse ROI sizing: PASS, and it fixes v3's REV problem

| side | fits in stack? | Mvoxel | Ni graph nodes |
|---|---|---|---|
| **12 µm** | yes | 67.9 | **298** ✓ |
| 14 µm | **no** | — | — |
| 16 µm | **no** | — | — |

The coarse stack is 456 slices × 29.14 nm = **13.3 µm deep**, so 14 and 16 µm
are geometrically impossible, not merely memory-limited. **12 µm is both the
only option and sufficient**: 298 nodes against the ≥ 150 bar, and against the
**48–72 nodes** that made v3's coarse results uninterpretable. **Coarse is
admissible at 12 µm.** Frozen.

## 3. A1 — O1 validity: FAIL

| anode | ROI | **pristine P_span** | P_span n=1 | vol loss n=1/3/5 | **TPB retention n=1** |
|---|---|---|---|---|---|
| fine | 8 µm | **0.9821** | **1.0000** | 0.031 / 0.113 / 0.232 | **15.24** |
| medium | 8 µm | **0.9713** | **1.0000** | 0.035 / 0.096 / 0.171 | **12.72** |
| coarse | 12 µm | **0.8878** | **1.0000** | 0.105 / 0.174 / 0.252 | **7.71** |

**Volume loss is monotone and sensible (criterion b passes).** Criteria (c) and
(d) fail:

- **P_span INCREASES from pristine to n = 1, in all three anodes** — 0.9821 →
  1.0000, 0.9713 → 1.0000, **0.8878 → 1.0000**. O1's first action makes the
  electrode *better connected than pristine*. The oxidative-expansion step
  (dilation by 1 voxel into pore) bridges the 1.8–11.2 % of Ni that is
  disconnected in the real structure, and the subsequent largest-component
  pruning then reports a perfect network.
- **TPB is multiplied 7.7–15.2×** at n = 1 (fine 4.48 → 68.2 µm⁻²) before
  falling. Dilating Ni into pore manufactures enormous Ni/YSZ/pore contact.

**Correction to my own automated verdict.** The script printed
"A1: PASS". **That verdict is wrong and is superseded by this report.** My
monotonicity check compared only n = 1, 3, 5 — all of which are exactly
1.0000 — and never compared against the **pristine** state, which is where the
violation is. The check as written could not have detected the failure it was
written to detect. The CSV is correct; the verdict line is not.

## 4. Why this was invisible until now, and why it matters

**On the synthetic platform, pristine Ni P_span was exactly 1.0000 by
construction** — gate G1-c *required* it. So the expansion step had nothing to
heal and no headroom to increase P_span, and its TPB inflation was hidden inside
a platform whose TPB was already ~8× too high.

**Real anodes are 1.8–11.2 % disconnected in the pristine state.** That is a
platform-fidelity gap nobody had measured: G1-c enforced perfect pristine
connectivity, which is not a property of real electrodes, and it masked an
operator artifact.

**Consequence for the interpretation of Step 2.** O1's C1 null on the synthetic
platform was obtained with an operator whose expansion step is benign only when
pristine connectivity is perfect. On real structures the same operator is
disqualified before the first bisection step. This does not retroactively
invalidate Step 2's null — the synthetic result is what it is — but it means
**O1 cannot be carried onto real data unchanged**, and the rate-vs-topology
question cannot be answered with it as frozen.

## 5. What I did not do

- Did not adjust `p_erode` or `expand_vox` (A1(e), explicit).
- Did not run the bisection, the Spearman tests, or C3-real.
- Did not re-run the v3 audit with corrected extents (it is a separate
  pre-registered artifact; the correction is recorded here for the advisor to
  rule on).
- Did not silently accept the script's "A1: PASS" line.

## 6. Recommendation

**The ruling needed is on O1's expansion step.** Three options, for decision:

1. **Drop the expansion step for real-data work** (`expand_vox = 0`), making O1
   pure stochastic surface erosion. This is a *definition* change requiring an
   amendment, and it would no longer be Project 1's frozen D4 — the operator
   would need a new name and its own validity gate.
2. **Keep O1 as frozen and accept that pristine P_span rises**, scoring the
   transition from the n = 1 state rather than from pristine. Cheap, but it
   means the reported "damage" begins from an artificially healed structure.
3. **Use a different operator for real data.** The audit (`c41e4fa`) showed real
   min-cuts are 1–2 % of throats; a throat-targeted operator is now defensible
   on real structures in a way it was not on the lattice, where critical edges
   were at chance overlap with narrow necks.

**My recommendation is option 1**, with the new operator explicitly named and
gated, because the expansion step's physical justification (NiO→Ni volume
change) applies to *reoxidation*, whereas the degradation being modelled is Ni
loss under reduction — and because option 2 begins every measurement from a
structure that is more connected than the real pristine one.

**A3's result should be adopted regardless of that ruling:** coarse at 12 µm
with 298 nodes fixes the v3 sub-REV limitation and makes three-class comparison
possible for the first time.
