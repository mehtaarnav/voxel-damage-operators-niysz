# T5 coupling decision experiment — report

Run 2026-08-10, after gates G0 and G−1. Answers execution-spec Q1/Q2. Full
data: `t5_coupling_experiment.csv` (per-structure), `t5_coupling_experiment_agg.csv`
(aggregated), `t5_coupling_experiment.png` (figure). Code:
`scripts/next/t5_coupling_experiment.py`.

## Setup

5×5×5 cubic lattice of Ni spheres (R=8 vox, pitch=20 vox → base gap 4 vox,
domain 81×128×128 voxels, 20 nm isotropic voxels), 300 nearest-neighbour pairs
each given a random BASE neck width (uniform integer 2–6 voxels). Two
widening modes tested, target neck width ∈ {8,10,12,14} voxels, 5 seeds each,
41 structures + 5 base replicates:

- **selective:** the narrowest ~20% of necks (by base width) set to the
  target width; the rest untouched.
- **uniform:** every neck's width shifted by the same increment.

All measurements (`ws_d_volweighted_nm`, `d_cPSD_r50max_nm`, `neck_p10_nm`,
`neck_p50_nm`) come from the SAME `cmlib.api` functions used everywhere else
in this project (`compute_particle_stats`, `extract_network` +
`compute_network_metrics`) — this is what the analysis pipeline actually
sees, not the generator's "intended" geometry.

## Headline result: NOT satisfiable as built — two specific, fixable reasons found

**Zero of the 8 tested points (either mode) satisfy neck p10 ≥1.5× AND any
size measure within ±5% AND Ni volume fraction within ±5%.** But the failure
has two distinct, well-understood, and separately fixable causes — this is
not evidence the decoupling is infeasible.

### Finding 1 — Φ_Ni is not preserved (dominant effect)

Widening a neck by simply adding Ni voxels adds Ni **mass**. Summed over 300
pairs, this is not negligible:

| mode | target_w | neck p10 (measured) | watershed size dev. | c-PSD size dev. | **Φ_Ni deviation** |
|---|---|---|---|---|---|
| selective | 8 vox | 2.0× | +4.0% | 0.0% | **+11.3%** |
| selective | 14 vox | 2.0× | +18.0% | +3.0% | **+57.4%** |
| uniform | 8 vox | 3.0× | +10.2% | 0.0% | **+30.8%** |
| uniform | 14 vox | 6.0× | +34.3% | +10.7% | **+128.1%** |

Even the mildest tested point overshoots the ±5% Φ_Ni tolerance more than
2×. This is exactly the risk the original Phase-1 "generator requirements"
anticipated ("adjust background phases to preserve volume fractions") — this
scaffold deliberately did not implement that compensation, since it is
generator machinery, not decision-experiment scaffolding. **Fix (not yet
built or tested): shrink Ni particle radius slightly as necks are added, so
added neck volume is offset by reduced particle volume, holding Φ_Ni fixed.**

**Confirms part of the original R1 rationale directly:** at matched neck p10
movement, the **c-PSD size measure is consistently far more stable than the
watershed measure** (0.0–10.7% vs 4.0–34.3% deviation across all 8 points) —
watershed inflation under neck-widening is real and c-PSD is the better
neck-insensitive companion measure, exactly as `cmlib/particles.py` argues.

### Finding 2 — "selective bottom-20%" cannot move p10 by construction (a percentile-definition trap, found empirically)

The **measured** neck p10 for "selective" is **stuck at exactly 2.0× across
all four target widths** (80 nm regardless of whether the widened necks are
set to 160, 200, 240, or 280 nm):

| target_w (intended) | intended-array p10 | **measured (SNOW) p10** |
|---|---|---|
| 8/10/12/14 vox | 76 nm (constant) | 80 nm (constant) |

This is not a bug in the measurement — it is a mathematical property of the
construction. Moving the bottom 20% of a distribution UP AND OUT of the low
end means the new 10th-percentile rank is filled by material that was
originally at the (10+20)=30th percentile of the UNTOUCHED population — the
target width you assign to the widened 20% never enters the p10 calculation
at all, for ANY target value. ("uniform" mode does not have this problem: it
shifts the entire distribution including the lower tail, and its measured p10
tracks the applied widening cleanly — 3×, 4×, 5×, 6× for target_w=8,10,12,14.)

**Fix (not yet built or tested): to move a target percentile P, the widened
population fraction must be ≥P (not an arbitrary larger fraction chosen for
other reasons), or the widening should target the percentile directly by
construction** (e.g., iteratively select exactly enough of the narrowest necks
to guarantee the intended p10 shift) rather than a fixed "bottom X%" rule.

## What this experiment does NOT show

- It does not show decoupling is impossible. "Uniform" mode moves the
  measured neck p10 cleanly (3×–6×); the barrier is Φ_Ni compensation and the
  selective-mode percentile-targeting bug, both diagnosed with concrete fixes,
  not fundamental limits.
- It does not test the FIXED versions of either mechanism (Φ_Ni compensation,
  correct percentile targeting) — those require generator code, which is
  explicitly out of scope for this decision experiment.
- The lattice here (regular cubic, 6 necks/particle, 300 total necks in a
  small domain) may itself be more neck-volume-dense than a realistic random
  packing (typical coordination number ~6–8 is similar, but a REAL packing's
  neck count relative to particle count, and hence the Φ_Ni sensitivity to
  widening, has not been measured and may differ).

## Recommendation (for your decision — no code written past this point)

Per the pre-registration (`out/next/preregistration.md` §8) and execution-spec
Q1/Q2, this is neither a clean "satisfiable" nor "not testable" — it is
**"testable, contingent on two specific, already-anticipated generator fixes
that are not yet built or verified."** Options:

1. **Build the two fixes into `cmlib/synth.py`** (Φ_Ni compensation via
   particle-radius shrinkage; percentile-correct selective widening) and
   re-run this SAME T5-style coupling check before committing to the full
   Family B design — this is the path the original execution spec anticipated
   and is the one I'd recommend, since both fixes are well-understood and
   small.
2. **Relax the Φ_Ni tolerance** (documented deviation, per preregistration
   discipline — not silently) if ±5% proves unreachable even with
   compensation.
3. **Scope the primary claim to the c-PSD size measure only** (finding 1
   already shows c-PSD is far more robust to uncompensated widening than
   watershed) and treat watershed as a secondary, expected-to-move measure —
   consistent with `out/next/EXECUTION_SPEC.md` Q2's recommendation (a)+(b).
4. **Conclude Path B now** (decoupling not practically testable in this
   framework) — not recommended given neither finding is a hard blocker, both
   have named fixes, and "uniform" mode already demonstrates p10 CAN be moved
   substantially through simple, correctly-targeted geometric construction.

No generator, damage model, Family C, or large experiment has been built or
run. Awaiting your decision on the above before any further implementation.
