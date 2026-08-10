# Platform v2 — Ni generator qualification report

Run 2026-08-10. Ni generator only: no YSZ/pore placement, no D4/damage, no
Family C, no real-dataset calibration (all per instruction). Code:
`scripts/platform_v2/{design_probe,qualification_run}.py`, generator in
`cmlib/synth.py`. Data: `qualification_run.csv`,
`qualification_deviations.csv`, `qualification_base_validity_log.csv`,
`qualification_gating_log.txt`, `qualification_run.png`. Design rationale:
`design_memo.md`.

## Headline

**Gates P2-A and P2-B PASS. Gate P2-C FAILS on the c-PSD criterion**, for a
specific, diagnosed, and (in hindsight) predictable reason that is a genuine
property of the lower-Φ_Ni geometry, not a bug and not seed noise.

## Gate P2-A: composition and topology — **PASS**

| criterion | result |
|---|---|
| ≥5 base seeds | 5 ✓ |
| Φ_Ni near 0.250 | mean **0.2502**, range 0.2474–0.2561 ✓ |
| SNOW mean degree in 3.5–4.5 | mean **4.193**, range 4.184–4.206, all 5 in-band ✓ |
| single connected Ni cluster | n_clusters=1 for all seeds ✓ |
| initial P_span intact | 1.000 for all seeds ✓ |
| no severe node-count loss | 97–98 nodes, no loss at base ✓ |

**Topology finding (the open question), stated plainly: plain
nearest-neighbour lattice adjacency already meets the coordination target.
No topology modification was required or built.** Raw topological mean degree
2×224/96 = 4.667; SNOW-measured 4.184–4.206. Degree distribution at base:
median 5, range 1–7, 14/98 nodes (14.3%) at degree 1, **0 nodes at degree 0**.
The earlier "coordination is too low, needs face-diagonal bonds" framing was
an artifact of the ungrounded 6–8 target — confirmed, and no machinery was
built on it.

## Gate P2-B: base neck distribution — **PASS (with one honest note)**

| criterion | result |
|---|---|
| p50/p10 in range / above floor | mean **2.90**, range 2.50–3.00 — all ≥2.5 floor ✓, but mean sits *below* the 3.0–4.3 target |
| base p10 resolved ≥3 vox | mean **5.2 vox** ✓ (exceeds the preferred ≥4) |
| validity/rejection log recorded | ✓ all 5 seeds accepted on first attempt, no rejections |

The p50/p10 shortfall is small and every seed clears the pre-registered
validity floor, but the *mean* (2.90) missing the stated 3.0–4.3 target is
reported rather than rounded up. Seed 4 sits exactly at 2.50, i.e. exactly on
the floor.

## Gate P2-C: lower-tail decoupling — **FAIL**

Per-seed (no aggregate-mean criterion, per `preregistration.md` §0c/B):

**Intermediate (nominal 1.45×) — 0/5 pass:**

| seed | achieved p10 | Φ_Ni dev | c-PSD dev | p50 ratio | n_nodes | P_span | verdict |
|---|---|---|---|---|---|---|---|
| 0 | 1.33 | +0.01% | +0.00% | 1.00 | 98/98 | 1.000 | fail (p10) |
| 1 | 1.33 | +0.02% | **−11.33%** | 1.00 | 97/97 | 1.000 | fail (p10, c-PSD) |
| 2 | 1.33 | +0.01% | +0.00% | 1.00 | 97/97 | 1.000 | fail (p10) |
| 3 | 1.33 | +0.02% | **−10.48%** | 1.00 | 97/97 | 1.000 | fail (p10, c-PSD) |
| 4 | 1.33 | +0.02% | −0.95% | 1.00 | 99/98 | 1.000 | fail (p10) |

**High (nominal 2.0×) — 1/5 pass:**

| seed | achieved p10 | Φ_Ni dev | c-PSD dev | p50 ratio | n_nodes | P_span | verdict |
|---|---|---|---|---|---|---|---|
| 0 | 2.00 | +0.01% | −3.13% | 1.00 | 99/98 | 1.000 | **PASS** |
| 1 | 2.00 | +0.02% | **−12.40%** | 1.00 | 98/97 | 1.000 | fail (c-PSD) |
| 2 | 2.00 | +0.02% | **−9.90%** | 1.00 | 101/97 | 1.000 | fail (c-PSD) |
| 3 | 2.00 | +0.01% | **−10.76%** | 1.00 | 97/97 | 1.000 | fail (c-PSD) |
| 4 | 2.00 | +0.01% | **−9.40%** | 1.00 | 98/98 | 1.000 | fail (c-PSD) |

### Two distinct failure causes, both diagnosed

**(1) Intermediate bin: quantization, not a mechanism failure.** All 5 seeds
landed on the 1.33× rung, below the 1.45× nominal target. Every other
criterion is essentially perfect on those seeds (Φ_Ni ~0.01–0.02%, p50 ratio
exactly 1.00, no node loss, P_span intact). This is the same discrete-rung
behaviour already documented in Family B — achievable p10 values are set by
the base mixture's own value grid, and 1.45× simply falls between rungs.
**Fixable by targeting an achievable rung (1.33×) rather than 1.45×**, since
achieved ratio is the scientific variable anyway.

**(2) High bin: c-PSD deviation is real, systematic, and mechanistic — this
is the item-6 headroom risk materialising.** c-PSD deviation is **negative in
every failing case** (−9.4% to −12.4%), never random-signed. Cause, confirmed
directly from the saved data:

| seed | r_base | r_final | shrink | c-PSD base → high |
|---|---|---|---|---|
| 0 | 12.10 | 11.34 | 94% | 414 → 401 nm |
| 1 | 12.10 | 11.42 | 94% | 452 → 396 nm |
| 2 | 12.10 | 11.28 | 93% | 440 → 396 nm |
| 3 | 12.10 | 11.31 | 93% | 444 → 396 nm |
| 4 | 12.10 | 11.47 | 95% | 444 → 402 nm |

Mass conservation itself is **excellent** (net residuals −31 to +44 voxels on
~30,000 added/removed; Φ_Ni deviation ≤0.02% everywhere). But it is achieved
by shrinking sphere radius 5–7%, and **c-PSD measures particle body size, so
it registers that shrink directly.** At the Family B pilot's Φ_Ni≈0.33
(R=14, spheres tangent) there was enough Ni volume that the same widening
needed proportionally less radius compensation. At Φ_Ni=0.250 the necks carry
~47% of total Ni volume (spheres alone give only 0.132), so widening them
demands a proportionally larger radius give-back.

**This is precisely the headroom concern raised in instruction item 6 — and
it appears not as an r_lo-floor blowout (radii stayed at 93–95% of base,
nowhere near the 1.0-voxel floor, no seed flagged) but as a c-PSD gate
failure.** The mass-conservation mechanism is working exactly as designed;
the tension is that Φ_Ni conservation and c-PSD particle-size conservation
are **not simultaneously satisfiable** at this Φ_Ni when necks carry half the
Ni volume.

### Node counts, per-seed at the high-ratio point (item 6)

No seed lost meaningful nodes — 97–101 vs base 97–98 (seed 2 *gained* 4).
Nothing approached the r_lo floor; no seed was flagged. The aggregate ≥95%
gate would have passed cleanly, which is exactly why the per-seed reporting
was required: **the real failure was invisible to the node-count gate and
showed up only in c-PSD.**

## What this means, and what I have not done

Gate P2-C fails, so **platform v2 is not yet qualified** and I have not
proceeded to YSZ/pore placement or D4 (per instruction, and independently
warranted).

The failure is not "the generator is broken" — Φ_Ni, coordination,
percolation, p50 stability, and node counts are all excellent. It is a
**genuine geometric tension** newly exposed by the (correct) move to
Φ_Ni=0.250: at this composition, necks are ~47% of the Ni volume, so
lower-tail widening cannot be mass-compensated without measurably shrinking
particle bodies, which the c-PSD gate — correctly, since it is the
neck-insensitive size control — refuses.

Options, none of which I have taken unilaterally since each changes a frozen
target or the generator's structure:

1. **Reduce the neck volume fraction** (more spheres/smaller pitch so spheres
   carry more of Φ_Ni, or a narrower base neck distribution), restoring
   compensation headroom. Changes the geometry, needs a re-run of P2-A/B.
2. **Compensate somewhere other than sphere radius** — e.g. remove volume
   from *untouched* wide necks rather than from particle bodies, leaving
   c-PSD (a body-size measure) untouched. Structurally new mechanism.
3. **Re-examine whether c-PSD ±5% is the right ceiling** at Φ_Ni=0.250 given
   that a 5–7% body shrink is arguably the *physically honest* consequence of
   conserving mass while widening necks. This would be relaxing a frozen
   gate and I am not proposing it unilaterally.
4. **Target the 1.33× rung** for the intermediate point regardless (fixes
   failure cause 1 cheaply, independent of the c-PSD question).

Recommendation: option 1 or 2 (keep the frozen gate, fix the generator),
with option 4 folded in either way. Awaiting review before implementing any
of them.

---

## Verification addendum (2026-08-10, post-review) — both checks answered before any option was implemented

### Q2: was "spheres-alone Φ_Ni=0.132" measured or derived? **MEASURED DIRECTLY.**

`design_probe.py` calls
`rasterize(centres, pairs, R_VOX, np.zeros(n_pairs), shape)` — all neck
widths set to zero — and `rasterize` skips every neck via `if w <= 0:
continue`. That is a true spheres-only rasterization, pre-neck-union, not a
subtraction. Φ_Ni(spheres only) = **0.1323** at R=12.1 against 0.2502 with
necks, so necks carry **~47% of total Ni volume**. The c-PSD diagnosis rests
on a directly measured number.

### Q1: is the ±5% c-PSD ceiling grounded in real variability? **NO — real variability is roughly 2× looser.**

c-PSD had **never been computed on the real data** (the real-data study used
watershed sizing and SNOW chamber diameters; `cpsd_r50max` was added later for
the synthetic study). Computed here for the first time, on the same 8 µm ROI
tiling as phase3/4, `sizes=100`
(`scripts/platform_v2/real_cpsd_variability.py`, `real_cpsd_variability.csv`):

| anode | n_ROI | mean | sd | CV | full spread | per-ROI dev from own mean | ROIs beyond ±5% |
|---|---|---|---|---|---|---|---|
| medium | 6 | 1061.3 nm | 106.1 | **10.0%** | 26.9% of mean | −9.7% … +17.3% | **4/6** |
| coarse | 10 | 1264.2 nm | 85.3 | **6.7%** | 23.2% of mean | −8.2% … +15.0% | **5/10** |

Pooled CV across all 16 ROIs: **11.4%**. Anode-to-anode (medium→coarse):
**+19.1%**.

**The ±5% ceiling is tighter than the real material's own ROI-to-ROI
variability.** In the real anodes, 9 of 16 ROIs deviate from their *own
anode's* mean by more than ±5% — i.e. two ROIs cut from the same real anode
routinely differ in c-PSD by more than the synthetic gate permits between a
base structure and its widened counterpart. **Per the pre-registered
instruction, this reopens option 3 as a legitimate fix rather than a
shortcut.**

### Additional finding: the c-PSD measure itself was under-resolved

While verifying, found that `CPSD_SIZES_DEFAULT = 25` quantizes c-PSD into
bins ~6% wide in diameter — comparable to the ±5% gate itself. Only **8
distinct c-PSD values** appeared across all 15 platform-v2 structures.
Convergence test on two seeds:

| sizes | seed 0 dev | seed 1 dev |
|---|---|---|
| 25 | −3.13% (**"pass"**) | −12.40% |
| 50 | −8.05% | −12.21% |
| 100 | −8.04% | −7.57% |
| 200 | −9.18% (**fail**) | −6.37% |

Seed 0's lone P2-C "PASS" was a **quantization artifact**. At converged
resolution every seed sits at roughly −6% to −9%. **The P2-C c-PSD failure is
therefore more uniform than first reported (0/5, not 1/5), while the
individual numbers were less trustworthy than they looked.**
`CPSD_SIZES_DEFAULT` raised 25 → 100.

### Where this leaves the four options

- Real deviation is ~−6 to −9% (converged), real ROI-to-ROI variability is
  CV 6.7–10.0% with 9/16 ROIs beyond ±5%. The synthetic deviation is
  **within the real material's own natural spread**.
- **Option 3 is now evidence-backed**, not a loosening of standards: a
  ceiling near the real CV (~10%) would be grounded in exactly the way the
  ±5% figure never was. This mirrors the earlier coordination-target
  correction, where a round number turned out to be ungrounded once the real
  data was pulled.
- Options 1/2 (change the generator) remain available and would still be
  needed if the goal were to hold c-PSD tighter than the real material does —
  but that is no longer obviously the right goal.
- The cheap fix is applied regardless: `TARGET_RATIOS` 1.45 → **1.33**, the
  rung the generator actually produces (Family B design lesson: achieved
  ratio is the scientific variable, nominal is a label).

**No option has been implemented.** Reporting both findings first, as
instructed.

---

## Convergence verification (2026-08-10) — **FAILED. Nothing locked.**

Per instruction, before locking `CPSD_SIZES_DEFAULT` or any revised ceiling:
(a) confirm the real-ROI and synthetic c-PSD used matching resolution, and
(b) run one more ladder step (sizes=300) to confirm the metric has stopped
moving (<0.5 pp between successive resolutions).

### (a) Resolution mismatch — confirmed, and it was real

The 16 real ROIs used **sizes=100**. The platform-v2 qualification run used
**sizes=25** (the default at the time it executed; the default was raised only
afterwards). The two sides of the comparison in the previous addendum were
**not resolution-matched.** Flagged correctly on review.

### (b) sizes=300 convergence — **the metric has NOT stopped moving**

| sizes | seed 0 dev | step | seed 1 dev | step |
|---|---|---|---|---|
| 25 | −3.13% | — | −12.40% | — |
| 50 | −8.05% | 4.92 pp | −12.21% | 0.18 pp |
| 100 | −8.04% | **0.01 pp** | −7.57% | 4.65 pp |
| 200 | −9.18% | 1.15 pp | −6.37% | 1.19 pp |
| 300 | −9.56% | 0.38 pp | −7.55% | **1.17 pp** |

**There is no resolution at which both seeds change <0.5 pp.** Worse, the
approach is **non-monotone**: each seed produces a spuriously small step
(seed 0: 0.01 pp at 100; seed 1: 0.18 pp at 50) and then jumps by >1 pp at
the next step. Either seed inspected alone would have given a false
"converged" reading — the same trap as the original sizes=25 artifact, one
level up.

Absolute values confirm it is genuine non-monotonicity, not slow decay
(seed 1 base: 452.46 → 453.44 → 453.91 → 454.14 → **458.08**, a +3.94 nm jump
at the last step after settling to +0.23; seed 1 high moves **−1.68 nm**,
reversing direction).

**Diagnosis.** `local_thickness` assigns every voxel a radius from a discrete
grid of `sizes` values spanning the distance-transform range; `cpsd_r50max`
then takes the **median** of that binned field. Changing `sizes` moves the
grid, so the median can hop between bins in either direction. A
median-of-a-binned-distribution has no guarantee of smooth convergence, which
is exactly what is observed.

### Consequence: the ceiling proposal and P2-C re-gate are NOT delivered

Per the instruction ("before locking ... any revised gate threshold"), and
because a threshold cannot be responsibly grounded on a measurement that has
not demonstrated convergence:

- `CPSD_SIZES_DEFAULT` is **left at 100 but explicitly marked unvalidated**
  (it is still far better than 25, which was demonstrably broken — but it is
  not a converged choice).
- **No revised ceiling is proposed**, and **P2-C has not been re-gated.**
- The resolution-matched recompute of all 15 structures was **stopped**
  rather than completed at the unvalidated sizes=100.

### What the numbers do support, stated conservatively

Measurement uncertainty from resolution is roughly **±1.2 pp** (largest
successive step at the top of the ladder). Real ROI-to-ROI spread is CV
6.7–10.0%. So the metric is far too noisy for a ±5% gate (noise is ~1/4 of
the ceiling) but would have ~8× headroom against a ceiling near the real
~10% spread. At sizes=300 the two tested seeds read −9.56% and −7.55%:
**both fail a ±5% ceiling and both pass a ±10% one**, and that verdict does
not flip within ±1.2 pp. This is suggestive, not established — it rests on
2 seeds at one resolution, with a metric that has not converged.

### Recommended next step (not taken)

Rather than pushing `sizes` higher on an unstable estimator, **replace the
gate statistic with one that is not a median of a binned field** — e.g. the
volume-weighted mean of the local-thickness field, or a direct
opening-based granulometry — then re-run this same convergence ladder
against it. A gate is only as trustworthy as the stability of the number it
tests, and that has not yet been demonstrated for `cpsd_r50max`.

---

## Step 2 — candidate-metric correlation pre-check (2026-08-10)

Three size measures on the platform-v2 structures. (b) and (c) computed for
all 15 structures (10 widened points); (a) is the expensive 4-point anchor
reused from the ladder. Data: `cpsd_candidate_precheck.csv`.

| structure set | raw-EDT dev | generator-radius dev | local-thickness dev (sizes=300) |
|---|---|---|---|
| intermediate (1.33×, 5 seeds) | −0.79 to −1.09% | −1.70 to −2.13% | not computed |
| high (2.0×, 5 seeds) | −1.77 to −2.47% | −5.18 to −6.79% | −9.56%, −7.54% (2 seeds) |

**Correlation, raw-EDT vs generator radius (n=10): r = +0.9934.** Direction
and ranking track almost perfectly. **But the scale does not:** offset
+2.44 pp (sd 1.56), i.e. raw EDT reports roughly **one third** of the
body-size change that the generator actually applied.

**Raw-EDT vs true local thickness (2-point anchor):** offsets **+7.38 pp**
and **+5.52 pp** — large and inconsistent. Absolute values differ ~3×
(raw-EDT ≈ 150 nm vs local thickness ≈ 458 nm), because they measure
different things: mean distance-to-surface over all Ni voxels, versus the
radius of the largest sphere *covering* each voxel.

### Verdict: raw EDT is NOT a legitimate stand-in — but the check found something better

**Raw EDT fails the step-2 criterion.** It does not track true local
thickness in magnitude (understates ~4×), and with only 2 anchor points
there is no basis to claim otherwise.

**The more useful finding is the third column.** `generator_radius_deviation`
is *exact* — it is the radius the generator actually applied, not an estimate
— and it is free for every synthetic structure. Against that ground truth:

- **local thickness (sizes=300) overstates** the true body shrink:
  −9.56% / −7.54% vs a true −6.25% / −5.58%.
- **raw EDT understates** it: −2.18% / −2.02% vs the same truth.

So the two image-based candidates bracket the ground truth from opposite
sides, and **neither is faithful**. That is worth knowing before any of them
is adopted as a gate: the incumbent metric was not merely unstable, it was
also biased high against a known answer.

### Step 3 decision — flagged, not silently taken

Per instruction, opening-based granulometry is **not** being adopted as a
free fallback: it samples the same discrete radii that broke
`cpsd_r50max`, so it would need the identical 5-point non-monotone ladder
before it could be trusted, at comparable cost.

Recommendation (not implemented, awaiting review):

1. **Use `generator_radius_deviation` as the synthetic-side gate.** It is
   exact, free, has no resolution parameter, and cannot exhibit the
   binned-median instability. It is the honest answer to "did particle
   bodies shrink, and by how much" for a structure we generated.
2. **Accept that it cannot be the cross-comparison metric** — a real anode
   has no generator radius. The real-vs-synthetic comparison needs its own
   separately-validated measure, or should be reframed.
3. **Worth flagging about the original gate design:** the ±5% ceiling
   compares a *controlled within-structure perturbation* (base → widened,
   same seed) against a ceiling that would be grounded in *between-ROI
   natural variability* of real anodes. Those are different quantities.
   Even with a perfect metric, that comparison may not be the right one.
   Raising this rather than quietly building on it.

**Nothing locked; P2-C still not re-gated.** Steps 4–5 (recompute both
sides at matched resolution, propose ceiling, re-gate) remain blocked on
choosing a metric that demonstrably converges.

---

## Resolution (2026-08-10): path (b) — provisional generator gate, image-based size-comparability deferred

### 1. Provisional internal self-consistency gate — **PASSES on its own terms**

Scope of this gate, stated precisely: it asks *"is the mass-conservation
mechanism behaving sanely?"* — **not** *"are these structures size-comparable
to real tomography?"* (see §2). It uses the exact generator diagnostic, with
the qualification in `design_memo.md` about what that number is and is not
ground truth for.

Per-seed, high ratio (2.0×), R_base = 12.1 vox, r_lo floor = 1.0 vox:

| seed | r_final | radius dev | headroom above r_lo floor | net voxel residual | Φ_Ni dev |
|---|---|---|---|---|---|
| 0 | 11.341 | −6.25% | 93.2% | −15 | +0.01% |
| 1 | 11.423 | −5.58% | 93.9% | +44 | +0.02% |
| 2 | 11.276 | −6.79% | 92.6% | +37 | +0.02% |
| 3 | 11.311 | −6.50% | 92.9% | −31 | +0.01% |
| 4 | 11.471 | −5.17% | 94.4% | −5 | +0.01% |

Intermediate (1.33×): −1.70% to −2.13%, headroom 97.7–98.1%.

**Verdict: the mechanism is not broken.** Radius shrink is a tight, consistent
5.17–6.79% at high ratio (~2% at intermediate); every seed retains >92% of its
headroom above the r_lo floor; Φ_Ni is conserved to ≤0.02%; net voxel
residuals are ±44 on ~30,000 voxels moved. The item-6 headroom concern is
answered — nothing is near the floor, and no seed is an outlier.

### 2. What this does NOT close

P2-C's original purpose was **size-comparability with real tomography** — that
is why an image-based metric was specified rather than a generator parameter.
A generator diagnostic cannot serve that purpose: real anodes have no
generator radius. Two image-based candidates have now failed
(`cpsd_r50max`: non-monotone, no converged resolution; raw EDT: fails to
track local thickness in magnitude, ~3× compressed).

### 3. Explicit choice: path (b)

A reliable *measured* particle-size proxy is not achievable at this domain
size and resolution with the tools in hand. Therefore:

- P2-C's **mass-conservation claim** is gated on the exact generator
  diagnostic above — and **passes**.
- P2-C's **real-data size-comparability claim** is **explicitly deferred to
  the real-dataset calibration phase**, which was out of scope for platform
  v2 in every version of this spec.
- Opening-based granulometry is **not** built. It samples the same discrete
  radii that broke `cpsd_r50max` and would need the identical 5-point
  non-monotone ladder. There is **no consumer for a validated image-based
  size metric until real-data calibration begins**, so building one now is
  effort ahead of need — the same sunk-cost pattern already avoided twice in
  this thread (the ungrounded 6–8 coordination target; the topology machinery
  that proved unnecessary).

**Platform v2 is qualified for its Ni-generator purpose** (P2-A pass, P2-B
pass, P2-C mass-conservation pass), with size-comparability carried forward
as a **named open obligation** of the calibration phase — not a silently
dropped requirement.

### 4. The variance-equivalence argument — written out, and it does not survive

No ceiling is proposed under either path, because the logic a ceiling would
rest on does not hold as stated:

> The ±5% ceiling would bound a **within-structure perturbation** — same seed,
> same geometry, base vs widened, differing only by a deliberate intervention
> — using a bound derived from **between-ROI natural variability** of real
> anodes, where different ROIs are different regions of different material
> with independent microstructural realisations. These are different sources
> of variance, and there is no general reason the first should be bounded by
> the second. Between-ROI spread reflects sampling heterogeneity; within-
> structure deviation reflects the magnitude of a controlled intervention. A
> structure could be perfectly comparable to real material in its size
> distribution while responding to widening by more than the between-ROI CV,
> or vice versa. The defensible use of the real spread is as a
> **measurement-noise reference** — a deviation far below real ROI-to-ROI
> scatter is unlikely to be *resolvable* against real data — not as a *bound*
> on how much a controlled perturbation may change a structure.

Under that corrected reading, the real CV (6.7–10.0%) says the ~6% generator
shrink sits **at the edge of what would be resolvable** against real
between-ROI scatter: useful context, not a pass/fail criterion. Any future
ceiling should be justified from what change is *acceptable* on physical
grounds, not from what is *typical between unrelated ROIs*.
