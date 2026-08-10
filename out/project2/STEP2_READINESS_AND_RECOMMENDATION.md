# Step 2 readiness assessment and recommendation

Written 2026-08-11 in response to the request for a scientific recommendation on
the next high-leverage move. **No operator was implemented and no damage was
run.** The numbers below come from three cheap measurements on the already-
qualified structures (`scripts/project2/` probes; Ni surface area, TPB, and the
YSZ margin-to-threshold read off `killtest_scope.csv`).

**Bottom line: run Step 2 as planned, but add two pre-checks and one control arm
first. Do not start with O4.** Reasons below, in order of how much they change
the decision.

---

## 1. Three findings that change how Step 2 should be read

### 1.1 C1 is a genuine prediction, but its margin is thin

Under a surface-erosion operator, Project 1 established by direct measurement
that collapse is **surface-area-mediated** (the failure-step voxels were
one-voxel-deep surface material). Specific Ni surface area on the qualified
structures:

| analog | Ni surface / Ni volume | Φ_Ni |
|---|---|---|
| fine | **14.25 µm⁻¹** | 0.3207 |
| medium | 13.69 | 0.2500 |
| coarse | 11.31 | 0.2287 |

Fine is highest, so **C1 (fine loses Ni percolation first) is predicted by the
geometry, not tuned into it.** That is a real, falsifiable prediction and the
strongest scientific asset going into Step 2.

**But the spread is only 1.26×**, against the ~1.49× implied by the real
particle-size ratio — compressed by the 3× size compression and the
proportionally fat necks already recorded in `DESIGN_MEMO` §1.4. Project 1's
integer bisection resolves nothing below **1.0 damage round**, and that
resolution ceiling is exactly what produced its null. **A 1.26× driver may not
move the transition by a full round.** This is the single largest risk to
Step 2, and it is a resolution risk, not a physics risk.

### 1.2 TPB now matches real in ordering *and* ratio — C3 is scoreable

Measured on the qualified structures (particulate YSZ, seed 0):

| analog | synthetic TPB | real TPB | ratio |
|---|---|---|---|
| fine | 27.80 µm⁻² | 3.624 | 7.7× |
| medium | 17.53 | 2.109 | 8.3× |
| coarse | 11.75 | 1.473 | 8.0× |

Absolute magnitude is still ~8× high — the same order as the random field's
7.3–19.2× — **but two things improved that matter for C3**: the over-prediction
is now *uniform* across analogs (7.7–8.3× rather than 7.3–19.2×), and the
fine/coarse ratio is **2.37 against a real 2.46**. Ordering and ratio both match.

**C3 is therefore scoreable ordinally**, which was not safe to assume before.
This had not been measured on the particulate structures at all — the TPB caveat
on record was for the random-field placement.

### 1.3 C2's margin is real but the three-level version is not resolvable

Margin between the operating `p_sinter` and each analog's measured spanning
onset (from `killtest_scope.csv`):

| analog | p_operating | spanning onset | margin |
|---|---|---|---|
| fine | 0.955 | 0.15–0.25 | **~0.76** |
| medium | 0.523 | 0.15–0.25 | ~0.33 |
| coarse | 0.416 | below 0.15 | **~0.29** |

**C2 as specified (coarse before fine) has a ~2.6× margin and should resolve.**
The three-level version (coarse < medium < fine) has medium and coarse only
~16 % apart and will probably not — which is acceptable, since the real data's
medium-vs-coarse is itself unresolved (Project 1, carried caveat).

**The pre-loading concern, stated plainly.** `p_sinter` was calibrated to match
real *pristine* Q. One could argue C2 is then built in: the analog that starts
closest to threshold fails first. I think that is **partly** true and must be
controlled, not waved away — but it is not a tautology, because the calibration
targeted pristine Q while the margin-to-threshold is a different quantity that
emerged from it (and note medium, not coarse, has the highest pristine
coordination-normalised fragility, so the ordering was not forced). **The clean
way to separate them is §3.**

---

## 2. Recommendation, ranked

### R1 (do first, ~1 hour) — pre-register the C1 resolution risk and fix it if cheap

Before running Step 2, decide how a sub-round C1 effect will be handled. Project
1 spent its entire budget discovering that an integer-bisection outcome cannot
see effects below one round. Two options:

- **preferred:** add a continuous secondary outcome alongside the bisection —
  retained Ni `P_span` at a *fixed* damage intensity, rather than the intensity
  at which spanning is lost. This is free (the damage run already computes it),
  has no resolution floor, and would have detected a sub-round effect in
  Project 1.
- record it as a pre-registered secondary, **never** as a substitute for the
  frozen primary.

This is the highest-value hour available, because it insures against the exact
failure mode that ended Project 1.

### R2 (do first, already mostly done) — TPB baseline

§1.2 is the baseline. Freeze these pristine TPB values now so the C3 comparison
at the transition state has a pre-recorded reference.

### R3 (add to Step 2) — a matched-`p_sinter` control arm for O3

Run O3 on a fourth structure set: **three analogs at identical `p_sinter`**
(suggest the medium value, 0.523), so grain size varies while pristine YSZ
fragility does not.

- If coarse-worst **persists** at matched pristine fragility, C2 is mechanistic
  and the result is strong.
- If coarse-worst **vanishes**, C2 was pristine-state loading, and the honest
  claim shrinks to "YSZ fragility ordering is an input, and the operator
  preserves it."

Cost: 5 seeds × 3 analogs at one extra `p_sinter` ≈ the same as one Step 2 arm.
**Without this control, a C2 pass is not safely attributable to the mechanism.**
I regard this as the difference between a publishable causal claim and an
over-claim.

### R4 — treat C3 as the discriminating test, not the tie-breaker

C1 is predicted by surface area; C2 is partly loaded by construction. **C3 (TPB
divergence) is the only criterion that is neither**, and it is the one the field
has no mechanistic explanation for: fine retains TPB *best* (0.799) while losing
Ni percolation *worst* (0.680). If Step 2 reproduces that, it is the most novel
result available here — more so, in my judgement, than the ordering results.

---

## 3. On O4 (the coupled mechanism) — worth doing, but reframe the claim

**The physics is sound.** NiO → Ni is a ~70 % volume change (molar volumes
≈ 11.2 vs 6.59 cm³/mol), and reoxidation-driven expansion cracking the YSZ
scaffold is the accepted redox-degradation pathway in the Ni-YSZ literature
(Sarantaridis & Atkinson's redox-cycling review is the standard entry point;
Simwonis et al. for Ni coarsening). *These citations are from memory and must be
verified from source before they appear in a manuscript — I have no offline
access to them in this session.* O4 also fits our architecture cleanly, because
we have an explicit YSZ contact graph to sever, and a smoothed local-Ni-fraction
field is a defensible strain proxy under the no-multi-physics constraint.

**But the headline as drafted would over-claim.** "A single mechanochemical
coupling explains the divergent degradation" implies the coupling produces both
halves of the divergence. It cannot, on our own numbers: strain ∝ local Ni
content, and **fine has the most Ni** (Φ = 0.322 vs 0.229). So the coupling
alone predicts the *most* YSZ damage in fine — the opposite of C2. Any C2 pass
under O4 will be carried by the coarseness-dependent YSZ contact strength, which
is a generator **input**, not an output of the coupling.

**The defensible claim** is therefore: *"the divergence arises from a single
Ni-driven strain field acting on a coarseness-dependent YSZ contact strength"* —
one driver, two susceptibilities. That is still novel and still connects the two
literatures, but it is honest about which half is assumed. With R3's matched
control, the assumed half becomes measurable rather than asserted.

**Sequencing:** O4 after Step 2, not instead of it. O4's interpretation depends
entirely on knowing what the independent operators do first — that is what makes
"the coupled mechanism is a stronger explanation" a comparison rather than an
assertion.

---

## 4. What I would *not* do

- **Do not start with O4.** Without O1/O2/O3 baselines there is nothing to
  compare against, and its C2 result would be uninterpretable.
- **Do not add operators if Step 2 fails.** The pre-registration's anti-tuning
  rules apply, and Project 1's value came precisely from reporting the null with
  a mechanism rather than searching for a passing operator.
- **Do not treat an ordinal 3-level C2 pass as required.** Real medium-vs-coarse
  is unresolved; demanding three levels invites fitting to noise.

---

## 5. Honest statement of the main risk

If I am wrong about anything here, it is most likely **§1.1**: that the 1.26×
specific-surface spread is too small to move an integer-bisection transition, and
Step 2 returns a C1 null that is a resolution artifact rather than a physics
result — the same failure that ended Project 1, in a new costume. **R1 is the
mitigation and it costs about an hour.** I would not run Step 2 without it.

Second most likely: the ~8× TPB over-prediction interacts with C3 in a way the
ratio agreement hides, e.g. because synthetic TPB lives on smooth sphere caps
rather than on sintered necks, so it may degrade differently under damage. Worth
watching, not worth blocking on.
