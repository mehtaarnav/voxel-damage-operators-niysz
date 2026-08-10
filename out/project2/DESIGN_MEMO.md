# Project 2 — Mechanism Fingerprinting of Ni-YSZ Redox Degradation
## Design Memo and Pilot Spec

> **UPDATE 2026-08-10 — Step 0 has since been executed.** The premise audit in
> §0 is resolved: the gate returned **outcome (a)**, coarse worst on YSZ, with
> `coarse_post` YSZ failing to percolate entirely (`P_span = 0.0000`). My
> expectation of outcome (c) (saturation) was **wrong**. Two gates below are
> amended as a result — **G1-d is revised and G1-i is added** — see
> `STEP0_REPORT.md` §5, which takes precedence over §4.2 of this memo where they
> conflict. The rest of the memo stands as written.

**Status: design memo only.** No damage code was written, no operator was
implemented, no dataset was generated. The feasibility numbers in §1 and §2 come
from a throwaway geometry-and-percolation probe using the *existing* frozen
Platform-v2 generator (rasterization + `add_ysz_pore` + `percolation_summary` +
`extract_ni_network`, all unmodified). Nothing was written to `out/` by that
probe; every configuration below is reproducible from the parameters tabulated
here.

**Awaiting approval before any pilot implementation.**

---

## 0. Premise audit — read this before §1

Project 1 closed on a lesson that applies immediately here: *before optimising
against a benchmark, establish that the benchmark is the correct one.* The
target signature for Project 2 is stated in the brief as three ordinal facts.
Two of them are verified in our own committed data. **The third — the one that
carries the whole YSZ half of the project — is not measured anywhere in this
repository.**

### 0.1 What is verified

**Fact 1 — fine loses Ni percolation worst.** Verified, robust
(`out/phase6/phase6_comparison_table.csv`):

| anode | `P_span_retained` | `P_reach_retained` | published *P* retained |
|---|---|---|---|
| fine | **0.680** | 0.857 | 0.821 |
| medium | 0.855 | 0.979 | 0.916 |
| coarse | 0.947 | 0.942 | 0.924 |

Carried caveat from Project 1: **only "fine is worst" is robust; medium vs
coarse is unresolved** and flips with the definition. Project 2 must therefore
score against a *two-level* Ni ordering (fine < {medium, coarse}), not a
three-level one. Requiring a synthetic operator to reproduce a three-level Ni
ordering would be fitting to noise.

**Fact 2 — fine retains TPB best.** Verified, and it runs *opposite* to Ni
percolation retention: 0.799 / 0.746 / 0.590 (this work), 0.743 / 0.586 / 0.608
(published). This divergence is the most valuable single feature of the target
signature, because it is the one an operator cannot reproduce by accident.

### 0.2 What is **not** verified — YSZ percolation was never measured

`phase5_percolation.py` computed percolation for the **Ni phase only**. There is
no YSZ percolation number, no YSZ node count, and no YSZ largest-component
fraction anywhere in this repository, for any of the six stacks. The brief's
Fact 3 ("coarse suffers YSZ network degradation / YSZ percolation loss") is
currently a **literature paraphrase, not a measurement we hold.**

What we *do* hold is phase volume fractions, pristine and degraded, for all
three anodes (`out/phase2/phase2_volume_fractions.csv`) — and they do not
obviously support the stated ordering:

| anode | Φ_Ni pre → post | ΔΦ_Ni | Φ_YSZ pre → post | **ΔΦ_YSZ** | Φ_pore pre → post |
|---|---|---|---|---|---|
| fine | 0.3222 → 0.2216 | **−0.1006** | 0.4213 → 0.3123 | **−0.1090** | 0.2541 → 0.4661 |
| medium | 0.2497 → 0.2332 | −0.0165 | 0.3881 → 0.3764 | **−0.0117** | 0.3622 → 0.3901 |
| coarse | 0.2293 → 0.2439 | **+0.0146** | 0.3838 → 0.3243 | **−0.0595** | 0.3869 → 0.4318 |

Three things follow, and all three change the design:

1. **By YSZ volume loss, fine is worst (−0.109), not coarse (−0.060).** Coarse
   is second, medium is nearly untouched (−0.012). Volume loss is not
   percolation loss — a backbone can fragment at near-constant volume, which is
   precisely the mechanism the brief posits for coarse — but we cannot assert
   the coarse-worst YSZ ordering from the data we have. **It has to be measured.**
2. **Ni volume fraction *increases* in coarse (+0.0146)** while its Ni
   percolation retention is the best (0.947). Any operator that destroys Ni by
   removing Ni voxels will move Φ_Ni monotonically down and can never reproduce
   this. Either the sign is a specimen effect (see 3), or Ni redistribution
   rather than Ni removal is doing the work in coarse.
3. **Pre and post are different specimens, not the same volume re-imaged**:
   pristine = Rx36/Rx37/Rx38, degraded = Rx41-1/Rx41-2/Rx41-3. Every Δ in the
   table above is therefore confounded with specimen-to-specimen variation. This
   is a property of the dataset, not a processing choice, and it bounds how
   seriously any Δ*magnitude* can be taken — which is an additional, independent
   reason the ordinal-only constraint in the brief is the right call.

### 0.3 Required consequence: a Step 0 before the pilot

**The pilot must begin by measuring YSZ percolation on the six real stacks.**
This is cheap (the pipeline is phase-agnostic — see §2 — and the Ni run took
~10 s per stack), it is the only way to know what signature we are asking the
operators to reproduce, and without it §4's Go/No-Go criterion for the YSZ half
is unfalsifiable.

Step 0 has three possible outcomes, and the honest position is to pre-commit to
all three now:

- **(a) Coarse is worst on YSZ percolation retention.** The brief's premise is
  confirmed; the pilot proceeds exactly as specified in §4.
- **(b) Fine is worst on YSZ too** (which the volume fractions weakly suggest).
  Then there is **no divergent signature to fingerprint** — fine is simply worst
  at everything except TPB — and the project's scientific question collapses to
  a much weaker one. **This is a legitimate kill outcome for Project 2 as
  framed**, and it should kill it before any operator is built, not after.
- **(c) YSZ percolation is saturated at ~1.0 in all six stacks** (plausible:
  Φ_YSZ is 0.31–0.42, well above the 6-connectivity site threshold 0.3116, and
  the degraded minimum 0.3123 sits essentially *on* it). Then YSZ percolation is
  the wrong outcome variable and the YSZ half must be re-specified onto a
  sensitive one — largest-component fraction, component count, or YSZ effective
  conductance — before operators are designed.

Outcome (c) is the one I consider most likely, on the arithmetic: fine's
degraded Φ_YSZ of 0.3123 is 0.0007 above the random-site threshold, but real YSZ
is spatially correlated rather than random, and correlated media percolate *well
below* the random threshold. So I expect YSZ still to span in all six stacks,
with the interesting signal in the *component structure*, not the spanning
boolean. §4 is written so this does not stall the pilot.

---

## 1. Synthetic analog feasibility

### 1.1 Verdict

**Yes — all three analogs are feasible inside a ~192³ domain, and I have built
and measured them.** Every analog hits its real Φ_Ni target to within 0.5 %,
lands on the real per-anode Φ_YSZ essentially exactly, percolates in Ni and in
YSZ at t = 0 as a single Ni cluster, and comfortably clears the ≥30-particle REV
floor.

**But there is one fidelity gap that is forced by voxel resolution rather than
by tuning, and it is the single most important number in this memo** — see
§1.4. It does not block the pilot; it bounds what the pilot can claim.

### 1.2 Measured configurations (probe output, structure seed 0)

Lattice from `platform_v2_lattice_geometry`; neck widths from the frozen
mixture (`frac_weak = 0.20`, weak ∈ [4,6] vox, normal ∈ [12,20] vox, validity
floor p50/p10 ≥ 2.5) scaled by a per-analog `neck_scale` bisected to hit Φ_Ni;
YSZ from `add_ysz_pore` at the **real per-anode** YSZ-share-of-remainder.

| | **FINE** | **MEDIUM** | **COARSE** |
|---|---|---|---|
| lattice `nlat_z × nlat_xy²` | 8 × 6² | 7 × 5² | 6 × 4² |
| `pitch` (vox) | 24 | 30 | 36 |
| `R` (vox) | 10.5 | 12.1 | 14.0 |
| `margin` (vox) | 8 | 8 | 8 |
| `jitter_frac` | 0.15 | 0.15 | 0.15 |
| `neck_scale` (solved) | 0.747 | 0.807 | 0.988 |
| **domain shape** | 169 × 180 × 180 | 181 × 190 × 190 | 181 × 188 × 188 |
| domain voxels | 5.48 M | 6.53 M | 6.40 M |
| **particles** | **288** | **175** | **96** |
| neck pairs | 732 | 430 | 224 |
| D_particle (nm) | 420 | 484 | 560 |
| neck p10 / p50 (nm) | 80 / 220 | 80 / 240 | 100 / 300 |
| **Φ_Ni achieved** | **0.3207** | **0.2500** | **0.2287** |
| Φ_Ni target (real) | 0.322 | 0.250 | 0.229 |
| Φ_Ni deviation | **−0.41 %** | **−0.02 %** | **−0.12 %** |
| `ysz_frac_of_rest` | 0.6237 | 0.5173 | 0.4981 |
| **Φ_YSZ achieved** | **0.4237** | **0.3880** | **0.3841** |
| Φ_YSZ target (real) | 0.421 | 0.388 | 0.384 |
| **Ni percolates (t=0)** | **yes**, P_span 1.0000, 1 cluster | **yes**, 1.0000, 1 cluster | **yes**, 1.0000, 1 cluster |
| **YSZ percolates (t=0)** | **yes**, P_span 0.9997 | **yes**, 0.9989 | **yes**, 0.9987 |
| YSZ components (t=0) | 167 | 246 | 190 |
| build time / seed | ~16 s | ~17 s | ~7 s |

All three max dimensions (180 / 190 / 188) sit inside 192. Φ_Ni is hit by
bisecting `neck_scale` at fixed lattice — a 9-iteration 1-D solve reusing the
same pattern as the existing `match_radius_for_mass_conservation`.

**Coarse REV alternative, also measured.** If 96 particles is judged too thin
for the coarse analog, `6 × 5²` at `pitch = 34`, `R = 14.0` gives **150
particles**, Φ_Ni = 0.2295 (+0.20 %), Ni and YSZ both percolating — at the cost
of a 171 × **214** × 214 domain, which **exceeds 192** on two axes (9.8 Mvoxel).
Recommendation in §1.5.

### 1.3 Resolution and percolation limits — the actual binding constraints

Taking the brief's two stated limits in turn, plus two more the probe exposed:

**(i) REV / particle-count floor — satisfied with margin.** 288 / 175 / 96
particles against the stated ≥30 floor. Coarse is the tightest at 96 and is
still 3.2× the floor. Effective interior counts (excluding the half-spheres on
the two z-boundary layers) are 252 / 150 / 80. **Not a constraint at these
parameters.**

**(ii) Initial percolation — satisfied, and by construction for Ni.** All Ni
analogs give a *single* connected cluster with P_span = 1.0000. This is
guaranteed by the generator's design (neck bars span centre-to-centre, so
connectivity survives any radius shrink) and by the boundary layers sitting
exactly on the z-faces. YSZ spans in all three at P_span ≥ 0.9987 without any
tuning. **Not a constraint.**

**(iii) Neck resolution floor — the real constraint, and it binds.** See §1.4.

**(iv) Domain-size arithmetic.** `nz = (nlat_z−1)·pitch + 1` and
`n_xy = nlat_xy·pitch + 2R + 2·margin`. The z-formula carries no margin (the
spanning guarantee requires boundary layers on the faces), so the two axes
cannot be squared by a single particle count — the existing generator already
decouples `nlat_z` from `nlat_xy` for exactly this reason, and the configurations
above use that.

### 1.4 The forced fidelity gap: neck proportion vs voxel resolution

The analogs match real Φ_Ni, real Φ_YSZ, and the real *ordering* of particle
size. They do **not** match the real proportion between neck size and particle
size, and this is not a tuning failure:

| | model neck_p50/D | real neck_p50/D | model neck_p10/D | real neck_p10/D |
|---|---|---|---|---|
| fine | 0.524 | **0.259** | 0.190 | **0.061** |
| medium | 0.496 | **0.358** | 0.165 | **0.109** |
| coarse | 0.536 | **0.368** | 0.179 | **0.120** |

**Our necks are proportionally ~1.5–2× fatter than real, and — worse for this
project — the model trend is flat/non-monotone (0.524 / 0.496 / 0.536) where the
real trend increases monotonically with coarseness (0.259 → 0.358 → 0.368).**

**Why it cannot be tuned away at 192³.** To reproduce the real fine anode's
`neck_p10/D = 0.0605` while keeping the narrowest necks resolved at the
pre-registered ≥4-voxel floor requires

  D ≥ 4 / 0.0605 = **66 voxels = 1320 nm**

i.e. essentially the real particle diameter (1148 nm). Combined with a ≥6
particles-per-axis REV requirement, the fine analog then needs

| analog | real neck_p10/D | min D for p10 ≥ 4 vox | min axis at 6 particles |
|---|---|---|---|
| fine | 0.0605 | 66.1 vox (1322 nm) | **~397 vox** |
| medium | 0.1086 | 36.8 vox (737 nm) | ~221 vox |
| coarse | 0.1197 | 33.4 vox (669 nm) | ~200 vox |

**Fine is the binding case and it needs a ~400³ domain — about 8× the voxel
count of the current design.** Medium and coarse are nearly feasible at 192–224³
already. Our current analogs sit at D ≈ 21–28 voxels, i.e. **3× compressed**
relative to real (model/real = 0.366 / 0.335 / 0.327), and the compression is
what forces the fat necks: at D = 21 vox, a 4-voxel neck *is* 19 % of the
particle.

**Is ~400³ actually out of reach?** No — 397³ = 62.6 Mvoxel, which is *within*
the ~120–150 Mvoxel ceiling the real-data study operated under. It costs roughly
8–10× the build and damage time per seed, and SNOW extraction is the expensive
step (§2.3). So this is a real, affordable structural option, not a fantasy.

### 1.5 Recommendation

**Run the pilot at the compressed 192³ scale in the table of §1.2, and treat the
~400³ faithful-proportion build as a named, costed follow-on rather than a
pilot prerequisite.** Reasoning:

- Project 2's claims are **ordinal by constraint**. Uniform 3× size compression
  is close to *common-mode* across the three analogs (model/real spans only
  0.327–0.366, a 12 % spread), so it does not differentially advantage any
  analog in an ordinal comparison.
- The neck-proportion mismatch, by contrast, is **not** common-mode — it is
  flat where the real trend rises. That is a genuine differential bias and it
  must be stated as a limitation in the artifact, with its direction: our fine
  analog has necks that are *proportionally too robust* relative to real fine,
  which biases **against** reproducing "fine loses Ni percolation worst". The
  bias works against the hypothesis, so a positive pilot result is conservative;
  a null is partially confounded with this bias. **That asymmetry must be
  written into the Go/No-Go interpretation (§4.5), not discovered afterwards.**
- Spending the 8× up front, before we know from Step 0 whether a divergent
  signature even exists in the real data, would repeat Project 1's pattern of
  building a platform ahead of the question.

**The structural alternative, stated explicitly so the choice is yours:** build
the pilot directly at ~400³ for fine and ~224³ for medium/coarse, accepting
non-uniform domains and ~8× compute on the fine arm, in exchange for
neck-proportion fidelity and a monotone neck/D trend. My recommendation is
against it *for the pilot*, and for it *if the pilot passes Go* — at which point
the fidelity gap becomes the leading threat to the main result and is worth
paying to close.

---

## 2. YSZ network pipeline

### 2.1 Verdict: the existing pipeline already works on YSZ, unmodified

Both of the relevant modules are **phase-agnostic in their bodies and only
Ni-specific in their names**:

- **`cmlib/percolation.py`** — `percolation_summary`, `percolation_report`,
  `percolating_mask`, `label_phase` all take a bare boolean mask. **Nothing to
  add.** YSZ percolation, YSZ component count, and YSZ largest-component
  fraction (`P_largest`) are available today by passing `ysz_mask`. Confirmed
  by direct execution: §1.2's YSZ rows were produced this way.
- **`cmlib/pnm.py::extract_ni_network`** — despite the name, its body takes a
  generic mask and calls `ps.networks.snow2` on it. **Confirmed to run on YSZ**
  (Platform-v2 medium analog, 181 × 190 × 190):

| phase | SNOW nodes | throats | throat p10 | p50 | max | time |
|---|---|---|---|---|---|---|
| Ni | 175 | 392 | 80 nm | 240 nm | 320 nm | 26 s |
| **YSZ** | **1,889** | **4,890** | **40 nm** | **113 nm** | 312 nm | **48 s** |

The Ni node count recovering the particle count exactly (175 nodes / 175
particles) is a good sanity signal for the extraction.

### 2.2 Four minimal additions actually needed

None is large; all are mechanical. Listed with why the existing code cannot
serve as-is.

**A2.1 — `extract_network` must be able to run on a non-spanning phase.**
`extract_ni_network` calls `percolating_mask(...)` and **returns `(None, diag,
{})` if the phase does not span** (`reason = "no spanning Ni cluster"`). For
Project 2 that is fatal: the entire YSZ question is what the network looks like
*after* it stops percolating, and O3's whole purpose is to drive it there. The
minimal fix is a `restrict_to` argument — `"spanning"` (current behaviour,
default, so no existing result changes), `"largest"`, or `"all"` — with Project
2 using `"largest"` for post-damage analysis. Rename to `extract_network` with
`extract_ni_network` retained as a thin alias so nothing in Project 1 breaks.

**A2.2 — `apply_d4` hard-codes YSZ as immutable.** It takes `ysz_mask` only to
protect it (dilation is restricted to original pore, and "YSZ is never modified
by damage" is asserted at the call site). O3 damages YSZ, so the operator
signature must generalise to *(target phase, protected phase, host pore)*. This
is a **new operator module (`cmlib/damage2.py`), not an edit to `damage.py`** —
`apply_d4` is frozen by Project 1's pre-registration and must stay
bit-reproducible. Same reasoning for `rebuild_ternary`, which assumes an
unchanged YSZ mask.

**A2.3 — YSZ morphology must be allowed to vary with analog coarseness.**
This is the substantive one, and it is a **prerequisite for O3 having any lever
at all.** `add_ysz_pore` builds YSZ by thresholding a Gaussian-smoothed random
field at a **fixed `smooth_sigma_vox = 3.0`**. That means: in the §1.2 analogs,
YSZ morphology is *identical in length scale* across fine, medium and coarse —
only its volume fraction differs (0.424 / 0.388 / 0.384). **A YSZ operator
applied to three YSZ phases that differ only in Φ cannot produce a
coarseness-ordered YSZ signature except through Φ.** The fix is to scale σ with
the analog's particle size. Measured lever, medium analog, Φ_YSZ held at 0.3880:

| `smooth_sigma_vox` | YSZ components | YSZ EDT p50 | EDT p90 | percolates | P_span |
|---|---|---|---|---|---|
| 1.5 | 615 | 28.3 nm | 44.7 nm | yes | 0.9992 |
| 3.0 *(current)* | 246 | 40.0 nm | 80.0 nm | yes | 0.9989 |
| 5.0 | 98 | 56.6 nm | 121.7 nm | yes | 0.9987 |
| 7.0 | 47 | 66.3 nm | 161.2 nm | yes | 0.9977 |
| 10.0 | 38 | 82.5 nm | 220.0 nm | yes | 0.9942 |

σ is a clean, monotone control on YSZ length scale at **fixed Φ_YSZ**, and YSZ
percolates throughout — so it is safe to use. Proposed rule: **σ ∝ analog
particle diameter**, anchored so medium keeps its current σ = 3.0, giving
σ ≈ 2.6 / 3.0 / 3.5 for fine / medium / coarse. This is a *placement*
parameter, not a damage parameter, and must be frozen and pre-registered before
any operator is run.

**A2.4 — a YSZ metric set.** Trivial but must be named now so it cannot be
chosen after seeing results: for YSZ we record `percolates`, `P_span`,
`P_reach`, `P_largest`, `n_clusters`, `n_phase_voxels`, plus SNOW `n_nodes`,
`n_throats`, and throat-size percentiles. All are already returned by the
existing functions.

### 2.3 One honest limitation on YSZ SNOW throats

The YSZ phase is a **smoothed random field, not a designed sintered backbone**.
Its 1,889 watershed regions and 4,890 throats are therefore a property of the
placement field, not of a designed structure, and its throat p10 of **40 nm =
2 voxels sits below the ≥4-voxel resolution floor** we hold ourselves to for Ni
necks. Consequences, which shape §3:

- YSZ throat *identity* is not a controlled design variable the way Ni neck
  width is. There is no YSZ analogue of the max-clip intervention.
- **O3's primary variant must therefore be YSZ surface erosion** (well-defined,
  resolution-independent), with **targeted thin-throat removal as a secondary
  variant** whose throat population is acknowledged to be generator-defined.
  Building O3 the other way round would rest the YSZ conclusion on a
  sub-resolution measurement.
- SNOW on YSZ costs ~48 s per volume against ~26 s for Ni, and it is the
  dominant per-evaluation cost. §4 budgets accordingly, and does **not** call
  SNOW inside the damage-intensity bisection loop.

---

## 3. Phenomenological operator definitions

All three are **voxel/graph operators only** — no strain solver, no
multi-physics, per constraint 2. Each is specified algorithmically here and
**not implemented**. Common conventions: 6-connectivity throughout; damage
intensity is an integer `n_rounds`; the *only* free intensity variable per
operator is `n_rounds`, with all other parameters frozen before any comparison;
removed solid becomes pore; each operator is applied to a *pristine* structure,
never composed with another unless the composition is itself pre-registered
(§3.4).

### O1 — Ni surface erosion (adapted D4)

**Mechanism claim:** redox cycling roughens and thins Ni globally in proportion
to exposed surface area; Ni agglomerates and sheds disconnected islands.

**Algorithm** (identical to frozen `apply_d4` except where marked):

```
input: ni_mask, ysz_mask, n_rounds, p_erode, expand_vox, seed
1. pore0 <- NOT ni_mask AND NOT ysz_mask
2. ni <- ni_mask OR (dilate(ni_mask, STRUCT6, expand_vox) AND pore0)
      # oxidative expansion, may claim pore, never claims YSZ
3. repeat n_rounds times:
       boundary <- ni AND NOT erode(ni, STRUCT6)
       ni <- ni AND NOT (boundary AND uniform(shape) < p_erode)
4. keep only the largest 6-connected component of ni; the rest become pore
5. YSZ unchanged (assert array equality)
output: damaged ni_mask, diagnostics
```

**Parameters:** `p_erode = 0.35`, `expand_vox = 1` — inherited frozen from
Project 1, with their provenance caveat carried verbatim (E0 spike origin, never
re-derived). **Do not re-derive them for Project 2 either**; inheriting them
unchanged is what makes O1 a genuine control rather than a fourth tuned model.

**Predicted fingerprint** (stated *before* the pilot, so it is falsifiable):
O1 is surface-area-mediated — Project 1 established this by direct measurement
of the failure step. Specific surface area scales as ~1/particle size, so **fine
should lose Ni fastest**, reproducing Fact 1. O1 does **not** touch YSZ, so it
must produce **no YSZ signature at all**. It therefore cannot alone reproduce the
divergence, by construction — it is the Ni-half candidate and the null control
for the YSZ half.

**Known liability:** O1 removes Ni voxels monotonically and so cannot reproduce
the coarse Φ_Ni *increase* noted in §0.2.

### O2 — Ni contact/neck failure

**Mechanism claim:** redox strain concentrates at Ni-Ni contacts; narrow necks
sever, disconnecting Ni regardless of how much Ni volume survives.

**Algorithm:**

```
input: ni_mask, n_rounds, p_sever, neck_percentile_cut, seed
1. G, extras <- extract_network(ni_mask, restrict_to="largest")   # SNOW
2. w <- throat inscribed diameters; T <- percentile(w, neck_percentile_cut)
3. candidates <- throats with w < T
4. repeat n_rounds times:
       for each surviving candidate throat t, independently with prob p_sever:
           sever t: set to pore every voxel of the watershed interface
                    between its two regions, dilated by 1 voxel to guarantee
                    6-connectivity separation
5. keep only the largest 6-connected component; the rest become pore
6. YSZ unchanged
output: damaged ni_mask, diagnostics (n_severed, voxels_removed)
```

**Parameters to freeze before use:** `neck_percentile_cut` (proposal: 25th),
`p_sever` (proposal: 0.25). Intensity variable: `n_rounds`.

**Predicted fingerprint:** fine has the narrowest necks in absolute terms
(real: 69.5 vs 205.3 nm) and the most necks per volume, so **fine should lose Ni
percolation worst and fastest**. Unlike O1, O2 achieves this at *minimal volume
loss* — severing interfaces removes few voxels. **That is O2's discriminating
signature against O1: same Ni ordering, very different ΔΦ_Ni.** This is the one
place where a magnitude is diagnostic rather than a fitting target, and it is
diagnostic *between operators*, not against real numbers.

**Pre-registered interpretation limit, carried from Project 1's §3 D1 rule.**
O2 selects on neck width and then observes that neck-poor structures fail —
**this is the near-tautological D1 pattern.** Project 1 froze the rule that an
effect appearing under a width-selecting operator but not under a
non-selecting one scores NEGATIVE. Project 2 must not quietly evade that rule
by relabelling D1 as O2. **O2 is admissible here only because the object of the
test has changed**: we are not asking "does neck width matter" (assumed by
construction) but "does a neck-severing mechanism reproduce the *cross-analog
ordinal signature* including the TPB divergence" — a fingerprint that O2 can
fail. This distinction must be written into the pre-registration explicitly; if
it cannot be defended there, O2 should be dropped rather than smuggled in.

### O3 — YSZ backbone failure

**Mechanism claim:** redox volume change mechanically damages the YSZ skeleton;
the YSZ network fragments and loses percolation.

**Primary variant, O3a — YSZ surface erosion.** Structurally identical to O1
with the phase roles swapped, and **no oxidative expansion step** (YSZ does not
redox-expand; its damage is mechanical):

```
input: ysz_mask, ni_mask, n_rounds, p_erode_ysz, seed
1. repeat n_rounds times:
       boundary <- ysz AND NOT erode(ysz, STRUCT6)
       ysz <- ysz AND NOT (boundary AND uniform(shape) < p_erode_ysz)
2. removed YSZ becomes pore
3. Ni unchanged (assert array equality)
4. do NOT apply largest-component pruning to YSZ
output: damaged ysz_mask, diagnostics
```

**The island-removal step is deliberately omitted for YSZ**, and this is a
substantive asymmetry, not an oversight: disconnected Ni is *electrically* dead
and correctly deleted; disconnected YSZ is still physically present and still
ionically relevant to its local neighbourhood. Deleting it would manufacture
volume loss that the mechanism does not imply. Consequence: for YSZ, the outcome
of interest is **fragmentation** (`n_clusters`, `P_largest`, `P_span`) at
near-constant Φ_YSZ — exactly the signature §0.2 says we cannot see in volume
fractions alone.

**Secondary variant, O3b — targeted thin-YSZ-throat removal.** Same algorithm
as O2 with the YSZ SNOW network as input. **Explicitly demoted to secondary**
per §2.3: the YSZ throat population is generator-defined and its p10 sits at
2 voxels, below our resolution floor. Run it as a robustness check on O3a's
direction, never as the primary evidence.

**Predicted fingerprint:** at fixed `p_erode_ysz`, surface-mediated YSZ erosion
again favours coarse-grained YSZ. With A2.3's σ-scaling in place, **fine's
finer-scale YSZ should fragment first** — which, note, is the *opposite* of the
brief's Fact 3 and consistent with the volume-fraction table in §0.2. **If O3a
robustly produces fine-worst YSZ, then either the premise or the operator is
wrong, and Step 0 is what tells us which.** I am flagging this now because it is
the most likely way the pilot fails, and I would rather have predicted it than
explained it afterwards.

### 3.4 Composition

The brief asks which operator "or combination" reproduces the signature. Two
rules, to be frozen before the pilot:

1. **Single operators are evaluated first and alone.** A combination is
   considered only if no single operator passes, and only in the specific
   pairing {O1 or O2} × O3a — Ni mechanism plus YSZ mechanism — because that is
   the only composition with a mechanistic rationale.
2. **A composition introduces a free mixing parameter** (relative intensity of
   the two operators). That is a fitting knob, and fitting is forbidden by
   constraint 1. It is admissible **only** at a fixed, pre-registered ratio
   (proposal: equal `n_rounds` for both), never tuned to improve the match.

---

## 4. Pilot kill-test plan

**Duration: 2 weeks. Purpose: kill or confirm, not to produce the main result.**

### 4.1 Step 0 (day 1–2) — measure the target signature on the real data

Run YSZ percolation on all six real stacks using the existing phase-agnostic
`percolation_summary` (no code change needed beyond a phase selector in the
Phase 5 driver). Record for each stack: `percolates`, `P_span`, `P_reach`,
`P_largest`, `n_clusters`, `n_phase_voxels`. Compute retained YSZ percolation
per anode.

**Also re-run the truncated Phase 5 Ni table** (see the defect note in
`out/writeup/REPRODUCIBILITY_MANIFEST.md` §2.2) so both phases land in one
complete, committed CSV.

**Step 0 gate:**
- **(a) coarse worst on YSZ retention** → proceed to §4.2 as written.
- **(b) fine worst on YSZ retention** → **STOP.** No divergent signature exists;
  report the finding and re-scope Project 2 before building anything.
- **(c) YSZ percolation saturated (all six ≥ 0.99)** → proceed, but **substitute
  the YSZ outcome variable** with `P_largest` and `n_clusters`, re-deriving §4.4's
  YSZ criterion on the substituted variable **before** any synthetic run.

Nothing else starts until Step 0 is committed.

### 4.2 Step 1 (day 3–5) — analog qualification

Build **5 structure seeds per analog class** (fine / medium / coarse) at the
§1.2 parameters, with A2.3's σ-scaling. This is 15 structures; at ~7–17 s per
build that is minutes, so 5 seeds costs nothing and avoids Project 1's
under-seeding argument entirely.

Qualification gates, all per-seed, all recorded before any damage:

| gate | criterion |
|---|---|
| G1-a | Φ_Ni within ±2 % of the real anode's value |
| G1-b | Φ_YSZ within ±2 % of the real anode's value |
| G1-c | Ni percolates, single cluster, P_span = 1.000 |
| G1-d | YSZ percolates, P_span ≥ 0.99 |
| G1-e | particle count ≥ 30 |
| G1-f | Ni neck p10 ≥ 4 voxels (resolution floor) |
| G1-g | particle-size ordering fine < medium < coarse holds on measured SNOW node size |
| G1-h | **YSZ length-scale ordering fine < medium < coarse** on YSZ EDT p50 (verifies A2.3 actually did something) |

**G1-h is the one that can fail.** If YSZ length scale does not order with
coarseness after σ-scaling, O3 has no lever and the YSZ half must be redesigned
before proceeding.

### 4.3 Step 2 (day 6–10) — operator application

For each of the 3 operators × 3 analog classes × 5 structure seeds × **3 damage
seeds** (independent of structure seeds), find the damage intensity by the same
**integer bisection on `n_rounds`** used in Project 1 (bracket [1, 20],
expand-only, narrow to width ≤ 1). 135 bisections.

Two targets are bisected, separately:
- **Ni percolation-loss intensity** `n*_Ni` (for O1, O2)
- **YSZ percolation-loss intensity** `n*_YSZ`, or if Step 0 returned (c), the
  intensity at which `P_largest` first drops below 0.90 (for O3a)

**Cost control, learned from §2.3:** SNOW is *not* called inside the bisection
loop. Percolation is a cheap `ndi.label`. SNOW is called only twice per
structure — once at t = 0 and once at the transition state — for the network
metrics and TPB.

At the transition state, record for both phases: Φ, `P_span`, `P_reach`,
`P_largest`, `n_clusters`, SNOW nodes/throats, and TPB density.

**Frozen before Step 2 runs:** `p_erode = 0.35`, `expand_vox = 1` (O1,
inherited); `neck_percentile_cut = 25`, `p_sever = 0.25` (O2);
`p_erode_ysz = 0.35` (O3a); σ-scaling rule (A2.3). No parameter may be adjusted
after any result is seen. Damage-seed averaging is mandatory, per Project 1
§0g/1 — **no single-damage-seed comparison is permitted**, and any conclusion
resting on a group difference smaller than 1.0 round must be reported as
unresolved rather than interpreted.

### 4.4 Step 3 (day 11–14) — scoring and Go/No-Go

Each operator is scored on **ordinal reproduction only** — never on numerical
proximity to any real Δ (constraint 1). Three criteria:

**C1 — Ni ordering.** Fine loses Ni percolation at a **strictly lower** damage
intensity than both medium and coarse: `n*_Ni(fine) < n*_Ni(medium)` and
`n*_Ni(fine) < n*_Ni(coarse)`. Deliberately a **two-level** test (§0.1) —
medium vs coarse is unresolved in the real data and is **not** scored.

**C2 — YSZ ordering.** Coarse loses YSZ percolation (or `P_largest`, per Step 0
outcome (c)) at a **strictly lower** intensity than fine. *Direction to be set
by Step 0's measured result, not by the brief's premise.*

**C3 — divergence.** The operator (or pair) produces C1 and C2 **simultaneously**,
and reproduces the TPB divergence: retained TPB ordered **fine > medium >
coarse**, opposite to retained Ni percolation.

**Resolution requirement on every criterion:** a "strictly lower" claim counts
only if the group-mean difference is **≥ 1.0 damage round** — the same
interpretability threshold Project 1 froze, and for the same reason (integer
bisection cannot resolve less). Differences below 1.0 round are recorded as
**unresolved**, never as weak support.

### 4.5 Go / No-Go criteria

| outcome | criteria met | decision |
|---|---|---|
| **GO** | some single operator or the pre-registered pair meets **C1 + C2 + C3** | Proceed to full Project 2. First act: close the §1.4 fidelity gap at ~400³ before any headline claim. |
| **CONDITIONAL GO** | **C1 + C2** met but **C3** (TPB divergence) fails | Proceed **only** on a re-scoped question: "which operator reproduces the percolation orderings", with the TPB divergence recorded as unreproduced and named as the leading open problem. Not a mechanism fingerprint. |
| **PARTIAL — Ni only** | **C1** met, **C2** fails under all operators | **No-Go on the YSZ half.** The YSZ placement model (smoothed random field) is the prime suspect, not the operators — see §2.3. Report, and re-scope to a Ni-only mechanism question or invest in a genuine sintered-YSZ generator first. |
| **PARTIAL — YSZ only** | **C2** met, **C1** fails | Investigate the §1.4 bias first: our fine analog's proportionally fat necks bias *against* C1. Re-test C1 at ~400³ **before** concluding No-Go. |
| **NO-GO** | no operator meets C1 or C2, or all differences are < 1.0 round | Stop. Report as a second mechanism-specific null. Do not add operators, tune parameters, or extend seeds to rescue it. |
| **KILL AT STEP 0** | Step 0 returns outcome (b) | Stop before any synthetic work. The premise is not in the data. |
| **KILL AT STEP 1** | G1-h fails | Stop the YSZ half; O3 has no lever until the YSZ generator is redesigned. |

**Two anti-patterns pre-committed against**, both drawn from Project 1's
experience:

1. **No post-hoc operator addition.** If all three fail, that is the result. A
   fourth operator invented after seeing the failures is fitting by search, and
   the ordinal-only constraint does not protect against it.
2. **No parameter rescue.** `p_erode`, `p_sever`, `neck_percentile_cut` and σ
   are frozen at Step 2. If an operator misses, it misses at its frozen
   parameterization, and that is what gets reported — exactly as `p_erode = 0.35`
   was reported in Project 1 rather than swept to find an effect.

**Interpretation asymmetry that must travel with any result (from §1.5):** the
compressed analogs give our fine structure proportionally fatter necks than real
fine, biasing **against** C1. A **pass** on C1 is therefore conservative and
strong. A **failure** on C1 is partially confounded with the fidelity gap and is
**not** a clean null — it routes to the "PARTIAL — YSZ only" investigation path
above rather than to NO-GO. This asymmetry is stated now, before results exist.

### 4.6 Effort summary

| step | days | compute |
|---|---|---|
| 0 — real-data YSZ percolation | 1–2 | ~10 s × 6 stacks + I/O on 0.48–1.11 Gvoxel stacks |
| 1 — analog qualification, 15 structures | 3–5 | ~15 × 20 s build + ~15 × 75 s SNOW (both phases) ≈ 25 min |
| 2 — 135 bisections, ~6 evals each | 6–10 | dominated by ~810 damage evaluations + 90 SNOW calls; est. 8–14 h wall, parallelisable by seed |
| 3 — scoring, memo | 11–14 | negligible |

---

## 5. What this memo commits to, and what it does not

**Committed:** the analog parameters in §1.2 (measured, not proposed); the four
`cmlib` additions in §2.2; the three operator algorithms in §3; the gates and
criteria in §4. All of these should be frozen in a Project 2 pre-registration
**before** any operator code is written, following Project 1's practice.

**Not committed, and deliberately deferred:** the ~400³ faithful-proportion
build (§1.5); any real-data calibration of operator parameters (forbidden by
constraint 1 anyway); any composition beyond the single pre-registered pairing
(§3.4); a genuine sintered-YSZ generator (only if §4.5's "PARTIAL — Ni only"
branch is reached).

**Open risks, ranked by how likely they are to end the project:**

1. **Step 0 returns (b) or (c)** — the premise is not in our data. Highest
   probability, cheapest to check, and checked first by design.
2. **G1-h fails** — σ-scaling does not make YSZ morphology order with
   coarseness, so O3 has no lever.
3. **The §1.4 neck-proportion bias** confounds a C1 failure.
4. **O2's D1-tautology exposure** (§3) — if the pre-registration cannot defend
   the distinction, O2 is dropped and the Ni half rests on O1 alone.

---

*Stopping here per instruction. No operator implemented, no pilot run, no
`cmlib` modified. Awaiting approval.*
