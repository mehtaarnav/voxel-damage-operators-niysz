# Impact note — porespy `snow2` parallel-mode bug, effect on the delivered real-data study

Filed 2026-08-10, during the follow-on synthetic decoupling study. This note
exists because `REPORT.md`'s conclusions must never be silently altered — this
is the explicit, standalone disclosure required before any such finding is
considered "documented, not silently absorbed."

## What was found

`porespy.networks.snow2`'s own default is `parallel_kw={}` (an empty dict,
which the function treats as truthy — i.e. "use parallel/chunked
processing"), not `None`. `cmlib.pnm.extract_ni_network` (the function behind
every SNOW network extraction in this project) did not pass `parallel_kw` at
all, so **every SNOW extraction in the delivered real-data study silently ran
in chunked mode**, including the "Image was cropped to ..." warnings visible
in the original Phase 3 logs — those warnings are a symptom of chunked mode,
not a separate issue.

Found while writing a Phase-0 unit test for the follow-on synthetic study
(`scripts/next/phase0_validate_synthetic_pipeline.py`, test T4): an idealized
two-cube dumbbell, run through the default (chunked) path, **crashed** inside
porespy's `regions_to_network` with `IndexError: too many indices for array:
array is 1-dimensional, but 2 were indexed` — even though the two watershed
regions were confirmed genuinely 6-connectivity face-adjacent. Forcing
`parallel_kw=None` (serial) fixed the crash and recovered the exact analytic
answer. Full reproduction: `probe_snow2_parallel_bug.py`.

`cmlib.pnm.extract_ni_network` now defaults to `parallel_kw=None` (serial)
for all new work. **This was not applied retroactively to the delivered
real-data outputs** (`out/phase3/*.csv`, `out/phase4/*.csv`, `REPORT.md`) —
those stand exactly as reported.

## Which original outputs were affected

Everything downstream of `cmlib.pnm.extract_ni_network`, i.e. every SNOW-based
network quantity in **Phase 3 and Phase 4** of the real-data study:
`n_pores`, `n_throats`, `neck_p10_nm`, `neck_p25/50/90_nm`, `lambda2_raw`,
`lambda2_norm`, `mincut`, `g_eff`, `mean_degree` — across all 21 pristine
network extractions and the `r_max` sensitivity sweep. **NOT affected:**
Phase 0 (percolation validation, no SNOW involved), Phase 2 (volume
fractions, plain voxel counting), Phase 4b (TPB density, voxel-edge counting,
no SNOW), Phase 4d (watershed particle sizing — uses
`skimage.segmentation.watershed` directly, not porespy's SNOW/`snow2`), and
Phase 5 (percolation on full stacks, `scipy.ndimage.label`, no SNOW).

## Quantified impact (one ROI checked; not exhaustively re-verified)

Re-ran extraction on `coarse_pre` ROI `z0y0x0` (the exact ROI reported in
`out/phase3/phase3_snow_8.0um_rmax4.csv`) under both modes:

| | pores | throats | neck p10 | neck p50 |
|---|---|---|---|---|
| **reported (chunked, undisclosed at the time)** | 62 | 145 | 199 nm | 516 nm |
| serial (`parallel_kw=None`) | 65 | 153 | 188 nm | 530 nm |
| relative shift | +4.8% | +5.5% | −5.5% | +2.7% |

This is the **only ROI checked**. Re-checking all 21 would require rerunning
the full Phase 3/4 extraction pipeline in serial mode — a non-trivial
computation (the original run took on the order of tens of minutes) that is
out of scope for this note and not needed to answer the question that
matters: does this change any conclusion in `REPORT.md`?

## Whether any ordering changed, and why `REPORT.md`'s conclusions still stand

**No ordering used in `REPORT.md`'s verdict changes.** Reasoning, checked
against the magnitude of the effect above (~3-6% per quantity on the one ROI
tested):

- **λ₂ raw's "coarse > medium > fine" ordering, and its reduction to a
  node-count artifact** (`REPORT.md` §3, "λ₂ raw is not measuring
  connectivity — it is measuring node count"): the per-anode λ₂ means differ
  by **7.65×** (fine 0.458, medium 1.228, coarse 3.504), and node counts by a
  similar multiple (509 vs 165 vs 82). A ~5% perturbation to individual ROI
  values cannot plausibly reverse a several-fold, multi-anode separation.
  Unaffected.
- **min-cut / effective-conductance's "medium > fine > coarse" ordering**
  (reported as a clean **mismatch** against the outcome ordering
  "coarse > medium > fine" — a full permutation difference, not a close
  call): a ~5% shift on individual values cannot turn a full permutation
  mismatch into a match. The FAIL verdict for these two metrics is
  unaffected.
- **neck-width p10's "coarse > medium > fine" ordering**: separation between
  anodes is large (69.5 / 156.9 / 205.3 nm mean), and `REPORT.md` **already
  discloses** that this metric's medium-vs-coarse pair is not robust — the
  `r_max` sensitivity sweep (§"Robustness to the watershed marker parameter")
  already shows medium and coarse swapping at `r_max=2` (184 vs 179 nm, a 3%
  gap). A further ~5% parallel/serial shift is smaller than that
  already-disclosed swap magnitude, so it does not introduce a new risk to
  the reported conclusion — it is consistent with, not contradictory to, the
  already-stated caveat that "the resolvable content on both sides is
  'fine is worst', and medium vs coarse is noise in predictor and outcome
  alike."
- **The core verdict** ("the metrics that reproduce the outcome ordering give
  the same ranking as mean Ni particle size; the metrics not reducible to
  coarseness fail") does not depend on the exact numeric values, only on
  which anode ranks where, which is unaffected by a ~5% per-ROI shift given
  the multi-fold separations involved.
- Phase 2, 4b, 4d, 5 conclusions (volume fraction, TPB density, particle
  size, percolation retention — the parts of `REPORT.md` doing the most
  load-bearing work) do not use SNOW at all and are entirely unaffected.

**Caveat on this reasoning:** it is a bound derived from one ROI's measured
shift plus the already-disclosed spread/sensitivity magnitudes, not an
exhaustive re-verification of all 21 ROIs in serial mode. If a future session
needs to re-certify the exact numbers (as opposed to the qualitative verdict,
which this note is confident stands), that re-run should be treated as a
deliberate, disclosed task — not folded silently into unrelated work.

## Disposition

`REPORT.md` is unchanged. This note is the disclosure. `cmlib/pnm.py`'s
docstring carries a shorter version of the same finding and points here.
