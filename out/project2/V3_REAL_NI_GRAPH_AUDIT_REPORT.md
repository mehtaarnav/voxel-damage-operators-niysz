# Ni Platform v3 (Option B) — real Ni graph audit: **lattice diagnosis CONFIRMED; ordering clause not evaluable (coarse sub-REV)**

Run 2026-08-11 under `PREREG_NI_PLATFORM_V3.md` (frozen `50b0bde` **before** the
script was written). Audit metrics unchanged from `267efd3`. Code:
`scripts/project2/v3_real_ni_graph_audit.py`. Data:
`v3_real_ni_graph_audit.csv`, `v3_real_fragility_curves.csv`.

Read-only. No damage operator. `cmlib/damage.py`, `cmlib/synth.py`,
`cmlib/project2.py` untouched. **No O4, no O5v2.**

---

## Verdict

**Criteria 1–4 PASS. Criterion 5 splits: the `f_crit ≤ 0.05` half passes
decisively in all three anodes; the `fine ≤ coarse` ordering half fails, but the
coarse ROIs are sub-REV (48–72 nodes) and that clause is therefore reported as
NOT EVALUABLE rather than failed.**

**The audit's central diagnosis is confirmed.** Real Ni networks are disordered,
heterogeneous, and **2–4× more fragile** than the lattice platform, with minimum
cuts that are a *handful of edges* rather than a full cross-section.

---

## 1. The lattice was the problem — direct comparison

| | **synthetic lattice** | **real (this audit)** |
|---|---|---|
| min-cut, fine | 36 edges = **6²**, exactly | **16** edges of 1016 |
| min-cut fraction, fine | **0.0492** | **0.0210** |
| min-cut fraction, medium | 0.0581 | 0.0233 |
| min-cut fraction, coarse | 0.0714 | 0.0123 |
| **seed/ROI variance of the fraction** | **exactly 0.0000** | 0.0154–0.0217 / 0.0108–0.0382 / 0.0120–0.0174 |
| min cut = full cross-section? | **yes, always** | **no** — 16, 6, 2 edges |
| voxel-validated break possible for fine? | **no (unbreakable)** | n/a (graph-level, see §4) |
| coordination sd | fixed by lattice | **2.04–2.41** |
| CV of connected-pair distance | ~0 (periodic) | **0.38–0.45** |

**Every acceptance criterion for disorder is met:**

1. **Pair-distance distribution is non-lattice** — CV 0.379–0.447 across all nine
   ROIs. A lattice has sharp peaks and CV ≈ 0.
2. **Coordination variance is non-trivial** — sd 1.82–2.85, far above the ≥ 1.0
   bar, against a lattice whose interior degree is fixed at 6.
3. **Min-cut fraction varies ROI to ROI** — up to 3.5× within an anode (medium
   0.0108–0.0382), against the lattice's *exact* zero variance.
4. **Min cut is nowhere near a cross-section** — 16 / 6 / 2 edges out of
   1016 / 262 / 115.

**This is the measurement the previous self-prompt lacked.** A hard-sphere
generator (Option A) now has concrete targets: mean degree ≈ 3.4–3.9 with
sd ≈ 2.0–2.4, connected-pair-distance CV ≈ 0.41, and a min-cut fraction of
~1–2 %.

## 2. Full results

| anode | ROI | nodes | edges | min-cut | frac | mean deg | sd deg | CV dist | overlap w/ lower-quartile necks |
|---|---|---|---|---|---|---|---|---|---|
| fine | 0 | 430 | 761 | 16 | 0.0210 | 3.54 | 2.03 | 0.392 | 0.375 |
| fine | 1 | 541 | 1041 | 16 | 0.0154 | 3.85 | 2.41 | 0.445 | 0.375 |
| fine | 2 | 512 | 1016 | 22 | 0.0217 | 3.97 | 2.59 | 0.408 | 0.227 |
| medium | 0 | 137 | 262 | 10 | 0.0382 | 3.82 | 2.65 | 0.395 | 0.200 |
| medium | 1 | 164 | 257 | 6 | 0.0233 | 3.13 | 1.89 | 0.416 | 0.167 |
| medium | 2 | 162 | 278 | 3 | 0.0108 | 3.43 | 2.31 | 0.439 | 0.333 |
| coarse | 0 | **66** | 162 | 2 | 0.0123 | 4.91 | 2.85 | 0.447 | 0.000 |
| coarse | 1 | **72** | 115 | 2 | 0.0174 | 3.19 | 1.82 | 0.379 | 0.000 |
| coarse | 2 | **48** | 83 | **1** | 0.0120 | 3.46 | 2.04 | 0.416 | 1.000 |

## 3. Why the ordering clause is not evaluable

Criterion 5 required `f_crit(fine) ≤ f_crit(coarse)`. Measured medians: fine
0.0210, medium 0.0233, **coarse 0.0123** — coarse appears *most* fragile.

**I am not interpreting that**, for a reason recorded in the pre-registration
before the run: **the coarse anode is sub-REV at 8 µm ROIs.** The numbers make
the problem concrete — coarse ROIs contain **48–72 nodes and 83–162 edges**, and
their minimum cuts are **1–2 edges**. A min-cut of one edge in a 48-node graph is
a finite-size artifact: the network barely spans the ROI, so a single throat
carries it. The `overlap` column shows the same instability — 0.000, 0.000,
1.000 across three ROIs, which is what a statistic computed on a 1–2 element set
looks like.

Fine (430–541 nodes) and medium (137–164 nodes) are usable; coarse is not.
**Per the constraint "if the audit reveals a boundary artifact, stop and report
before interpreting class differences", the fine-vs-coarse comparison is
withheld.**

The fragility curves show the same artifact: coarse ROIs lose spanning under
random removal at 10 % (span rate 0.778) while fine survives **50 % random
removal in every ROI and every repeat**. That is not evidence that coarse
microstructure is fragile; it is evidence that a 48-node graph is fragile.

## 4. Honest limitation carried from the pre-registration

**This pass is graph-level only.** Removing a SNOW throat from a real
segmentation has no unambiguous voxel realisation — there is no generator neck
to zero out — so the voxel validation that governed the synthetic audit could
not be applied. The synthetic audit's most important single result (fine's graph
cut failing to disconnect the voxels) has no counterpart here. Graph-level cuts
are therefore an **upper bound on fragility**: the true voxel min-cut can only be
larger.

## 5. What this settles, and what it does not

**Settled.** The synthetic lattice's planar, invariant, full-cross-section
minimum cut is an artifact of the lattice, not a property of Ni networks. Real Ni
networks fail at **1–2 % of throats**, with cuts that vary by ROI and are small,
localised edge sets. Four Ni operators failed on a platform whose topology could
not fail the way real Ni does. **The diagnosis that motivated v3 is correct.**

**Not settled.** Whether pristine min-cut fragility predicts the *measured*
degradation ordering (fine worst, retention 0.680). Fine's min-cut fraction is
not the smallest among usable anodes (fine 0.0210 vs medium 0.0233 — a 10 %
difference on n = 3 ROIs, which is noise). **There is currently no evidence that
pristine topological fragility explains C1**, and that is the open question worth
attacking next.

## 6. Go/No-Go

**Go for Option B as the platform** — real ROI graphs meet every disorder
criterion, with the coarse REV limitation binding.

**Two things must happen before any damage operator runs on it:**

1. **Fix the coarse REV problem.** Larger ROIs for coarse (Project 1 found
   coarse needs ≳ 10 µm and is memory-bound at ~120–150 Mvoxel), or accept
   fine-and-medium-only comparisons and state it.
2. **Decide the voxel-validation question** (§4). Damage on real ROIs *can* be
   voxel-validated — an operator that removes voxels has an unambiguous voxel
   result — so this limitation applies to the *audit*, not to future damage work.

**No-Go, still, on O4, O5v2, and any Ni operator on the lattice platform.**

## 7. Correction to my own earlier recommendation

My previous self-prompt proposed bond dilution on the regular lattice. That was
wrong for the reason given in the objection, and this audit now quantifies why:
the property the lattice lacks is not bond sparsity but **positional disorder** —
CV(pair distance) 0.41 and coordination sd 2.2 in real networks, both structurally
unreachable by deleting bonds from a periodic grid.
