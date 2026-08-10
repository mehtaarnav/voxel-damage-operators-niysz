# Ni Platform v3 — design choice and pre-registration amendment

**Frozen 2026-08-11, before any v3 code is written.** Amends
`PREREGISTRATION_V2_1.md`. All anti-tuning rules carry over. `cmlib/damage.py`
and `cmlib/synth.py` remain untouched.

---

## 1. The objection is accepted in full

Bond dilution on a regular lattice is **statistically homogeneous at large
scales**. Diluting bonds shrinks the minimum cut but does not make it
non-planar: every cross-section remains statistically identical, so the cut
stays a plane and merely contains fewer surviving bonds. The audit's finding —
min cut = one full lattice cross-section, 36/25/16, **zero seed-to-seed
variance** — is a property of the *lattice*, and bond dilution does not leave
the lattice. My previous self-prompt proposed a fix that could not address the
defect it was aimed at. Accepted and superseded.

## 2. Committed choice: **Option B — real-data graph extraction**, with Option A as the named fallback

**Primary path: extract the Ni particle/neck graph directly from the
Holzer/Pecho FIB-SEM segmentation.** The data is present (`data/extracted`,
2.1 GB, all six stacks) and the extraction machinery is already built and
validated (`cmlib/pnm.py`, SNOW, with the serial-mode bug fixed and verified
against analytic cases in Project 1).

**Why B over A, stated as a decision rather than a preference:**

1. **Option A has no target to aim at.** Building a hard-sphere packing means
   choosing a coordination distribution, a contact-distance threshold, and a
   size distribution — and we have never measured what any of those are in a
   real Ni network. Building a disordered generator without those numbers would
   repeat the exact mistake this project already made twice: constructing a
   platform before measuring the thing it is supposed to reproduce. Step 0 fixed
   that for the YSZ side by measuring first; the Ni side needs the same.
2. **B answers the audit's open question directly.** The audit claims the planar
   minimum cut is the defect. That claim is only testable against a real Ni
   network's cut structure. If real Ni graphs also have near-planar cuts at
   comparable fractions, the diagnosis is wrong and Option A would have been
   built on a false premise.
3. **B removes every synthetic confound at once** — the ~3× size compression,
   the fat-neck proportion, the 8× TPB excess, the lattice periodicity, and the
   sphere-overlap channel — none of which can be fixed simultaneously in a
   generator without another full qualification cycle.
4. **Cost.** B reuses existing code and data. A is a new generator plus a full
   G1-a…G1-i re-qualification.

**Option C is explicitly declined**, including as a diagnostic probe. Its stated
purpose was to test "whether heterogeneity helps at all", but §1 establishes
that a diluted lattice does not produce the relevant heterogeneity, so a null
from C would be uninformative and a pass would be misleading. Spending the run
on B instead is strictly better.

**Fallback rule, frozen:** if B fails on REV or memory grounds (§5), fall back
to **Option A**, using the distributions measured in B as its explicit targets.
The sphere-overlap fix (`2R + 2·jitter·pitch < pitch`) remains mandatory for any
future synthetic Ni platform.

## 3. What is measured (frozen before running)

For each pristine real anode (fine / medium / coarse), on non-nested ROIs:

- **Ni graph** from `cmlib.pnm.extract_ni_network` (SNOW, serial, `sigma = 0.4`,
  `r_max = 4` — the Project-1 frozen parameters, unchanged). Nodes = Ni
  chambers, edges = throats, edge attribute = inscribed diameter and
  cross-sectional area.
- **Positional-disorder diagnostics:** pair-distance distribution between
  connected chamber centroids (must be checked for lattice-like sharp peaks);
  coordination-number distribution and its variance.
- **The frozen audit, unchanged:** unweighted S–T max-flow min-cut with virtual
  source/sink on the two faces of the transport axis, excluded from cuts;
  min-cut fraction of total edges; per-ROI variance; overlap between critical
  edges and the lower-quartile throat population; random and low-throat-area
  fragility curves.

**Transport axis** = axis 2 (x), the real-data convention from Phase 5 / Step 0.

## 4. Frozen decision criteria

The v3 platform (real graphs) is accepted as the basis for further damage work
if **all** hold:

1. **Pair-distance distribution is clearly non-lattice** — no sharp periodic
   peaks.
2. **Coordination-number distribution has non-trivial variance** (sd ≥ 1.0).
3. **Minimum-cut fraction has non-zero ROI-to-ROI variance** — unlike the
   lattice's exact zero.
4. **Minimum cut is not a full cross-section**, i.e. its edge count is
   materially below the number of edges crossing a mid-plane.
5. **`f_crit` ≤ 0.05 in at least one anode**, and **fine's `f_crit` ≤ coarse's**
   (the threshold-1 condition the lattice platform failed).

**If criteria 1–4 hold but 5 fails**, report that real Ni networks are also
robust to sparse neck removal, which would mean **percolation loss in real
anodes is not neck-mediated at all** — a substantive result that would redirect
the project away from every neck-based operator.

**If criteria 1–3 fail**, the SNOW graph is not capturing real disorder and the
extraction, not the platform, is at fault. Stop and report.

## 5. Scope, cost control, and honest limits

- ROIs sized as in Project 1 (8 µm cubes), 3 non-nested ROIs per anode, pristine
  only. Project 1's REV study found the **coarse anode is sub-REV at this ROI
  size** (Φ_Ni +13 % at 10 µm vs full stack, one ROI with no spanning cluster).
  That limitation carries over and is reported, not worked around; min-cut
  fractions for coarse are correspondingly the least reliable.
- Voxel validation of graph cuts is **not** performed in this pass. Removing a
  SNOW throat from a real segmentation has no unambiguous voxel realisation (no
  generator neck to zero out), so the graph result stands alone and is labelled
  as graph-level only. This is a real weakening relative to the synthetic audit
  and is stated rather than hidden.
- No damage operator is implemented. No O4, no O5v2.

## 6. Constraints

No parameter tuning. No new damage operators. `cmlib/damage.py`,
`cmlib/synth.py` unmodified. SNOW parameters frozen at Project-1 values. No
post-hoc redefinition of cut metrics — the audit metrics are those already
frozen in `PREREG_NI_VULNERABILITY_AUDIT.md` (`267efd3`).
