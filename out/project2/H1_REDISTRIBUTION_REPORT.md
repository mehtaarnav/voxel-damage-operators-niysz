# H1 — Ni redistribution: **UNTESTABLE on this dataset.** No support/falsify verdict issued.

Run 2026-08-11 under `PREREG_H1_REDISTRIBUTION.md` (frozen `5dda4ac` **before**
any H1 code was written). Slab thickness 0.5 µm, axis 2, orientation and layer
rules as frozen. Ordinal only; nothing fitted. Data:
`h1_depth_profiles.csv`, `h1_layer_summary.csv`.

## Verdict

**The pre-registered untestable branch fires.** H1 requires an
electrolyte-adjacent functional layer and a neighbouring support to be present
*in both stacks of a pair*. The gate — a sustained |ΔΦ_YSZ| ≥ 0.05 between the
two ends — is met in only one pair of three, **and not in the fine anode, which
is the anode H1 is about.**

| anode | pristine \|grad\| | degraded \|grad\| | pair evaluable? |
|---|---|---|---|
| **fine** | **0.0042** ✗ | **0.0240** ✗ | **NO** |
| medium | 0.0059 ✗ | 0.0552 ✓ | **NO** (pristine fails) |
| coarse | 0.1035 ✓ | 0.1157 ✓ | yes |

Clauses (a) and (b) cannot be evaluated for fine at all, so the ordinal clause
(c) — *fine > medium > coarse* — cannot be evaluated. **No support verdict and
no falsification verdict is issued.** Per the frozen rule this is not to be
worked around by redefining the interface post hoc, and it has not been.

## What the profiles actually show

Depth profiles were computed for all six stacks (38–48 slabs each, 0.5 µm).
Φ_Ni varies strongly with depth in every stack — e.g. coarse_post spans
**0.045–0.494** — but in fine and medium_pre that variation is **not organised
as a monotone interface gradient**. It is spatial heterogeneity, which is a
different thing and cannot be read as electrolyte-to-support transport.

The likeliest explanation is simply that these sub-volumes were sectioned within
a single layer. The stacks are 19–24 µm across the through-thickness axis, while
an anode support is hundreds of µm thick; there is no reason a given tomogram
must straddle the electrolyte interface, and five of six do not clearly do so.

## The one evaluable pair, reported descriptively only

For coarse — the sole pair passing the gate, and the anode for which H1 predicts
the *weakest* effect:

| state | functional-layer Φ_Ni | support Φ_Ni | FL − SUP |
|---|---|---|---|
| pristine | 0.1920 | 0.2648 | −0.0727 |
| degraded | 0.2046 | 0.2851 | −0.0805 |

The functional layer is Ni-poor relative to the support in both, by a similar
margin. **No paired ΔNi(x) is reported**, because pristine and degraded are
different specimens (Rx38 vs Rx41-3) sectioned independently: aligning x = 0
between them is arbitrary, and a difference computed across that alignment would
not mean what H1 needs it to mean. This was stated in the pre-registration
before the numbers existed.

## Limitations restated

**Pristine and degraded are different specimens.** Even had the interface been
present in all six stacks, ΔNi(x) would compare two different pieces of
material, and only the ordinal clause (c) would have carried any protection —
because a specimen artifact has no reason to order itself fine > medium >
coarse. That protection is unavailable here, since fine is the class that fails
the gate.

## What would make H1 testable

Not achievable with this dataset, and listed only so the boundary is explicit:
tomograms that deliberately straddle the electrolyte interface, and ideally the
same specimen imaged before and after — neither of which this dataset provides.
**H1 is neither supported nor rejected. It is unasked, because the data cannot
pose the question.**
