# Pre-registration — H1: electrolyte-to-support Ni redistribution

**Frozen 2026-08-11 before any H1 code is written or run.** Uses the existing
frozen codebase; `cmlib/damage.py`, `cmlib/synth.py` untouched. Ordinal only:
**no curve fitting, no diffusion coefficients, no kinetics.**

## H1
Fine anodes show greater electrolyte-to-support Ni redistribution than coarse:
(a) depth-dependent Ni depletion in the electrolyte-adjacent functional layer;
(b) corresponding Ni enrichment in the neighbouring support;
(c) ordinal strength fine > medium > coarse.

## Confound that must be stated first, not discovered later

**Pristine and degraded stacks are DIFFERENT SPECIMENS** — Rx36/Rx37/Rx38 versus
Rx41-1/Rx41-2/Rx41-3 (`phase2_volume_fractions.csv`). ΔNi(x) therefore compares
*two different pieces of material*, not the same volume before and after. Any
depth-dependent difference is confounded with specimen-to-specimen variation and
with where each specimen happened to be sectioned.

**Consequence, frozen:** a positive result can only ever be *suggestive* here,
never causal, and the ordinal requirement (c) is the only part that is even
partially protected — because a specimen artifact has no reason to order itself
fine > medium > coarse. **This is recorded as a limitation on the headline
before the result is known.**

## Protocol (frozen)

1. **Through-thickness axis = array axis 2 (x)**, the Phase-5/Step-0 transport
   convention.
2. **Slab thickness = 0.5 µm**, chosen now, applied to every stack. Slabs are
   full cross-sections; the last partial slab is discarded.
3. **Orientation.** x = 0 is placed at the **electrolyte-adjacent end**,
   identified objectively as the end whose first 10 % of slabs has the **higher
   mean Φ_YSZ**. The choice made for each stack is reported.
4. **Layer identification.** Functional layer = the electrolyte-adjacent region
   over which Φ_YSZ exceeds its stack median; support = the remainder. Reported
   per stack. If Φ_YSZ shows no sustained gradient, see the untestable branch.
5. **ΔNi(x) = Φ_Ni,degraded(x) − Φ_Ni,pristine(x)** per matched anode, on the
   common slab index range.
6. **ΔNi_connected(x)** computed on the spanning Ni cluster only, same slabs
   (this is the transport-relevant quantity and is reported alongside).

## Untestable branch — pre-registered

**If a stack does not contain an electrolyte/electrode interface**, x = 0 is
undefined and (a)/(b) cannot be evaluated. Operational test: the stack must show
a **sustained monotone Φ_YSZ gradient of at least 0.05 in absolute volume
fraction** between its two ends. If it does not, H1 is reported as **UNTESTABLE
on this dataset**, with the profiles shown, and no support/falsify verdict is
issued. **This is a legitimate outcome and is not to be worked around** by
redefining the interface post hoc.

## Decision rules

- **SUPPORT H1** if (a) fine shows the strongest electrolyte-adjacent Ni
  depletion, (b) fine shows the strongest support Ni enrichment, and (c) the
  effect is spatially coherent (same sign across adjacent slabs, not a
  single-slab spike).
- **FALSIFY H1** if fine does not show the strongest redistribution, or if
  depletion appears with no corresponding support enrichment.
- Both outcomes are reported. A falsification is reported as *"the Ni
  redistribution hypothesis is rejected for this dataset"*, with the
  different-specimen confound restated.

## Constraints
Ordinal comparison only. No fitting of any kind. Existing frozen codebase.
