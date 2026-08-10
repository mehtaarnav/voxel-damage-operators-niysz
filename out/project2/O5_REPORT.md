# O5 — volume-conserving Ni agglomeration: **STOP. Operator mis-specified; not a valid test of the hypothesis.**

Run 2026-08-11 under `PREREG_AMENDMENT_O5.md`, frozen at `cb2ca49` **before**
any O5 execution. Code: `cmlib.damage2.apply_o5`. Data: `step2_o5_scope.csv`,
`step2_O5_ceiling_check.csv`. `cmlib/damage.py`, `cmlib/synth.py` unmodified.

---

## Verdict

**Two findings, and the second overrides the first.**

1. **O5 produces no Ni percolation loss** at any intensity in [1, 20], in 15 of
   15 structures (P_span 0.9988 / 0.9971 / 0.9929). Per the inherited boundary
   rule: stop and report, bracket not expanded.
2. **The frozen implementation does not realise the mechanism it was written to
   test.** Ni surface area **increases** by 3.3 % under O5. Agglomeration is by
   definition surface-area-reducing. **The operator roughens; it does not
   agglomerate.**

**Therefore this is not a falsification of the agglomeration hypothesis.** It is
a null on a mis-specified operator, and reporting it as a hypothesis null would
be a false claim. The frozen decision rule ("TPB retention < 0.50 ⇒ falsified")
is **not invoked**, because its precondition — a Ni-percolation transition at
which to measure TPB retention — never occurs, and because the operator failed a
validity check that the amendment did not think to impose.

---

## 1. What the diagnostics show

Averaged over 5 structure seeds per class, n = 20:

| analog | Ni P_span | volume error | **Ni surface ratio** | TPB ratio | Ni–YSZ interface ratio |
|---|---|---|---|---|---|
| fine | 0.9988 | **0.000000** | **1.014** | 2.02 | 0.676 |
| medium | 0.9971 | **0.000000** | **1.032** | 2.10 | 0.596 |
| coarse | 0.9929 | **0.000000** | **1.053** | 2.26 | 0.587 |

**What worked exactly as designed:** volume conservation is *exact* — zero voxels
gained or lost, in every structure, at every intensity. The conservation gate
(≤ 0.5 %) passes with room to spare, and this was the hard part of the
specification.

**What is wrong:** surface area should *fall* under agglomeration and instead
rises (ratio 1.03). TPB *doubles* rather than being retained-then-lost, and the
Ni–YSZ interface *drops* 33 % — the opposite pairing from the real anodes, where
TPB is retained. These three numbers together say the operator is producing a
rough, convoluted Ni surface, not compact agglomerates.

## 2. Why the unit test passed and the real structures did not

The unit test (two spheres joined by a thin free-standing neck) passed
convincingly: the neck went **63 → 0 voxels** at exactly conserved volume. That
is genuine Rayleigh-type behaviour, and it is why the operator was allowed to
run.

It generalised badly for a specific, identifiable reason. The curvature proxy is
a **uniform filter of the Ni mask** — i.e. a local Ni *volume fraction*. In an
isolated dumbbell, the only low-density region is the neck, so the proxy is
effectively a curvature measure. In a dense, three-phase, close-packed structure
the same quantity is dominated by **local packing**, not by local surface shape:
low-density surface sites are scattered over the whole structure rather than
concentrated at necks. Removing from scattered low-density sites and adding at
scattered high-density front sites produces protrusions and pits — roughening —
which is exactly what the surface-area ratio records.

**The proxy is valid only in the dilute limit.** That limitation was not visible
in the unit test and was not anticipated in the amendment.

## 3. The missing validity gate — the transferable lesson

The amendment froze a mechanism *description* ("matter flows from high mean
curvature to low"), a parameter, and a success criterion, but **no check that the
implementation realises the description**. A one-line gate would have caught
this before any bisection ran:

> **Agglomeration validity gate: total Ni surface area must DECREASE
> monotonically under the operator.** If it does not, the operator is not an
> agglomeration operator and no result from it may be interpreted.

This is the same class of error as Step 2's saturated `n_secondary`: a
pre-registered rule that was correct in intent and wrong in the specific value it
selected. Both are cheap to prevent and expensive to discover late.

## 4. What I did not do

- **I did not tune `MOVE_FRAC`.** It remains 0.05 as frozen.
- **I did not extend the bracket** past n = 20.
- **I did not silently swap the curvature proxy** and re-run. Changing the
  operator definition after seeing results is precisely the pattern the
  anti-tuning rules exist to prevent, even when the change is a defensible bug
  fix. That decision is the advisor's, not mine.
- **I did not report a hypothesis null.** The hypothesis remains untested.

## 5. Recommendation

**A corrected O5 requires a new amendment, not a patch.** The minimum change:

1. Replace the density proxy with a move rule that is **explicitly
   surface-area-reducing** — e.g. accept a candidate voxel move only if it
   lowers the total Ni surface-face count (a greedy or Metropolis surface
   minimisation at conserved volume). This makes "agglomeration" a property the
   operator *enforces* rather than one it is hoped to exhibit.
2. Add the §3 validity gate as a **pre-run** check on one structure per class.
3. Keep everything else frozen: same 15 structures, same bisection, same primary
   criterion (TPB retention ≥ 0.50 at the Ni transition), same falsification
   rule.

**Standing caveat that now applies to three operators.** O2 (severing every
lower-quartile neck), and O5 (redistributing 5 % of the Ni surface per round for
20 rounds) both failed to disconnect the Ni network, and O1 disconnected all
three classes at an identical intensity. The common factor is the **redundancy
of the regular jittered lattice**: with 6-connected topology and ~470–730 necks
per structure, local damage has many parallel paths to route around. Before
investing in another Ni operator, it is worth asking whether the *platform*, not
the operator, is what prevents Ni percolation from being lost the way it is in
the real fine anode. That question is answerable with the existing structures by
measuring how many necks must be removed, in decreasing-importance order, to
break spanning — a pure graph calculation requiring no new operator.

## 6. Status

- O5: **run complete, stop-and-report.** No transition in [1, 20]; operator
  fails the agglomeration validity check.
- Agglomeration hypothesis: **untested.**
- Project 2 currently stands at: C1 unresolved (O1), C2 pristine-loaded (O3 +
  R3 control), C3 failed, O2 null, O5 invalid.
- **Not run:** O4. No parameter tuned. `damage.py` / `synth.py` untouched.
