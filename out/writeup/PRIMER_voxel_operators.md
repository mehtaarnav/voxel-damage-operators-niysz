# When the validity check destroys the physics

### A primer on voxel damage operators, Ni–YSZ anode degradation, and how a simulation can pass every test you set it and still be wrong

---

## Scope of this document

This is a teaching document. It assumes you know what a crystal, a grain
boundary and a diffusion coefficient are, and that you can read code. It does
not assume you have met solid oxide fuel cells, tomography, percolation theory,
Monte Carlo methods or phase-field models.

**What it covers.** The physical system and why it degrades; the experimental
puzzle that started the project; how a voxel damage operator works and what
metric you must use to judge one; the derivation of the area-change identity
that governs single-voxel swaps; the central result, which is that an operator
can satisfy a surface-area validity criterion in full while corrupting the
quantity you actually care about; and a separate treatment of why finite
temperature does not rescue the method.

**What it defers.** The full experimental protocol, the segmentation and
tomography details, the statistical pre-registration, and four of the six
failure modes in their derived form are in the manuscript
(`out/writeup/manuscript.tex`). Those four get a paragraph each here so you know
they exist and where to look.

**One thing to hold onto while reading.** The surface-area criterion discussed
below is not a strawman. It is the natural, defensible, physically motivated
thing to check. It is what I would have written. That is the entire point.

---

## 1. The physical system

### 1.1 Solid oxide fuel cells, briefly

A fuel cell converts chemical energy to electrical energy without burning
anything. A solid oxide fuel cell (SOFC) does this at high temperature —
typically 700–1000 °C — using a ceramic electrolyte that conducts oxide ions
(O²⁻) but not electrons.

Three layers matter:

- **Cathode.** Air arrives. Oxygen molecules pick up electrons and become oxide
  ions: O₂ + 4e⁻ → 2O²⁻.
- **Electrolyte.** A dense ceramic, usually yttria-stabilised zirconia (YSZ),
  which lets O²⁻ through and blocks electrons. Blocking electrons is what forces
  them through the external circuit, which is where you get your current.
- **Anode.** Fuel arrives, typically hydrogen. Oxide ions arriving through the
  electrolyte react with it: H₂ + O²⁻ → H₂O + 2e⁻.

Run the same device backwards and it becomes a solid oxide electrolysis cell
(SOEC), splitting steam into hydrogen. The degradation physics is closely
related, which is why the literature talks about "solid oxide cells" generally.

### 1.2 Why the anode is a composite

Look at that anode reaction again:

> H₂ + O²⁻ → H₂O + 2e⁻

For it to happen at a given point in space, three things must be simultaneously
present:

1. **Hydrogen gas** — so there must be open pore connected to the fuel supply.
2. **An oxide ion** — so there must be YSZ connected to the electrolyte.
3. **Somewhere for the electrons to go** — so there must be nickel connected to
   the current collector.

Nickel is the electronic conductor and the catalyst. YSZ is the ionic
conductor. Pore is the gas path. The reaction can only occur where all three
meet.

That locus is a **one-dimensional line** in three-dimensional space: the curve
along which the Ni surface, the YSZ surface and the pore all touch. It is called
the **triple-phase boundary (TPB)**. Its total length per unit volume — units of
µm⁻², i.e. length per volume — is the single most important microstructural
descriptor of an SOFC anode, because it is proportional to the number of
available reaction sites.

A real Ni–YSZ anode is therefore a three-phase interpenetrating composite:
roughly 25–32 % nickel, 25–40 % YSZ, the rest porosity, with all three phases
percolating.

### 1.3 Percolation

**Percolation** means: does a phase form a connected path all the way across the
sample?

This is not a subtlety. If the nickel network is broken, electrons produced deep
in the anode cannot reach the current collector, and that region is
electrochemically dead no matter how much TPB it contains. Ni percolation is a
hard requirement for the anode to function.

Two ways to measure it, and the difference matters later:

- **P_span** — the fraction of Ni voxels belonging to a cluster that touches
  *both* opposite faces of the volume along a chosen axis. A spanning path.
- **P_reach** — the fraction of Ni voxels reachable from *one* face.

P_reach ≥ P_span always, because spanning is a stronger condition. They can rank
samples differently, and in this dataset they do.

### 1.4 Why it degrades: nickel coarsening

Nickel at 950 °C is mobile. It has high surface energy (γ_Ni ≈ 2 J/m²) and,
critically, it **does not wet YSZ well** — the contact angle is obtuse, around
130° for an as-prepared electrode
([Björnsson et al. 2026](https://doi.org/10.1016/j.jpowsour.2025.239069)),
superseding the widely cited 120° from Trini et al. An obtuse contact angle means
nickel would rather ball up than spread out on the ceramic.

Left at temperature, the system does what all high-surface-energy systems do: it
reduces total interfacial energy. Nickel migrates by surface diffusion from
regions of high curvature to regions of low curvature. Small particles shrink,
large particles grow (Ostwald ripening). Narrow necks between particles can pinch
off. The nickel phase becomes coarser and more compact.

The consequences for cell performance:

- **TPB is lost.** Coarser nickel means less nickel surface, and TPB scales
  roughly with the amount of Ni/YSZ/pore junction available.
- **Percolation can be lost.** If enough necks pinch off, the nickel network
  fragments.

This is a well-documented, well-studied degradation mode
(Simwonis et al. 2000; Monaco et al. 2019). Redox cycling — accidentally
oxidising the nickel to NiO and reducing it back — accelerates it and also
cracks the YSZ scaffold (Sarantaridis & Atkinson 2007).

### 1.5 Rayleigh break-up, which will matter a great deal

There is a specific mechanism by which a neck fails. A solid cylinder connecting
two particles is **unstable** to perturbations whose wavelength exceeds its
circumference — the solid-state analogue of a water stream breaking into
droplets. This is the Rayleigh–Plateau instability, treated for solid-state
surface diffusion by Nichols and Mullins (Nichols 1965).

The feature to notice now, because the whole argument later turns on it: the
instability is **collective and long-wavelength**. It lowers energy only once a
correlated set of atoms has moved together. Any *small local* perturbation to a
straight cylinder increases surface area before the instability can pay it back.
The system must go uphill in energy before it goes downhill.

Hold that thought.

---

## 2. The central paradox

Three real Ni–YSZ anodes, from a published coarseness series
(Pecho et al. 2015), differing in nickel particle size:

| | mean Ni particle size | Φ_Ni | pristine TPB density |
|---|---|---|---|
| **fine** | 1148 nm | 0.322 | 3.624 µm⁻² |
| **medium** | 1445 nm | 0.250 | 2.109 µm⁻² |
| **coarse** | 1715 nm | 0.229 | 1.473 µm⁻² |

Each was subjected to the same degradation protocol and re-imaged by FIB-SEM
tomography. Measuring what fraction of each property survived:

| | Ni percolation retained (P_span) | TPB retained |
|---|---|---|
| **fine** | **0.680** | **0.799** |
| **medium** | 0.855 | 0.746 |
| **coarse** | 0.947 | 0.590 |

Read those two columns against each other. **The fine anode is the worst at
retaining its nickel network and the best at retaining its reaction sites.** The
coarse anode is the reverse.

That is strange. Both properties are supposed to be destroyed by the same
process. If nickel coarsens, you expect to lose connectivity *and* TPB together.
Instead they are anticorrelated across the series.

### 2.1 Three caveats you must carry, or you will draw a false conclusion

This is a teaching document, so here is the part that gets skipped in seminars.

**(a) Pristine and degraded are different specimens.** These are not one sample
imaged before and after. Pristine is Rx36/37/38; degraded is Rx41-1/2/3.
Different physical pieces of material. Every "retention" number above is a ratio
between two different specimens of nominally the same type. This is a real
limitation and it is why nothing in the project claims a kinetic rate.

**(b) Only two of the six orderings are defensible.** The clean three-level
ranking in the table is not established:

- Ni retention **flips with the metric**. P_span gives 0.680 / 0.855 / 0.947
  (coarse best). P_reach gives 0.857 / 0.979 / 0.942 (medium best).
- The **published** TPB retention series is **non-monotone**: 0.7434 / 0.5862 /
  0.6075 — coarse exceeds medium.

What survives both checks is exactly two statements:

> **Fine retains nickel percolation worst. Fine retains TPB best.**

Medium versus coarse is unresolved on both axes and should be treated as noise.
Any model that is scored on reproducing a three-level ranking can be *failed for
getting the physics right*, since the physics may well produce the non-monotone
published ordering.

**(c) Three samples is not a statistical result.** There are 3! = 6 possible
orderings. Since medium-vs-coarse is unresolved, the effective bar for a
predictor is "put fine last," which 2 of 6 orderings clear **by chance** — 33 %.
No correlation coefficient or p-value appears anywhere in this project, and none
should be computed from three points.

### 2.2 The confound that killed the original hypothesis

The project began by asking whether a graph-theoretic "connectivity margin"
computed on the pristine network predicts retained percolation better than plain
coarseness does.

It does not, and the reason is instructive. In this dataset coarseness varies
monotonically across all three samples and drags every candidate metric along
with it. The metrics that reproduce the outcome ordering — 10th-percentile neck
width (69.5 / 156.9 / 205.3 nm), algebraic connectivity λ₂ — produce *exactly the
same ordering as mean particle size*, so on three samples they are
indistinguishable from a coarseness proxy by construction. The metrics that are
genuinely topological rather than coarseness proxies — minimum-cut conductance,
effective network conductance — fail outright and match no outcome ordering.

Worse, λ₂ turned out to be measuring node count rather than connectivity: the
raw values span 7.65× across the anodes, but multiplying by node count collapses
that to 1.24× (232 / 200 / 247). Almost everything λ₂ appeared to say was a
restatement of how many watershed chambers each structure contains.

To separate a connectivity margin from a coarseness proxy you would need samples
where the two **disagree** — matched particle size and volume fraction, differing
neck-width distributions. Three samples that differ in everything at once cannot
do it, however carefully the graph metrics are computed.

That negative result is what motivated the rest: if a static predictor cannot
explain the anticorrelation, simulate the degradation process itself.

---

## 3. Simulating degradation on a tomogram

### 3.1 What the data actually is

FIB-SEM tomography gives you a stack of 2D images that, once segmented, becomes a
3D array where every voxel carries a label: Ni, YSZ, or pore. That is it. A big
integer array.

For this project:

| anode | ROI | voxels | voxel size |
|---|---|---|---|
| fine | 400 × 410 × 410 | 67.2 M | 20.0 × 19.5 × 19.5 nm |
| medium | 320 × 328 × 328 | 34.4 M | 25.0 × 24.4 × 24.4 nm |
| coarse | 400 × 412 × 412 | 67.9 M | 30.0 × 29.1 × 29.1 nm |

Note that the coarse ROI is 12 µm on a side while fine and medium are 8 µm. The
coarse structure is too coarse for an 8 µm cube to be representative — an 8 µm
coarse region yields only 48–72 network nodes, which was judged uninterpretable.
This size difference is a stated confound, not something to paper over.

### 3.2 Two ways to simulate

**Phase-field** treats the microstructure as continuous order-parameter fields
and evolves a free-energy functional by solving partial differential equations.
Interfaces are diffuse — smeared over several grid points. It handles topology
changes naturally and has an explicit, physical free energy including surface
energies and contact angle. It is the established approach for Ni–YSZ coarsening
(Chen et al. 2011; Trini et al. 2021; Yang et al. 2023; Hoffrogge et al. 2023).
It is also expensive: fourth-order dynamics on tens of millions of grid points.

**Voxel damage operators** work directly on the label array. Pick voxels, change
their labels according to some rule, repeat. Cheap, conceptually simple, and it
operates on exactly the data structure the experimentalist already has. There is
published precedent on the fabrication side: Zhou et al. (2023) built a Potts
kinetic Monte Carlo model of NiO–YSZ *sintering*, calibrated against FIB-SEM
data.

The attraction of the second approach is obvious, and this project took it. What
follows is an account of what that cost.

**A note on names, because the literature is loose here.** "Voxel-swap operator"
describes the move class: exchange labels between two voxels. It says nothing
about how you decide whether to accept a move. Two acceptance rules appear
below:

- **Greedy (zero-temperature)** — accept a move only if it does not increase
  some cost. Deterministic in spirit.
- **Metropolis (finite-temperature)** — accept improving moves always, and
  accept worsening moves with probability exp(−ΔE/kT). This is what "kinetic
  Monte Carlo" properly means.

They are different rules on the same move class, and conflating them makes
Section 7 unintelligible.

### 3.3 The metric problem, and why R_Ni exists

Before you can judge an operator you need a number that measures "how much
connected nickel is left." The obvious choice — P_span, the spanning fraction —
has a trap in it, and the trap is worth understanding because it is a general
one.

Real electrodes contain isolated nickel: little islands disconnected from the
main network. In this dataset, **1.2 % to 11.2 %** of the nickel is already
disconnected before any degradation. Per-ROI spanning fractions are 0.9821 /
0.9754 / 0.9877 for fine, 0.9713 / 0.9446 / 0.9554 for medium, 0.8878 / 0.9528 /
0.9157 for coarse.

Many operators include a housekeeping step: after damaging the structure, delete
the non-spanning clusters, because disconnected metal is electrochemically dead
anyway. This seems harmless. It is not.

P_span is defined as

$$P_{\text{span}} = \frac{\text{Ni voxels in the spanning cluster}}{\text{Ni voxels present now}}$$

If the operator deletes everything that is not in the spanning cluster, then the
numerator and denominator become the same set, and

$$P_{\text{span}} = 1.0000 \quad \text{identically, forever.}$$

Two operators in this project did exactly this and drove pristine P_span to
exactly 1.0000 at the first step. **The operator rewrote the denominator of the
metric being used to judge it.** The metric stopped measuring the structure and
started measuring the housekeeping.

The fix is to fix the denominator to something the operator cannot touch:

$$\boxed{R_{\text{Ni}}(n) = \frac{\text{Ni voxels in the spanning cluster after } n \text{ rounds}}{\text{Ni voxels in the } \textbf{pristine} \text{ structure}}}$$

Two properties, and one important non-property.

**Invariant to pruning.** Deleting isolated clusters changes neither the spanning
cluster (they were not in it) nor the pristine denominator (it is fixed). So
R_Ni is unchanged by pruning, and the circularity is gone.

**Monotone under removal operators.** If an operator only ever removes Ni voxels,
the spanning cluster can only shrink or fragment, so R_Ni is non-increasing. This
gives a meaningful damage curve.

**Not monotone in general.** This is not a technicality — we measured it. Under a
*volume-conserving swap* operator, R_Ni **rises**: fine 0.9821 → 0.9825, coarse
0.8878 → 0.8881. Moving a voxel from one place to another can bridge a gap and
connect previously isolated nickel. So "R_Ni is monotone" is true of removal and
false of rearrangement, and the difference is itself one of the findings below.

A mandatory sanity check comes free: at n = 0 the spanning cluster and the
pristine structure are the same, so **R_Ni(0) must equal pristine P_span
exactly**. If it does not, the implementation is wrong. Run that check before
anything else.

---

## 4. The mathematics of a single-voxel swap

### 4.1 Setting up

Work on a cubic lattice with **6-connectivity**: each voxel has six face
neighbours (±x, ±y, ±z). Diagonal contacts do not count as touching.

Model surface energy the simplest defensible way: **the exposed area of the
nickel phase is the number of faces where a Ni voxel touches a non-Ni voxel.**
Call that total $A$, measured in face units. The physical surface energy is
$\gamma_{\text{Ni}} A a^2$ where $a$ is the voxel edge, but for ranking moves the
face count is enough.

Coarsening is volume-conserving — nickel moves, it does not vanish. So the move
class is a **swap**: remove one Ni voxel at site $a$, add one Ni voxel at a pore
site $b$. Exactly one out, one in, so $\Delta V = 0$ by construction.

Let $n_N(x)$ be the number of Ni 6-neighbours of voxel $x$. **The voxel itself is
not counted.** (This sounds pedantic. It is not — see Appendix A.)

### 4.2 The area-change identity

Consider removing the Ni voxel at $a$. It has $n_N(a)$ Ni neighbours and
$6 - n_N(a)$ non-Ni neighbours.

- Its $6 - n_N(a)$ exposed faces disappear: area falls by $6 - n_N(a)$.
- Each of its $n_N(a)$ Ni neighbours was previously buried against it and is now
  exposed: area rises by $n_N(a)$.

Net change from removal:

$$\Delta A_{\text{remove}}(a) = n_N(a) - \bigl(6 - n_N(a)\bigr) = 2n_N(a) - 6$$

Now add a Ni voxel at pore site $b$ with $n_N(b)$ Ni neighbours, by the mirror
argument:

$$\Delta A_{\text{add}}(b) = \bigl(6 - n_N(b)\bigr) - n_N(b) = 6 - 2n_N(b)$$

Provided $a$ and $b$ are **not adjacent** — if they were, the bond between them
would be double-counted — the two contributions add:

$$\boxed{\Delta A = 2\bigl[n_N(a) - n_N(b)\bigr]}$$

This is exact on a 6-connected lattice, not an approximation. It has been
verified numerically against brute-force enumeration of every bond in the
structure, to zero error.

**What it means.** The entire surface-area criterion reduces to comparing two
integers:

$$\Delta A \leq 0 \iff n_N(a) \leq n_N(b)$$

In words: *remove nickel from where it has few nickel neighbours (a convex bump,
high curvature) and put it where it has many (a concave notch, low curvature).*
That is precisely curvature-driven surface diffusion, which is the physics we
wanted. The identity is a genuine success — it converts a continuum idea into an
exact lattice statement.

Two preconditions worth stating because they were learned the hard way:

1. **$a$ and $b$ must be non-adjacent**, else the shared bond is double-counted.
2. **Both must be interior voxels.** The derivation substitutes "6 neighbours";
   a voxel on the domain boundary has fewer, and the identity fails there.

### 4.3 The validity criterion, and why it is reasonable

How do you know your coarsening operator is actually coarsening?

Coarsening reduces surface area. That is what coarsening *is*, thermodynamically:
the system minimising interfacial energy at fixed volume. So the natural check is

> **Gate: specific surface area must decrease monotonically.**

with specific surface area $S_{\text{spec}} = A / V_{\text{Ni}}$, i.e. exposed
faces per nickel voxel.

Pause here. This criterion is correct physics, it is cheap to evaluate, it is
falsifiable, and it directly encodes the thermodynamic driving force. If you were
designing this study you would impose it. It was pre-registered before any run.

Everything that follows is about what that criterion fails to constrain.

### 4.4 The first consequence: the area barrier

Recall §1.5: Rayleigh break-up requires the system to go **uphill** before it
goes downhill. On the lattice, that means some accepted move must have
$\Delta A > 0$.

But the gate forbids exactly those moves.

So:

> **Any single-voxel-swap operator restricted to monotone area reduction cannot
> produce Rayleigh-type neck break-up.** The criterion and the mechanism cannot
> both be honoured.

This is an impossibility statement about the operator class, not a bug report
about an implementation. And it bites: any published study that both enforces
monotone surface-area reduction and claims to reproduce neck break-up is either
not producing the break-up or is violating its own validity criterion.

Empirically, on a model dumbbell — two spheres joined by a straight cylinder —
this bracketing is visible:

| rule | neck volume 63 → | S_spec (pristine 0.45052) |
|---|---|---|
| curvature-ranked, 26-connectivity | **63 → 57 → 15 → 0** | rises to 0.45195 at n=1 — **violates the gate** |
| greedy, ΔA ≤ 0 | 63 → 63 (unchanged) | satisfies the gate |

The operator that thins the neck breaks the criterion at its very first step. The
operator that respects the criterion never thins the neck. **Thinning and the
criterion were never observed together.**

Be precise about the quantifier: it is *not* true that "no operator thinned a
neck." One did. It failed the gate to do so.

**One honest qualification.** The argument above is airtight for a *strictly*
decreasing criterion, $\Delta A < 0$. Most implementations use the non-strict
form $\Delta A \leq 0$, which admits moves with $\Delta A = 0$ exactly. Those
moves form a large flat region — a **plateau** — that a monotone rule is free to
wander. That such wandering never reaches break-up is established here
empirically, over some 2.6 × 10⁵ accepted moves, not proved analytically. Keep
that distinction; it is the difference between a theorem and a strong
measurement.

---

## 5. What the operator actually does on real electrodes

Now run the greedy area-decreasing swap — accept iff $\Delta A \leq 0$, one move
at a time, neighbour field recomputed after each accepted move — on the three
real ROIs.

Volume is conserved exactly. YSZ is never touched. And:

| | S_spec before → after | ΔS | TPB before → after | TPB ratio | R_Ni before → after |
|---|---|---|---|---|---|
| **fine** | 0.15696 → 0.15136 | **−0.00560** | 4.4774 → 22.5346 | **5.03×** | 0.9821 → **0.9825** |
| **medium** | 0.11670 → 0.11419 | **−0.00251** | 1.8660 → 9.2329 | **4.95×** | 0.9713 → **0.9714** |
| **coarse** | 0.13991 → 0.13524 | **−0.00468** | 1.5365 → 5.6862 | **3.70×** | 0.8878 → **0.8881** |

Read the second column, then the fourth.

**The gate passes.** Specific surface area falls, strictly and monotonically, on
all three anodes, at exact volume conservation. By the criterion we set — the
correct, physically motivated criterion — this operator is behaving like a
coarsening process.

**The physics is destroyed.** Triple-phase boundary density, the quantity the
entire device depends on, **increases by a factor of 3.7 to 5.0.** Real nickel
coarsening *destroys* TPB. This operator manufactures it, wholesale, while
satisfying its validity check.

And R_Ni *rises*. Percolation is not even monotone under the operator: a swap can
bridge a gap and reconnect isolated nickel.

The neck, meanwhile, never thins.

This is the result. An operator can satisfy a correct, pre-registered,
thermodynamically motivated validity criterion in full, and simultaneously move
the microstructure in the opposite direction from the physics on the one metric
that matters most.

### 5.1 Why: the causal chain

Step through it.

**(1) The gate constrains one quantity: nickel/pore surface area.** That is what
$\Delta A = 2[n_N(a) - n_N(b)]$ measures. It counts faces on the Ni/pore
interface.

**(2) TPB is a different quantity.** It is a *one-dimensional junction* count —
the length of curve where Ni, YSZ and pore all meet. It is not a facet count on
the Ni/pore interface. Nothing in $\Delta A$ mentions YSZ at all. The gate is
blind to the three-phase junction.

**(3) A large fraction of accepted moves are exactly area-neutral.** Because
$\Delta A = 2[n_N(a) - n_N(b)]$ takes integer values and the criterion is
non-strict, every move with $n_N(a) = n_N(b)$ has $\Delta A = 0$ and is accepted.
Measured on the real ROIs:

| anode | accepted moves | of which ΔA = 0 | fraction |
|---|---|---|---|
| fine | 258 870 | 208 180 | **80.4 %** |
| medium | 87 945 | 76 152 | **86.6 %** |
| coarse | 175 055 | 141 572 | **80.9 %** |

Roughly four out of five accepted moves **cost nothing** under the criterion.

**(4) Unpriced moves are unconstrained moves.** A move with $\Delta A = 0$ is
invisible to the gate. The operator is free to make as many as it likes, in any
configuration, and the validity check will never object. This is a *null space*
of the criterion.

**(5) The null space scatters nickel across the junction.** Nothing steers those
moves. They redistribute nickel voxels essentially at random along the surface,
and every newly created Ni/pore facet that happens to sit next to YSZ adds triple
line. The junction count inflates.

That is the mechanism. **The gate prices surface area; TPB lives in the
directions the gate does not price.**

### 5.2 The obvious fix, and why it is not enough

If area-neutral moves are the problem, forbid them. Tighten the criterion to
strict inequality, $\Delta A < 0$, so every accepted move must *strictly* lower
the area.

| anode | accepted moves | ΔS | TPB ratio |
|---|---|---|---|
| fine | 31 368 | −0.00370 | **1.527×** |
| medium | 6 700 | −0.00147 | **1.325×** |
| coarse | 22 249 | −0.00329 | **1.340×** |

Inflation drops from 3.7–5.0× to **1.33–1.53×**. Large improvement. But it does
not vanish, and it is not close to 1.0.

So the plateau is the *dominant* driver, not the *exclusive* one. Even when every
single move strictly lowers surface area, the operator still manufactures a third
to a half again as much triple line as the electrode started with. R_Ni still
rises. The neck still does not thin.

**Tightening the criterion does not repair it.** And that forces the real
conclusion, which is more general than "the gate has a loophole":

> **Surface area is the wrong invariant to validate a coarsening rule against.**
> A rule can satisfy it in full — strictly, monotonically, at exact volume
> conservation — and still move the microstructure away from the physics it is
> meant to represent. The problem is not that the budget admits a null space. It
> is that the budget measures the wrong thing.

If you take one idea from this document, take that one.

### 5.3 A practical check that costs one line

An operator whose TPB retention **exceeds unity** is roughening, not coarsening.

Real degradation destroys TPB. If your simulated TPB goes up, your operator is
adding voxel-scale texture at the three-phase junction, and any TPB conclusion
drawn from it is reporting discretisation rather than the electrode. It costs one
line to check, and in this project it disqualified three operators.

Note that this bias cannot be removed by taking ratios to the pristine state or
by calibrating before damage. That voxelised TPB estimates carry a
resolution-dependent systematic bias is established (Jørgensen & Bowen 2014;
Shimura et al. 2016). What is different here is that the bias is a large,
non-monotone excursion *driven by the operator itself*, so it moves during the
run.

---

## 6. The other four failure modes

The area barrier (§4.4) and TPB manufacture (§5) get full treatment above. Four
more were measured. Each gets a paragraph; the derivations are in the manuscript.

**Largest-component pruning fixes the metric at unity.** Covered in §3.3, because
you needed it to understand R_Ni. Restated as a failure mode: an operator that
deletes non-spanning clusters drives P_span to exactly 1.0000 at the first damage
round, because it has made the numerator and denominator the same set. The
operator rewrites the metric used to judge it. Two operators in this project did
this. The general lesson: never let the operator modify the denominator of your
outcome metric.

**Qualification gates demanding perfect pristine connectivity are
unrepresentative.** It is tempting to require that your test structure be fully
connected before damage — it makes the initial condition clean. But real
electrodes carry 1.2–11.2 % disconnected nickel. A synthetic structure built to
be perfectly connected is not a harder test case, it is a *different* test case,
and operators tuned on it behave differently on real data. Worse, if the pristine
structure is forced to P_span = 1.0000 by the gate, the pruning artifact above
becomes invisible, because pristine was already at unity.

**A jittered regular lattice has a minimum cut that is a plane.** Synthetic
microstructures are often built as a perturbed lattice of spheres joined by
necks. The trouble: because every cross-section of such a lattice is
statistically identical, the minimum cut through the nickel network is *exactly
one full cross-section* — 36, 25 or 16 throats depending on the lattice — with
**zero seed-to-seed variance**. Real networks fail through small non-planar cuts
at 0.5–3 % of throats (min-cut fractions 0.0055–0.0067 for coarse, 0.0141–0.0196
for fine). Any conclusion about how a network fails, drawn on a lattice, is a
conclusion about the lattice.

**Curvature rank is not the same as area change.** A natural shortcut is to rank
candidate sites by a curvature proxy and move the most convex nickel to the most
concave pore site, on the reasoning that this is what curvature-driven diffusion
does. It is not equivalent to computing $\Delta A$. On the model dumbbell the
curvature-ranked operator *raised* specific surface area at the first step —
0.45195 against a pristine 0.45052 — before falling below pristine at n = 3. If
the ranking is a proxy for the thing you are gating on, you have no guarantee at
all. Equation §4.2 was derived precisely to replace the proxy with the exact
quantity.

**And a fifth, closely related to §5: erosion manufactures TPB even more
violently.** Stochastic single-voxel *removal* — a pure erosion operator, not a
swap — multiplies TPB density by **14.8, 12.3 and 7.4** at the first damage round
before collapsing it to 4.8–10.2 % of pristine by the eighth. Same mechanism as
§5: pitting the nickel surface at the voxel scale creates enormous amounts of new
junction. The swap result in §5 shows this is not specific to removal; a
volume-conserving, area-*decreasing* rule does it too, just less violently.

---

## 7. Why not just use finite temperature?

Everything in §4.4 rests on one thing: the gate forbids $\Delta A > 0$ moves, and
Rayleigh break-up needs them.

So there is an obvious objection, and it is a good one:

> *Do not use a greedy rule. Use Metropolis acceptance — a proper kinetic Monte
> Carlo. Accept improving moves always, and accept area-increasing moves with
> probability exp(−ΔE/kT). Thermal fluctuation carries the system over the
> barrier, exactly as it does in the real material at 950 °C. Problem solved.*

This section answers that objection. **It is a separate line of investigation
from §§4–6, and its conclusions stand on their own.**

### 7.1 A three-phase energy with a contact angle

To do this properly you need YSZ in the energy, not just as an obstacle. Write a
Potts-style bond energy over 6-connected neighbours:

$$E = \sum_{\langle ij \rangle} J(\sigma_i, \sigma_j), \qquad \sigma \in \{\text{Ni}, \text{YSZ}, \text{pore}\}$$

with same-phase bonds costing zero and three distinct bond energies:
$J_{NP}$ (Ni–pore), $J_{NY}$ (Ni–YSZ), $J_{YP}$ (YSZ–pore).

For the same volume-conserving swap — Ni at $a$ becomes pore, pore at $b$ becomes
Ni — count bonds as in §4.2 but now tracking both neighbour types. Writing
$n_P = 6 - n_N - n_Y$ and adding the two site contributions:

$$\Delta E = 2 J_{NP}\bigl[n_N(a) - n_N(b)\bigr] + \bigl[n_Y(a) - n_Y(b)\bigr]\bigl(J_{NP} + J_{YP} - J_{NY}\bigr)$$

Now bring in the physics. **Young's construction** relates the contact angle to
the three interfacial energies. For a droplet on a surface, the horizontal force
balance at the contact line gives $\gamma_{\text{sv}} = \gamma_{\text{sl}} +
\gamma_{\text{lv}}\cos\theta$. Translating to our three solid/vapour phases, with
bond energies standing in for interfacial energies:

$$\cos\theta = \frac{J_{YP} - J_{NY}}{J_{NP}}$$

so $J_{NP} + J_{YP} - J_{NY} = J_{NP}(1 + \cos\theta)$, and:

$$\boxed{\Delta E = J_{NP}\Bigl[\,2\bigl(n_N(a) - n_N(b)\bigr) + (1 + \cos\theta)\bigl(n_Y(a) - n_Y(b)\bigr)\Bigr]}$$

This has been verified numerically against brute-force bond enumeration, to zero
error, at θ = 0°, 90°, 130° and 180°. Check the limits, which is how you know a
derivation is right:

- **No YSZ present** ($n_Y = 0$): reduces exactly to $\Delta E = 2J_{NP}[n_N(a) -
  n_N(b)]$, the two-phase identity of §4.2. ✓
- **θ = 180°, perfect dewetting** ($\cos\theta = -1$): the YSZ term vanishes
  entirely, so YSZ becomes indistinguishable from pore as far as nickel is
  concerned. Correct — at 180° nickel does not wet YSZ at all. ✓
- **θ = 0°, perfect wetting** ($\cos\theta = +1$): the YSZ coefficient becomes 2,
  identical to the Ni coefficient, so a YSZ neighbour is exactly as favourable as
  a Ni neighbour. Correct. ✓

Same precondition as before: interior voxels only, since $n_P = 6 - n_N - n_Y$
assumes exactly six neighbours.

**An aside worth internalising.** If you specify $\gamma_{\text{Ni}}$,
$\gamma_{\text{Ni-YSZ}}$ *and* θ independently, you have over-determined the
system — Young's equation ties them. Taking commonly quoted values
$\gamma_{\text{Ni}} = 2.0$ J/m², $\gamma_{\text{Ni-YSZ}} = 1.5$ J/m² and θ = 130°
implies $\gamma_{\text{YSZ}} = 1.5 - 2.0\cos(130°) = 0.21$ J/m², roughly six
times below the literature value. Freeze two and derive the third; never all
three.

### 7.2 The athermal barrier

Now put a temperature in it. Map the bond energy to physical units:

$$J_{NP} = \gamma_{\text{Ni}} \times a^2$$

where $a$ is the voxel edge. At tomographic resolution $a \approx 20$ nm, so the
face area is $4 \times 10^{-16}$ m², and with $\gamma_{\text{Ni}} = 2.0$ J/m²:

$$J_{NP} = 2.0 \times (20 \times 10^{-9})^2 = 8.0 \times 10^{-16}\ \text{J}$$

Thermal energy at 950 °C = 1223 K:

$$k_B T = 1.381 \times 10^{-23} \times 1223 = 1.69 \times 10^{-20}\ \text{J}$$

Take the ratio:

$$\frac{J_{NP}}{k_B T} = \frac{8.0 \times 10^{-16}}{1.69 \times 10^{-20}} \approx 4.7 \times 10^{4}$$

**A single voxel face carries roughly forty-seven thousand kT of surface
energy.**

The Boltzmann factor for the smallest barrier-crossing move follows. Taking the
smallest three-phase barrier, $(1+\cos\theta) = 0.357$ at θ = 130°:

$$P = \exp\left(-0.357 \times 4.7\times10^4\right) = \exp(-1.7 \times 10^4) \approx 10^{-7300}$$

For a two-face two-phase move the exponent is around $9.5 \times 10^4$. Either
way the number underflows to **identically zero** in double precision.

Across the resolutions in this dataset:

| voxel | $J_{NP}$ | $k_BT/J_{NP}$ | acceptance of smallest barrier |
|---|---|---|---|
| 17.9 nm | 6.41 × 10⁻¹⁶ J | 2.6 × 10⁻⁵ | 0.000 |
| 20.0 nm | 8.00 × 10⁻¹⁶ J | 2.1 × 10⁻⁵ | 0.000 |
| 29.1 nm | 1.70 × 10⁻¹⁵ J | 9.9 × 10⁻⁶ | 0.000 |

### 7.3 What that means

**At the physical temperature, Metropolis acceptance of area-increasing moves is
numerically identical to forbidding them.** The proposed rescue does not
function.

The reason is physical, not numerical. A 20 nm nickel voxel has a volume of
8 × 10⁻²⁴ m³ and, at 9.14 × 10²⁸ atoms per m³, contains about **7 × 10⁵ atoms**.
Thermal fluctuation moves individual atoms; it does not move blocks of
three-quarters of a million of them in one step. The coarse-grained lattice has
discarded exactly the length scale at which thermal activation operates, so
putting a Boltzmann factor on a whole-voxel move asks the thermostat to do
something it cannot do.

To get any acceptance at all you need $kT/J_{NP} \approx 0.078$ — about **three
thousand times** the physical temperature.

So any voxel KMC at tomographic resolution that visibly crosses barriers is
running at an **effective temperature**, not the furnace temperature. That is not
automatically illegitimate — an effective temperature is a defensible
coarse-graining parameter standing in for sub-voxel fluctuation. But it must be
named as such. It has no literature anchor, it cannot be frozen from a
measurement, and a study that reports "simulated at 950 °C" while using it is
describing something other than what it did.

**And lifting the restriction is necessary, not sufficient.** In probe runs at a
range of effective temperatures on the model dumbbell, the operator did move —
and produced surface *roughening*, not neck thinning. Low effective temperature
thickened the neck (correct curvature-driven flow, material migrating to the
concave region); high effective temperature roughened everything. No setting
produced Rayleigh break-up. Permitting uphill moves is a precondition for the
mechanism, not a route to it.

---

## 8. What this means

### 8.1 The narrow conclusion

For single-voxel-swap operators on tomographic data:

- Validating with monotone surface-area reduction **structurally excludes**
  Rayleigh-type neck break-up. Analytically for a strict criterion; empirically,
  over ~2.6 × 10⁵ moves, for the non-strict one.
- Such operators **manufacture TPB** — 3.7–5.0× while passing the gate, still
  1.33–1.53× when the gate is tightened to strict inequality. TPB conclusions
  from this class of operator are not trustworthy.
- **Percolation is not monotone** under volume-conserving rearrangement, so
  metrics must be chosen accordingly and R_Ni's monotonicity claim is specific to
  removal.
- **Finite temperature does not rescue any of this** at physical values, because
  the voxel-face energy exceeds $k_BT$ by four to five orders of magnitude.

Be careful about the scope. This is not "voxel methods cannot work." Zhou et al.
(2023) built a working Potts KMC for NiO–YSZ *sintering* and validated it against
FIB-SEM data. The claim is about this operator class *combined with this validity
criterion*, applied to *degradation*.

### 8.2 The mechanism question is open

Nothing here explains the anticorrelation of §2.

Dewetting and retraction of nickel into compact particles would break long-range
connectivity while comparatively preserving the Ni/YSZ contact perimeter where
TPB lives, which is the right shape of explanation. It is an **untested
possibility, not a surviving candidate** — nothing in this work eliminates its
competitors, and no run either supports or excludes it. Saying "the remaining
mechanism" would be inheriting an elimination argument that was never made.

What can be said: no operator respecting monotone area reduction thinned a neck
at any intensity, seed, move budget or acceptance rule tried, strictly
area-lowering or finite-temperature. So the anticorrelation stands as an open
puzzle that this class of method is not equipped to pose, let alone answer.

If you want to pursue it, phase-field is the established route — it evolves a
continuous field by deterministic PDE, needs no thermal activation to cross a
barrier, and does not have a voxel-face energy scale to fall foul of. Note it is
not a guaranteed answer either: Hoffrogge et al. (2023) ran multiphase-field on
FIB-SEM reconstructions at several wetting angles and found TPB evolution
comparatively *insensitive* to the angle, which is not what a naive dewetting
story predicts.

### 8.3 The general lesson

Every convention in this project was individually reasonable.

Validate a coarsening operator by requiring surface area to fall — correct
thermodynamics. Delete disconnected nickel — it is electrochemically dead. Build
synthetic test structures fully connected — a clean initial condition. Rank sites
by curvature — that is what drives surface diffusion. Use a Metropolis rule at
the experimental temperature — standard statistical mechanics.

Each defensible in isolation. In combination:

- the area gate excludes the mechanism it is validating;
- the pruning step rewrites the metric doing the validating;
- the synthetic structure hides both;
- the curvature ranking does not deliver the quantity being gated on;
- and the temperature that would fix the first problem is four orders of
  magnitude too small to do anything.

**No individual decision was wrong. The combination was.** That is the failure
mode this project actually documents, and it is not specific to Ni–YSZ or to
voxels. It is what happens whenever a validity criterion is chosen for its
correctness in isolation rather than tested against the specific mechanism it is
supposed to certify.

The practical version, for your own work: **for every validity criterion you
impose, ask what it does *not* constrain, and check whether the physics you care
about lives there.**

---

## Appendix A: How this was found

*This appendix is deliberately separate from the argument above. It is included
because omitting it would misrepresent how the numbers in this document came to
be, and because it demonstrates the manuscript's own thesis inside the project's
own code.*

The area barrier was originally reported with a different explanation. The
greedy operator was recorded as accepting **zero moves at every intensity** — the
structure returned bit-identical — and the conclusion drawn was that a
sphere–neck–sphere body is already a local minimum under single-voxel swaps, so
the operator is frozen at the barrier.

That was wrong, and the cause was one line.

The neighbour count $n_N$ was computed by convolving the nickel mask with the
standard 6-connectivity structuring element from the image-processing library.
That element has **sum 7, not 6 — it includes the centre voxel.** Convolving with
it returns $n_N + 1$ on a nickel site and $n_N + 0$ on a pore site. The predicate
$n_N(a) \le n_N(b)$ therefore evaluated as $n_N(a) + 1 \le n_N(b)$, i.e.
$n_N(a) < n_N(b)$ — strictly stronger than specified, silently rejecting **every
area-neutral move**.

On the original test structure that alone produced the recorded zero: the
extremal comparison was 4 > 3 under the biased count and 3 = 3, hence admissible,
under the correct one. One asymmetric off-by-one, and an operator that accepts
96 % of proposals reads as an operator that accepts none.

Three further defects surfaced from pulling that thread:

**The acceptance rate was a degenerate return value.** The function returned
`accepted / max(proposed, 1)`, which yields 0.0 both when every move is rejected
*and when no move is ever proposed*. Two opposite meanings, one number. Recording
`proposed` and `accepted` separately — one line — would have exposed the whole
chain immediately.

**Batching invalidated the identity.** $\Delta A = 2[n_N(a) - n_N(b)]$ is exact
for **one** swap against a **current** neighbour field. The implementation applied
about 59 moves per round, all ranked once against a stale field. The moves
interact and the per-move guarantees do not sum: applied in batches the operator
*raises* the area it is meant to lower, 0.45052 → 0.46365. Everything in §5 uses
strictly sequential application for this reason.

**Tie-breaking was load-bearing and unspecified.** With ~80 % of moves
area-neutral, *which* equal-cost move you take determines how much area actually
comes off. Last-in-first-out ordering achieves ΔS = −0.00024 where uniform random
selection achieves −0.00955 on the same structure at the same 1220-move budget —
a factor of 40 from a choice nobody had written down. It is now part of the
operator specification, frozen as uniform random among equal-ΔA candidates.

And one found while verifying the fix for the previous one: the adjacency retry
rejected a candidate without removing it from the pool, so under a deterministic
selection policy it drew the same voxel repeatedly and halted the run after zero
moves — the same symptom as the original bug, from an unrelated cause.

Three verification passes, four defects, each in code that had already been
checked once. Two published figures also turned out to be overstated on
re-derivation from source: a TPB collapse described as reaching "a few percent"
of pristine actually reaches 4.8–10.2 %, and a synthetic surface-area rise quoted
as 1.3–6.0 % is 1.4–5.3 % in the recorded class means.

The reason this belongs in a teaching document is not confession. It is that the
project's central claim — *individually reasonable conventions become wrong in
combination* — was demonstrated first, and most expensively, inside the project's
own implementation. Using the library's default structuring element is
reasonable. Returning a rate rather than its two components is reasonable.
Batching for speed is reasonable. Breaking ties by list order is reasonable.

The composition of those four reasonable decisions produced a published result
that was exactly backwards.

**Practical habits that follow.** Record raw counts, never only their ratios.
Check that your library's structuring element contains what you think it does.
Derive the exact quantity rather than a proxy, then verify it numerically against
brute force. Write down every tie-breaking and ordering decision as part of the
specification. And re-derive inherited numbers from source before repeating them,
including your own.

---

## Appendix B: Notation

| symbol | meaning |
|---|---|
| $n_N(x)$ | number of Ni 6-neighbours of voxel $x$, **centre excluded** |
| $n_Y(x)$, $n_P(x)$ | likewise for YSZ and pore; $n_N + n_Y + n_P = 6$ (interior) |
| $A$ | exposed nickel area, in face units |
| $S_{\text{spec}}$ | $A / V_{\text{Ni}}$, specific surface area |
| $\Delta A$ | change in $A$ from one swap; $= 2[n_N(a) - n_N(b)]$ |
| $\Delta E$ | change in bond energy; three-phase form in §7.1 |
| $J_{NP}, J_{NY}, J_{YP}$ | bond energies: Ni–pore, Ni–YSZ, YSZ–pore |
| $\theta$ | Ni-on-YSZ contact angle, ≈130° as prepared |
| $P_{\text{span}}$ | fraction of Ni in a face-to-face spanning cluster |
| $P_{\text{reach}}$ | fraction of Ni reachable from one face |
| $R_{\text{Ni}}$ | spanning-cluster Ni voxels ÷ **pristine** Ni voxels |
| TPB | triple-phase boundary length per volume, µm⁻² |
| $\Phi_{\text{Ni}}$ | nickel volume fraction |

## Appendix C: Where the numbers live

| quantity | source |
|---|---|
| measured retention, all metrics | `out/phase6/phase6_comparison_table.csv` |
| area barrier, dumbbell | `out/project2/o5v2_area_barrier.csv` |
| erosion TPB excursion | `out/project2/c1real_rni_gate.csv` |
| synthetic redistribution ratios | `out/project2/O5_REPORT.md` |
| real-ROI gate, all five conditions | `scripts/project2/o7_gate_a1v2_real.py` |
| strict-inequality control | `scripts/project2/o7_strict_inequality.py` |
| tie-breaking sensitivity | `scripts/project2/o7_tiebreak_sensitivity.py` |
| three-phase identity verification | `scripts/project2/o7_derivation_check.py` |
| athermal barrier calculation | same file, final block |
| full audit trail | `out/project2/O7_O5V2B_RERUN_REPORT.md` |

## References

Björnsson, M., Hansen, K.V., Jørgensen, P.S., Chen, M., Hauch, A., Simonsen, S.B.
(2026). Establishing a framework for determining the wetting of Ni on YSZ in SOC
fuel electrodes. *J. Power Sources*. doi:10.1016/j.jpowsour.2025.239069

Chen, H.-Y. et al. (2011). Simulation of coarsening in three-phase solid oxide
fuel cell anodes. *J. Power Sources*.

Herring, C. (1950). Effect of change of scale on sintering phenomena. *J. Appl.
Phys.* **21**, 301. — the origin of the $t \propto L^4$ scaling for
surface-diffusion-controlled processes.

Hoffrogge, P.W. et al. (2023). Performance estimation by multiphase-field
simulations and transmission-line modeling of nickel coarsening in FIB-SEM
reconstructed Ni-YSZ SOFC anodes I: influence of wetting angle. *J. Power
Sources*.

Jørgensen, P.S., Bowen, J.R. (2014). On the accuracy of triple phase boundary
lengths calculated from tomographic image data.

Lorenz, C.D., Ziff, R.M. (1998). Precise determination of the bond percolation
thresholds for the simple cubic lattice.

Monaco, F. et al. (2019). Degradation of Ni-YSZ electrodes in solid oxide cells.

Nichols, F.A. (1965). On the spheroidization of rod-shaped particles of finite
length. — the solid-state Rayleigh instability.

Pecho, O. et al. (2015). 3D microstructure effects in Ni-YSZ anodes: prediction
of effective transport properties and optimization of electrode performance
(*Materials*); and the companion TPB study.

Sarantaridis, D., Atkinson, A. (2007). Redox cycling of Ni-based solid oxide fuel
cell anodes: a review.

Shimura, T. et al. (2016). Evaluation of triple phase boundary from tomographic
data.

Simwonis, D., Tietz, F., Stöver, D. (2000). Nickel coarsening in annealed
Ni/8YSZ anode substrates for solid oxide fuel cells.

Trini, M. et al. (2021). Towards the validation of a phase field model for Ni
coarsening in solid oxide cells. *Acta Materialia*.

Yang, Y. et al. (2023). Ni coarsening in Ni-yttria stabilized zirconia
electrodes: three-dimensional quantitative phase-field simulations supported by
ex-situ ptychographic nano-tomography. *Acta Materialia*.

Zhou, Y. et al. (2023). Kinetic Monte Carlo (KMC) simulation of sintering of
nickel oxide-yttria stabilized zirconia composites: model, parameter calibration
and validation. *Materials & Design*.
