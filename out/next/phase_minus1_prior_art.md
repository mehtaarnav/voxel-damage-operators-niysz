# Phase −1 — prior-art memo

Search date 2026-08-10. 13 queries run across the 5 clusters in
`out/next/EXECUTION_SPEC.md` §5 (not all 17 listed queries were run verbatim;
several near-duplicates were merged, noted below). All hits are web-search
summaries, not full-text reads — treat as a literature MAP, not a verified
literature review. Full query log in the transcript; the 13 queries actually
run are listed under each cluster below.

## Verdict

**No direct equivalent found.** No located study (i) holds Ni volume fraction
AND mean Ni particle size fixed while (ii) independently varying the
neck-size distribution or engineering a coarse percolating backbone, and
(iii) measures retained Ni electronic percolation under a redox-like damage
model, while (iv) tracking the TPB tradeoff. Gate G−1 exit criterion met
without narrowing the claim. Several adjacent threads exist and must be cited;
they are listed below with an explicit statement of how this study differs
from each.

## Closest prior work, and how we differ

### 1. Holzer/Pecho/Iwanschitz group (the group whose own data we reuse)
The transport paper (ma8095265, already our ground truth) states the loss
mechanism qualitatively: Ni percolation loss in the FINE anode is
coarsening-driven (particle growth), while YSZ percolation loss in the COARSE
anode is attributed to "weak bottlenecks... due to lower sintering activity of
the coarse YSZ" — i.e., a neck-strength argument, in their own words, for a
DIFFERENT phase (YSZ, not Ni) and a different failure mode (YSZ disintegration
under Ni-agglomeration-induced stress, not Ni-network attrition). They also
*propose* — but do not test — mixing fine and coarse NiO/YSZ powder fractions
across the anode's functional and current-collector layers as a redox-stability
strategy [ma8095265, Discussion]. This is the closest real-world design
intuition to our hypothesis.
**How we differ:** their neck-strength argument is about YSZ, ours is about Ni;
their "mixing" proposal varies particle size distribution (not neck width at
fixed particle size), and is never tested against a controlled comparison; and
all six of their real samples co-vary particle size and neck size, which is
exactly the confound our synthetic study exists to break.

### 2. TPB pathway / bottleneck-radius analysis (Yan et al., and related)
"Triple phase boundary specific pathway analysis for quantitative
characterization of solid oxide cell electrode microstructure"
(ScienceDirect, S037877531500066X) introduces, for each TPB site, the pathway
length AND bottleneck radius to reach it — a per-site neck metric close in
spirit to our neck-width-quantile approach — applied to a Ni/ScYSZ
*reduction protocol* study.
**How we differ:** characterizes existing microstructures; does not engineer
neck geometry, does not fix particle size, not a redox-degradation study.

### 3. Random-packing / particle-size-distribution electrode models
"Random-packing model for solid oxide fuel cell electrodes with particle size
distributions" (ScienceDirect, S0378775310017040) models coordination number,
percolation probability, and TPB from a sphere-packing generator parametrized
by particle size distribution and composition — methodologically the closest
template to a Family A/B sphere-packing generator.
**How we differ:** optimizes *pristine* electrode design (no damage model);
varies particle-size distribution broadly (including bimodal ratios), not
neck width at *fixed* particle size; no redox-cycling or aging step.

### 4. Stochastic 3D generative + physics-based aging models (Ulm group:
Neumann, Schmidt, and collaborators)
"A stochastic geometrical 3D model for time evolution simulation of
microstructures in SOC-electrodes" and "A time-continuous approach to
analyzing anode aging in solid-oxide fuel cells via stochastic 3D
microstructure modeling and physics-based simulations" (2026) combine a
Gaussian-random-field / excursion-set generative model with a physics-based
aging step, predicting coarsening and reporting a shift toward anisotropy
after 3800 h operation. This is the closest single combination of
"synthetic generator + degradation model" for Ni-YSZ-like electrodes found in
this search.
**How we differ:** their aging mechanism is thermal Ni coarsening over
long-term operation (particle growth), which is the SAME confound (particle
size and connectivity co-evolve) our study is designed to break, not a
redox-cycle-like damage operator; they do not hold particle size fixed to
isolate a neck effect; "constant particle size, varying neck distribution" is
not their design variable.

### 5. Powder-metallurgy / general porous-media neck-conductivity models
Discrete-element and resistor-network models of sintered porous-material
conductivity DO parametrize neck (inter-particle contact) size independently
of particle size (e.g., "Discrete element model for effective electrical
conductivity of spark plasma sintered porous materials," ScienceDirect
S2196438625005091), and note the neck/contact-area effect on conductivity
directly.
**How we differ:** not Ni-YSZ, not a three-phase (Ni/YSZ/pore) system with TPB,
not framed around redox damage, and typically an abstracted resistor-network
rather than a full segmented 3D microstructure with watershed particle/TPB
recomputation.

### 6. Constrictivity as a named quantity (same group as #1)
The transport paper's own M-factor model already names a constriction factor
β = (r_min/r_max)², i.e., neck size relative to particle "bulge" size, as one
of three transport-controlling parameters (with Φ_eff and τ) — so the
FIELD already has vocabulary for "neck relative to particle size" as a
transport predictor. It has never been used as an independently engineered,
controlled variable for a redox-retention comparison.

## Queries run

**Cluster 1 (neck/decoupling):** "Ni-YSZ anode neck size distribution redox
tolerance microstructure percolation"; "decouple neck size particle size
percolation porous electrode simulation sintering".

**Cluster 2 (bimodal/backbone):** "bimodal Ni particle size distribution SOFC
anode redox cycling percolation"; ""coarse backbone" fine decoration Ni-YSZ
infiltrated anode TPB percolation redox stability".

**Cluster 3 (TPB vs percolation tradeoff):** "triple phase boundary
percolation tradeoff Ni-YSZ redox cycling quantitative tomography";
""triple phase boundary" pathway bottleneck radius analysis SOFC electrode
quantitative characterization Yan".

**Cluster 4 (synthetic generation + redox damage):** "synthetic microstructure
generation Ni-YSZ anode percolation TPB sphere packing simulation";
"stochastic geometry model SOFC electrode redox degradation simulation
microstructure generator"; "time-continuous stochastic 3D microstructure
modeling anode aging solid oxide fuel cell physics-based simulation";
"random packing model solid oxide fuel cell electrode particle size
distribution percolation threshold neck".

**Cluster 5 (named groups):** "Shikazono Kishimoto synthetic microstructure
Ni-YSZ percolation TPB reconstruction"; "Cronin OR Wilson OR Barnett Ni-YSZ
redox cycling percolation network model X-ray tomography"; "constrictivity
bottleneck engineering redox stability optimization Ni-YSZ microstructure
design".

Not run (judged as unlikely to surface material beyond the above, given
overlapping results already returned): the "Neumann OR Schmidt" variant
(subsumed by cluster 4's stochastic-model queries, same group); "sintering
neck radius independent particle size percolation SOFC electrode simulation"
(subsumed by cluster 1's decoupling query, same near-duplicate).

## Caveat

This is a web-search literature map, not a systematic review: titles and
abstracts only, no full-text verification, and search-engine summaries were
used directly for some claims above (flagged where a source is named). Before
any Path-A (positive-result) writeup, the 4 papers named above as closest
prior work should be read in full and cited properly; this memo is sufficient
to clear gate G−1 (claim is not duplicated) but not sufficient as a
publication-ready literature review.
