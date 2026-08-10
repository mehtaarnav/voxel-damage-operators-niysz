# Phase 1 — ground-truth reference table

Every number below is transcribed from the two Pecho et al. papers. Source key:

- **`T-S4`** — Pecho et al., Materials 2015, 8(9), 5554-5585 (doi:10.3390/ma8095265), SUPPLEMENTARY Table S4, p. S2
- **`T-S2`** — same paper, SUPPLEMENTARY Table S2 (image dimensions)
- **`T-S1`** — same paper, SUPPLEMENTARY Table S1 (raw powder PSD)
- **`P-F7`** — Pecho et al., Materials 2015, 8(10), 7129-7147 (doi:10.3390/ma8105370), FIGURE 7 A/B - DIGITIZED FROM A BAR CHART (no table exists)

> **Provenance warning.** Neither paper's *main text* contains per-sample numeric tables. `T-S4` comes from the transport paper's supplementary PDF and is exact. `P-F7` is **digitized from a bar chart** and is good to about ±0.05 µm⁻²; its self-consistency check is below.

## Sample identity and acquisition

| sample      | zenodo_folder      | grain   | state    |   nx |   ny |   nz |   vx_nm |   vy_nm |   vz_nm | zenodo_label_note              |
|:------------|:-------------------|:--------|:---------|-----:|-----:|-----:|--------:|--------:|--------:|:-------------------------------|
| fine_pre    | 3_Rx36_Segmented   | fine    | pristine |  995 | 1304 |  733 |   19.53 |   19.53 |   20    | Ni white, YSZ gray, pore black |
| medium_pre  | 4_Rx37_Segmented   | medium  | pristine |  960 | 1110 |  610 |   24.41 |   24.41 |   25    | Ni white, YSZ gray, pore black |
| coarse_pre  | 5_Rx38_Segmented   | coarse  | pristine |  744 | 1417 |  456 |   29.14 |   29.14 |   30    | Ni gray, YSZ white, pore black |
| fine_post   | 6_Rx41-1_Segmented | fine    | degraded | 1171 | 1343 |  461 |   19.53 |   19.53 |   20.47 | Ni white, YSZ gray, pore black |
| medium_post | 7_Rx41-2_Segmented | medium  | degraded | 1318 | 1520 |  459 |   17.9  |   17.9  |   25    | Ni white, YSZ gray, pore black |
| coarse_post | 8_Rx41-3_Segmented | coarse  | degraded | 1368 | 1630 |  500 |   17.9  |   17.9  |   25    | Ni white, YSZ gray, pore black |

*(dimensions independently confirmed by `T-S2` and by the Zenodo `2_3D_Data_Info.xlsx`; the two agree exactly)*

## Ni — transport-relevant parameters (`T-S4`)

| sample      | grain   | state    |   Phi |     P |   Phi_eff |   beta |   tau |   M_pred |
|:------------|:--------|:---------|------:|------:|----------:|-------:|------:|---------:|
| fine_pre    | fine    | pristine | 0.322 | 0.985 |     0.317 |  0.275 | 1.219 |    0.071 |
| medium_pre  | medium  | pristine | 0.25  | 0.965 |     0.241 |  0.26  | 1.341 |    0.033 |
| coarse_pre  | coarse  | pristine | 0.229 | 0.959 |     0.22  |  0.22  | 1.605 |    0.011 |
| fine_post   | fine    | degraded | 0.222 | 0.809 |     0.179 |  0.188 | 1.375 |    0.019 |
| medium_post | medium  | degraded | 0.233 | 0.884 |     0.206 |  0.345 | 1.358 |    0.029 |
| coarse_post | coarse  | degraded | 0.244 | 0.886 |     0.216 |  0.372 | 1.673 |    0.011 |

## YSZ — transport-relevant parameters (`T-S4`)

| sample      | grain   | state    |   Phi |     P |   Phi_eff |   beta |   tau |   M_pred |
|:------------|:--------|:---------|------:|------:|----------:|-------:|------:|---------:|
| fine_pre    | fine    | pristine | 0.421 | 0.999 |     0.421 | 0.367  | 1.108 |    0.173 |
| medium_pre  | medium  | pristine | 0.388 | 0.986 |     0.383 | 0.095  | 1.176 |    0.071 |
| coarse_pre  | coarse  | pristine | 0.384 | 0.923 |     0.354 | 0.007  | 1.889 |    0.002 |
| fine_post   | fine    | degraded | 0.312 | 0.961 |     0.3   | 0.088  | 1.43  |    0.02  |
| medium_post | medium  | degraded | 0.376 | 0.869 |     0.327 | 0.042  | 1.353 |    0.022 |
| coarse_post | coarse  | degraded | 0.324 | 0.184 |     0.06  | 0.0001 | 1.1   |    0.001 |

## Pore — transport-relevant parameters (`T-S4`)

| sample      | grain   | state    |   Phi |     P |   Phi_eff |   beta |   tau |   M_pred |
|:------------|:--------|:---------|------:|------:|----------:|-------:|------:|---------:|
| fine_pre    | fine    | pristine | 0.254 | 0.988 |     0.251 |  0.271 | 1.324 |    0.037 |
| medium_pre  | medium  | pristine | 0.362 | 0.998 |     0.361 |  0.55  | 1.11  |    0.17  |
| coarse_pre  | coarse  | pristine | 0.387 | 0.999 |     0.386 |  0.563 | 1.103 |    0.19  |
| fine_post   | fine    | degraded | 0.466 | 0.999 |     0.466 |  0.547 | 1.073 |    0.26  |
| medium_post | medium  | degraded | 0.39  | 0.998 |     0.389 |  0.594 | 1.081 |    0.216 |
| coarse_post | coarse  | degraded | 0.432 | 0.992 |     0.428 |  0.487 | 1.082 |    0.22  |

## TPB densities (`P-F7`, digitized)

| sample      | grain   | state    |   TPB_total (um^-2) |   TPB_active (um^-2) |
|:------------|:--------|:---------|--------------------:|---------------------:|
| fine_pre    | fine    | pristine |                2.65 |                 2.38 |
| medium_pre  | medium  | pristine |                2.03 |                 1.79 |
| coarse_pre  | coarse  | pristine |                1.07 |                 0.78 |
| fine_post   | fine    | degraded |                1.97 |                 1.18 |
| medium_post | medium  | degraded |                1.19 |                 0.55 |
| coarse_post | coarse  | degraded |                0.65 |                 0.05 |

### Digitization self-consistency

Figure 7 panel C prints the relative change independently of panels A and B. Recomputing it from the digitized A/B values reproduces panel C to within 0.7 percentage points:

| kind   | grain   |   pre |   post |   delta_from_AB_pct |   delta_printed_panelC_pct |   abs_diff_pp |
|:-------|:--------|------:|-------:|--------------------:|---------------------------:|--------------:|
| total  | fine    |  2.65 |   1.97 |               -25.7 |                        -26 |           0.3 |
| total  | medium  |  2.03 |   1.19 |               -41.4 |                        -42 |           0.6 |
| total  | coarse  |  1.07 |   0.65 |               -39.3 |                        -39 |           0.3 |
| active | fine    |  2.38 |   1.18 |               -50.4 |                        -50 |           0.4 |
| active | medium  |  1.79 |   0.55 |               -69.3 |                        -70 |           0.7 |
| active | coarse  |  0.78 |   0.05 |               -93.6 |                        -93 |           0.6 |

**Digitization check: PASS**

## Raw powder particle size (`T-S1`)

| powder | d50 (um) |
|---|---|
| NiO | 0.62 |
| YSZ_fine | 0.47 |
| YSZ_medium | 3.33 |
| YSZ_coarse | 10.19 |

## Definitions used by the authors (these differ from the obvious choices)

- **`P` (percolation factor)** — *"describes the fraction of a phase, which forms a connected network, and it can be obtained from the MIP-PSD analysis"* (ma8095265, Sec. 3). This is a simulated-intrusion measure of the fraction **reachable from a boundary**, not a connected-component face-to-face spanning fraction. It is therefore an upper bound on a two-face spanning fraction.
- **`TPB total`** — all three-phase lines in the full cube, no connectivity check. Length obtained *"based on the skeletonization of TPB-voxels in each object"* (ma8105370, Methods) — **not** voxel-edge counting.
- **`TPB active`** — TPB lines where all three phases pass a connectivity check toward their relevant border, counted in a central sub-cube to suppress boundary truncation.
- **`beta`** — constriction factor (r_min/r_max)^2 from MIP-PSD vs c-PSD.
- **`tau`** — mean geodesic tortuosity.
