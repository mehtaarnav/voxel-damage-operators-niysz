"""
PHASE 1 ground-truth table, transcribed verbatim from the two Pecho et al. papers.

EVERY number below carries its source.  Nothing here is inferred or interpolated
except where explicitly marked `digitized`.

SOURCES
-------
[T-S4] Pecho, Stenzel, Iwanschitz, Gasser, Neumann, Schmidt, Prestat, Hocker,
       Flatt, Holzer, "3D Microstructure Effects in Ni-YSZ Anodes: Prediction of
       Effective Transport Properties and Optimization of Redox Stability",
       Materials 2015, 8(9), 5554-5585, doi:10.3390/ma8095265.
       SUPPLEMENTARY Table S4 ("Summary of transport-relevant microstructure
       parameters and the predicted M-factors (Mpred) obtained using Equation
       (4)").  File materials-08-05265-s001.pdf, page S2.
       -> exact printed values, 3 decimal places.

[T-S2] same paper, supplementary Table S2 ("Summary of image data for 3D
       quantitative microstructure analyses").  Voxel sizes / matrices.
       -> used only as a cross-check against the Zenodo 2_3D_Data_Info.xlsx.

[T-S1] same paper, supplementary Table S1: raw-powder laser-diffraction PSD.

[T-txt] same paper, main text Section 3, which independently restates several
       Table S4 numbers (porosity 25/36/39 vol%; Ni P 0.99->0.80 fine, 0.90
       coarse; YSZ P -> 0.18 coarse).  All agree with [T-S4].

[P-F7] Pecho, Mai, Muench, Hocker, Flatt, Holzer, "3D Microstructure Effects in
       Ni-YSZ Anodes: Influence of TPB Lengths on the Electrochemical
       Performance", Materials 2015, 8(10), 7129-7147, doi:10.3390/ma8105370.
       FIGURE 7, panels A (before), B (after), C (relative change).
       -> DIGITIZED FROM A BAR CHART.  No table of these values exists in the
       paper or its supplement.  Read to the nearest 0.01-0.05 um^-2.
       Self-consistency check: recomputing panel C from panels A and B
       reproduces the printed percentages to within 1 percentage point
       (see validate_tpb_selfconsistency() below), so the digitization is sound
       at roughly the +/-2% level.  Treat as +/-0.05 um^-2.

DEFINITIONS AS USED BY THE AUTHORS  (these differ from the obvious choices and
must be honoured when comparing)
---------------------------------
* Phi   : volume fraction of the phase in the analysed window.
* P     : "percolation factor ... describes the fraction of a phase, which forms
          a connected network, and it can be obtained from the MIP-PSD analysis"
          [ma8095265, Section 3].  This is a simulated mercury-intrusion
          (morphological opening + reachability) measure, NOT a
          connected-component face-to-face spanning fraction.  It measures the
          fraction of the phase REACHABLE FROM A BOUNDARY, so it is an upper
          bound on a two-face spanning fraction.  Our pipeline reports both.
* Phi_eff = Phi * P.
* beta  : constriction factor (r_min/r_max)^2 from MIP-PSD vs c-PSD.
* tau   : mean geodesic tortuosity.
* TPB total  : all three-phase boundary lines in the full cube, no connectivity
          check.  Length of each TPB line obtained "based on the skeletonization
          of TPB-voxels in each object" [ma8105370, Methods].  NOT voxel-edge
          counting.
* TPB active : TPB lines where the connectivity check is positive for each of
          the three phases toward its relevant border (e.g. Ni connected to the
          current-collector side), evaluated in a central sub-cube to suppress
          boundary-truncation effects.
"""

from __future__ import annotations

import pandas as pd

# --------------------------------------------------------------------------
# Sample identity.  Zenodo folder <-> paper row.
# Dimensions/voxel sizes from the Zenodo 2_3D_Data_Info.xlsx AND independently
# from [T-S2]; the two agree exactly.
# --------------------------------------------------------------------------
SAMPLES = [
    # key,        zenodo folder,        grain,    state,      nx,   ny,   nz,  vx,    vy,    vz
    ("fine_pre",   "3_Rx36_Segmented",   "fine",   "pristine",  995, 1304, 733, 19.53, 19.53, 20.00),
    ("medium_pre", "4_Rx37_Segmented",   "medium", "pristine",  960, 1110, 610, 24.41, 24.41, 25.00),
    ("coarse_pre", "5_Rx38_Segmented",   "coarse", "pristine",  744, 1417, 456, 29.14, 29.14, 30.00),
    ("fine_post",  "6_Rx41-1_Segmented", "fine",   "degraded", 1171, 1343, 461, 19.53, 19.53, 20.47),
    ("medium_post","7_Rx41-2_Segmented", "medium", "degraded", 1318, 1520, 459, 17.90, 17.90, 25.00),
    ("coarse_post","8_Rx41-3_Segmented", "coarse", "degraded", 1368, 1630, 500, 17.90, 17.90, 25.00),
]

# Phase-label convention stated in the Zenodo metadata, per sample (column AE).
# NOTE THE ODD ONE OUT: 5_Rx38 is documented with Ni and YSZ SWAPPED relative to
# every other stack.  This is verified empirically in Phase 2, not trusted.
ZENODO_LABEL_NOTE = {
    "fine_pre":    "Ni white, YSZ gray, pore black",
    "medium_pre":  "Ni white, YSZ gray, pore black",
    "coarse_pre":  "Ni gray, YSZ white, pore black",   # <-- INVERTED
    "fine_post":   "Ni white, YSZ gray, pore black",
    "medium_post": "Ni white, YSZ gray, pore black",
    "coarse_post": "Ni white, YSZ gray, pore black",
}

# --------------------------------------------------------------------------
# [T-S4] verbatim.  Order of columns in the paper:
#   Before Redox Cycling: Fine, Medium, Coarse | After Redox Cycling: Fine, Medium, Coarse
# --------------------------------------------------------------------------
_ORDER = ["fine_pre", "medium_pre", "coarse_pre",
          "fine_post", "medium_post", "coarse_post"]

_TS4 = {
    "Ni": {
        "Phi":     [0.322, 0.250, 0.229, 0.222, 0.233, 0.244],
        "P":       [0.985, 0.965, 0.959, 0.809, 0.884, 0.886],
        "Phi_eff": [0.317, 0.241, 0.220, 0.179, 0.206, 0.216],
        "beta":    [0.275, 0.260, 0.220, 0.188, 0.345, 0.372],
        "tau":     [1.219, 1.341, 1.605, 1.375, 1.358, 1.673],
        "M_pred":  [0.071, 0.033, 0.011, 0.019, 0.029, 0.011],
    },
    "YSZ": {
        "Phi":     [0.421, 0.388, 0.384, 0.312, 0.376, 0.324],
        "P":       [0.999, 0.986, 0.923, 0.961, 0.869, 0.184],
        "Phi_eff": [0.421, 0.383, 0.354, 0.300, 0.327, 0.060],
        "beta":    [0.367, 0.095, 0.007, 0.088, 0.042, 0.0001],
        "tau":     [1.108, 1.176, 1.889, 1.430, 1.353, 1.100],
        "M_pred":  [0.173, 0.071, 0.002, 0.020, 0.022, 0.001],
    },
    "Pore": {
        "Phi":     [0.254, 0.362, 0.387, 0.466, 0.390, 0.432],
        "P":       [0.988, 0.998, 0.999, 0.999, 0.998, 0.992],
        "Phi_eff": [0.251, 0.361, 0.386, 0.466, 0.389, 0.428],
        "beta":    [0.271, 0.550, 0.563, 0.547, 0.594, 0.487],
        "tau":     [1.324, 1.110, 1.103, 1.073, 1.081, 1.082],
        "M_pred":  [0.037, 0.170, 0.190, 0.260, 0.216, 0.220],
    },
}

# --------------------------------------------------------------------------
# [P-F7] digitized bar chart, um of TPB per um^3 == um^-2.
# --------------------------------------------------------------------------
_TPB = {
    "TPB_total":  {"fine_pre": 2.65, "medium_pre": 2.03, "coarse_pre": 1.07,
                   "fine_post": 1.97, "medium_post": 1.19, "coarse_post": 0.65},
    "TPB_active": {"fine_pre": 2.38, "medium_pre": 1.79, "coarse_pre": 0.78,
                   "fine_post": 1.18, "medium_post": 0.55, "coarse_post": 0.05},
}
# Panel C, printed relative changes (%), used only to validate the digitization.
_TPB_DELTA_PRINTED = {"total": {"fine": -26, "medium": -42, "coarse": -39},
                      "active": {"fine": -50, "medium": -70, "coarse": -93}}

# [T-S1] raw powder PSD, um
RAW_POWDER_D50 = {"NiO": 0.62, "YSZ_fine": 0.47, "YSZ_medium": 3.33,
                  "YSZ_coarse": 10.19}


def validate_tpb_selfconsistency(tol_pp=1.5):
    """Recompute Figure 7C from Figure 7A/B; confirms the digitization."""
    rows = []
    ok = True
    for kind in ("total", "active"):
        for grain in ("fine", "medium", "coarse"):
            pre = _TPB[f"TPB_{kind}"][f"{grain}_pre"]
            post = _TPB[f"TPB_{kind}"][f"{grain}_post"]
            calc = 100.0 * (post / pre - 1.0)
            printed = _TPB_DELTA_PRINTED[kind][grain]
            d = abs(calc - printed)
            ok &= d <= tol_pp
            rows.append(dict(kind=kind, grain=grain, pre=pre, post=post,
                             delta_from_AB_pct=round(calc, 1),
                             delta_printed_panelC_pct=printed,
                             abs_diff_pp=round(d, 1)))
    return pd.DataFrame(rows), ok


def ground_truth_frame() -> pd.DataFrame:
    """One row per sample; every published quantity we could extract."""
    rows = []
    for i, key in enumerate(_ORDER):
        meta = [s for s in SAMPLES if s[0] == key][0]
        r = {
            "sample": key,
            "zenodo_folder": meta[1],
            "grain": meta[2],
            "state": meta[3],
            "nx": meta[4], "ny": meta[5], "nz": meta[6],
            "vx_nm": meta[7], "vy_nm": meta[8], "vz_nm": meta[9],
            "zenodo_label_note": ZENODO_LABEL_NOTE[key],
        }
        for phase, d in _TS4.items():
            for q, vals in d.items():
                r[f"{phase}_{q}__T-S4"] = vals[i]
        r["TPB_total_um-2__P-F7_digitized"] = _TPB["TPB_total"][key]
        r["TPB_active_um-2__P-F7_digitized"] = _TPB["TPB_active"][key]
        rows.append(r)
    df = pd.DataFrame(rows)
    # derived sanity column: phase fractions must sum to 1
    df["phi_sum_check"] = (df["Ni_Phi__T-S4"] + df["YSZ_Phi__T-S4"]
                           + df["Pore_Phi__T-S4"])
    return df


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 100)
    df = ground_truth_frame()
    print(df.T.to_string())
    print("\nPhi sum check (should all be ~1.000):")
    print(df[["sample", "phi_sum_check"]].to_string(index=False))
    chk, ok = validate_tpb_selfconsistency()
    print("\nFigure-7 digitization self-consistency:")
    print(chk.to_string(index=False))
    print("PASS" if ok else "FAIL")
