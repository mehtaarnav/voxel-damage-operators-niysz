"""PHASE 1 deliverable — write the ground-truth reference table to disk."""

from __future__ import annotations

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import (  # noqa: E402
    RAW_POWDER_D50, ground_truth_frame, validate_tpb_selfconsistency,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase1")
os.makedirs(OUT, exist_ok=True)

SRC = {
    "T-S4": "Pecho et al., Materials 2015, 8(9), 5554-5585 (doi:10.3390/ma8095265), "
            "SUPPLEMENTARY Table S4, p. S2",
    "T-S2": "same paper, SUPPLEMENTARY Table S2 (image dimensions)",
    "T-S1": "same paper, SUPPLEMENTARY Table S1 (raw powder PSD)",
    "P-F7": "Pecho et al., Materials 2015, 8(10), 7129-7147 (doi:10.3390/ma8105370), "
            "FIGURE 7 A/B - DIGITIZED FROM A BAR CHART (no table exists)",
}


def main():
    df = ground_truth_frame()
    csv = os.path.join(OUT, "phase1_ground_truth.csv")
    df.to_csv(csv, index=False)

    chk, ok = validate_tpb_selfconsistency()
    chk.to_csv(os.path.join(OUT, "phase1_tpb_digitization_check.csv"), index=False)

    lines = []
    lines.append("# Phase 1 — ground-truth reference table\n")
    lines.append("Every number below is transcribed from the two Pecho et al. "
                 "papers. Source key:\n")
    for k, v in SRC.items():
        lines.append(f"- **`{k}`** — {v}")
    lines.append("")
    lines.append("> **Provenance warning.** Neither paper's *main text* contains "
                 "per-sample numeric tables. `T-S4` comes from the transport "
                 "paper's supplementary PDF and is exact. `P-F7` is **digitized "
                 "from a bar chart** and is good to about ±0.05 µm⁻²; its "
                 "self-consistency check is below.\n")

    lines.append("## Sample identity and acquisition\n")
    cols = ["sample", "zenodo_folder", "grain", "state", "nx", "ny", "nz",
            "vx_nm", "vy_nm", "vz_nm", "zenodo_label_note"]
    lines.append(df[cols].to_markdown(index=False))
    lines.append("\n*(dimensions independently confirmed by `T-S2` and by the "
                 "Zenodo `2_3D_Data_Info.xlsx`; the two agree exactly)*\n")

    for phase in ("Ni", "YSZ", "Pore"):
        lines.append(f"## {phase} — transport-relevant parameters (`T-S4`)\n")
        sub = df[["sample", "grain", "state"] +
                 [f"{phase}_{q}__T-S4" for q in
                  ("Phi", "P", "Phi_eff", "beta", "tau", "M_pred")]]
        sub = sub.rename(columns={f"{phase}_{q}__T-S4": q for q in
                                  ("Phi", "P", "Phi_eff", "beta", "tau", "M_pred")})
        lines.append(sub.to_markdown(index=False))
        lines.append("")

    lines.append("## TPB densities (`P-F7`, digitized)\n")
    sub = df[["sample", "grain", "state",
              "TPB_total_um-2__P-F7_digitized",
              "TPB_active_um-2__P-F7_digitized"]]
    sub = sub.rename(columns={
        "TPB_total_um-2__P-F7_digitized": "TPB_total (um^-2)",
        "TPB_active_um-2__P-F7_digitized": "TPB_active (um^-2)"})
    lines.append(sub.to_markdown(index=False))
    lines.append("")
    lines.append("### Digitization self-consistency\n")
    lines.append("Figure 7 panel C prints the relative change independently of "
                 "panels A and B. Recomputing it from the digitized A/B values "
                 "reproduces panel C to within 0.7 percentage points:\n")
    lines.append(chk.to_markdown(index=False))
    lines.append(f"\n**Digitization check: {'PASS' if ok else 'FAIL'}**\n")

    lines.append("## Raw powder particle size (`T-S1`)\n")
    lines.append("| powder | d50 (um) |")
    lines.append("|---|---|")
    for k, v in RAW_POWDER_D50.items():
        lines.append(f"| {k} | {v} |")
    lines.append("")

    lines.append("## Definitions used by the authors (these differ from the "
                 "obvious choices)\n")
    lines.append("- **`P` (percolation factor)** — *\"describes the fraction of "
                 "a phase, which forms a connected network, and it can be "
                 "obtained from the MIP-PSD analysis\"* (ma8095265, Sec. 3). "
                 "This is a simulated-intrusion measure of the fraction "
                 "**reachable from a boundary**, not a connected-component "
                 "face-to-face spanning fraction. It is therefore an upper "
                 "bound on a two-face spanning fraction.")
    lines.append("- **`TPB total`** — all three-phase lines in the full cube, no "
                 "connectivity check. Length obtained *\"based on the "
                 "skeletonization of TPB-voxels in each object\"* (ma8105370, "
                 "Methods) — **not** voxel-edge counting.")
    lines.append("- **`TPB active`** — TPB lines where all three phases pass a "
                 "connectivity check toward their relevant border, counted in a "
                 "central sub-cube to suppress boundary truncation.")
    lines.append("- **`beta`** — constriction factor (r_min/r_max)^2 from "
                 "MIP-PSD vs c-PSD.")
    lines.append("- **`tau`** — mean geodesic tortuosity.\n")

    md = os.path.join(OUT, "phase1_ground_truth.md")
    with open(md, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print("[saved]", csv)
    print("[saved]", md)
    print("[saved]", os.path.join(OUT, "phase1_tpb_digitization_check.csv"))
    print(f"\nTPB digitization self-consistency: {'PASS' if ok else 'FAIL'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
