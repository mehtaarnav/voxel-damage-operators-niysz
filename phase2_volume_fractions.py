"""
PHASE 2 — load the real data, assign phases, and gate on volume fraction.

GATE (user-specified): voxel-counted volume fractions must match the published
values to within ~10-15 % relative.  If not, STOP.

The label->phase assignment comes from cmlib.phases, which uses only the
dataset's own metadata (brightness ordering), never the published numbers, so
this comparison is a real test.

Also writes, for visual inspection:
  phase2_slices_<sample>.png   a mid-stack slice, phase-coloured, plus the
                               raw label image and a per-slice volume-fraction
                               profile through the stack.
"""

from __future__ import annotations

import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE, ground_truth_frame  # noqa: E402
from cmlib.io import label_histogram, slice_paths  # noqa: E402
from cmlib.phases import assign_labels, volume_fractions  # noqa: E402

import tifffile  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase2")
os.makedirs(OUT, exist_ok=True)

GATE_REL = 0.15   # 15 % relative, the loose end of the user's 10-15 % band


def per_slice_profile(folder, mapping, step=10):
    """Volume fraction of each phase as a function of slice index."""
    ps = slice_paths(folder)[::step]
    rows = []
    for i, p in enumerate(ps):
        a = tifffile.imread(p)
        c = np.bincount(a.ravel(), minlength=256)
        t = c.sum()
        rows.append({"z": i * step,
                     "Ni": c[mapping["Ni"]] / t,
                     "YSZ": c[mapping["YSZ"]] / t,
                     "pore": c[mapping["pore"]] / t})
    return pd.DataFrame(rows)


def save_slice_figure(key, folder, mapping, gt, mine, prof):
    ps = slice_paths(folder)
    mid = tifffile.imread(ps[len(ps) // 2])

    # phase-coloured rendering: 0 pore, 1 YSZ, 2 Ni, 3 unassigned
    disp = np.full(mid.shape, 3, dtype=np.uint8)
    disp[mid == mapping["pore"]] = 0
    disp[mid == mapping["YSZ"]] = 1
    disp[mid == mapping["Ni"]] = 2
    cmap = ListedColormap(["#101010", "#8a8a8a", "#f2f2f2", "#d62728"])

    fig = plt.figure(figsize=(15.5, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.05])

    ax0 = fig.add_subplot(gs[0])
    im0 = ax0.imshow(mid, cmap="viridis", interpolation="nearest")
    ax0.set_title(f"raw labels (values {sorted(np.unique(mid).tolist())})",
                  fontsize=9)
    ax0.set_xticks([]); ax0.set_yticks([])
    fig.colorbar(im0, ax=ax0, fraction=0.035)

    ax1 = fig.add_subplot(gs[1])
    ax1.imshow(disp, cmap=cmap, interpolation="nearest", vmin=0, vmax=3)
    ax1.set_title("assigned: black=pore  gray=YSZ  white=Ni  red=unassigned",
                  fontsize=9)
    ax1.set_xticks([]); ax1.set_yticks([])

    ax2 = fig.add_subplot(gs[2])
    for ph, col in (("Ni", "C3"), ("YSZ", "C0"), ("pore", "C2")):
        ax2.plot(prof["z"], prof[ph], lw=1.0, color=col, label=f"{ph} (this work)")
        ax2.axhline(gt[ph], color=col, ls="--", lw=1.2, alpha=0.75)
    ax2.set_xlabel("slice index z")
    ax2.set_ylabel("volume fraction in slice")
    ax2.set_title("per-slice profile; dashed = published (Table S4)", fontsize=9)
    ax2.legend(fontsize=8, frameon=False)
    ax2.grid(alpha=0.25)
    ax2.set_ylim(0, 0.75)

    sub = "  ".join(f"{ph}: mine={mine[ph]:.4f} pub={gt[ph]:.3f}"
                    for ph in ("Ni", "YSZ", "pore"))
    fig.suptitle(f"{key}  ({folder})   {sub}", fontsize=10)
    fig.tight_layout()
    dest = os.path.join(OUT, f"phase2_slices_{key}.png")
    fig.savefig(dest, dpi=140)
    plt.close(fig)
    return dest


def main():
    print("=" * 78)
    print("PHASE 2 — volume fractions from voxel counts vs published Table S4")
    print("=" * 78)

    gtf = ground_truth_frame().set_index("sample")
    rows = []
    for key, folder, grain, state, *_ in SAMPLES:
        note = ZENODO_LABEL_NOTE[key]
        h = label_histogram(folder)
        counts = h["counts"]
        mapping = assign_labels(counts, note)
        mine = volume_fractions(counts, mapping)
        gt = {ph: float(gtf.loc[key, f"{ph if ph != 'pore' else 'Pore'}_Phi__T-S4"])
              for ph in ("Ni", "YSZ", "pore")}

        print(f"\n--- {key} ({folder}) ---")
        print(f"  brightness note used : {note}")
        print(f"  label -> phase       : pore={mapping['pore']}, "
              f"YSZ={mapping['YSZ']}, Ni={mapping['Ni']}, "
              f"unassigned={mapping['unassigned']}")
        una = sum(counts[v] for v in mapping["unassigned"])
        print(f"  unassigned voxels    : {una:,} "
              f"({100.0*una/counts.sum():.3f} %)")
        row = {"sample": key, "folder": folder, "grain": grain, "state": state,
               "label_pore": mapping["pore"], "label_YSZ": mapping["YSZ"],
               "label_Ni": mapping["Ni"],
               "unassigned_labels": str(mapping["unassigned"]),
               "unassigned_frac": float(una / counts.sum())}
        worst = 0.0
        for ph in ("Ni", "YSZ", "pore"):
            rel = abs(mine[ph] - gt[ph]) / gt[ph]
            worst = max(worst, rel)
            flag = "OK" if rel <= GATE_REL else "**FAIL**"
            print(f"  {ph:5s}  this work = {mine[ph]:.4f}   "
                  f"published = {gt[ph]:.3f}   rel.diff = {100*rel:5.2f} %  {flag}")
            row[f"{ph}_mine"] = mine[ph]
            row[f"{ph}_published"] = gt[ph]
            row[f"{ph}_rel_diff"] = rel
        row["worst_rel_diff"] = worst
        row["gate_pass"] = worst <= GATE_REL

        prof = per_slice_profile(folder, mapping)
        dest = save_slice_figure(key, folder, mapping, gt, mine, prof)
        prof.to_csv(os.path.join(OUT, f"phase2_profile_{key}.csv"), index=False)
        print(f"  [fig] {os.path.basename(dest)}")
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase2_volume_fractions.csv"), index=False)

    print("\n" + "=" * 78)
    print("PHASE 2 GATE")
    print("=" * 78)
    print(df[["sample", "Ni_mine", "Ni_published", "Ni_rel_diff",
              "worst_rel_diff", "gate_pass"]].to_string(index=False))
    allpass = bool(df["gate_pass"].all())
    print(f"\n  worst relative difference over all samples/phases: "
          f"{df['worst_rel_diff'].max()*100:.2f} %")
    print(f"  gate tolerance: {GATE_REL*100:.0f} % relative")
    print(f"  RESULT: {'PASS' if allpass else 'FAIL'}")
    return 0 if allpass else 1


if __name__ == "__main__":
    sys.exit(main())
