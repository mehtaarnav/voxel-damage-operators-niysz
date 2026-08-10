"""
E0b — saturation bridge, per the pre-registered rule in e0_vertical_slice.py
instruction D: "If all structures retain ~1.0 or collapse to ~0.0, report
saturation and adjust only the damage-intensity range, not the scientific
criterion."

The pre-registered run (n_rounds in {2,5,10}, p_erode=0.35, expand_vox=1 --
UNCHANGED here) saturated: 100% of structures retained P_span=1.0 at
n_rounds<=5 and 100% collapsed to P_span=0.0 at n_rounds=10. This is a
SEPARATE, clearly-labelled supplementary sweep at n_rounds in {6,7,8} (the
gap between the saturated regimes) -- it does not alter, overwrite, or
retroactively edit e0_vertical_slice.csv, and it does not change p_erode,
expand_vox, or any other D4 parameter. Same structures, same reconstruction,
same damage model; only the intensity GRID is extended, exactly as the
pre-registered rule permits.
"""

from __future__ import annotations

import os
import sys
import time

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import e0_vertical_slice as e0  # noqa: E402

OUT = e0.OUT
N_ROUNDS_BRIDGE = [6, 7, 8]     # the only thing that changes


def main():
    print("=" * 78)
    print("E0b SATURATION BRIDGE (supplementary to the pre-registered E0 run)")
    print("=" * 78)
    print(f"bridging n_rounds={N_ROUNDS_BRIDGE}, p_erode={e0.P_ERODE} "
          f"(UNCHANGED), expand_vox={e0.EXPAND_VOX} (UNCHANGED)\n")

    pilot_df = pd.read_csv(os.path.join(ROOT, "out", "next",
                                        "familyB_pilot.csv"))
    spacing = (e0.VOXEL_NM,) * 3
    t0 = time.time()
    centres, pairs, shape = e0.reconstruct_geometry()
    levels = [("base", 1.0), ("lower_tail", 1.5), ("lower_tail", 2.0)]

    rows = []
    for mode, tr in levels:
        for seed in e0.STRUCT_SEEDS:
            row = pilot_df[(pilot_df["mode"] == mode)
                           & (np.isclose(pilot_df["target_ratio"], tr))
                           & (pilot_df["seed"] == seed)].iloc[0]
            threshold = 0.0 if mode == "base" else float(row["intended_T_vox"])
            ni_mask, build, widths = e0.reconstruct_ni_mask(
                centres, pairs, shape, seed, threshold)
            e0.verify_against_pilot_csv(pilot_df, mode, tr, seed, build)
            achieved = 1.0 if mode == "base" else float(row["achieved_ratio"])

            vol, ysz_mask = e0.add_ysz_pore(ni_mask, seed=seed * 1000 + 7)
            m_pre = e0.measure_ternary(vol, spacing)

            for n_rounds in N_ROUNDS_BRIDGE:
                for dseed in e0.DAMAGE_SEEDS:
                    final_ni, ddiag = e0.apply_d4(
                        ni_mask, ysz_mask, n_rounds, e0.P_ERODE,
                        e0.EXPAND_VOX, dseed)
                    vol_post = np.where(
                        ysz_mask, e0.LABELS["YSZ"],
                        np.where(final_ni, e0.LABELS["Ni"], e0.LABELS["pore"])
                    ).astype(np.uint8)
                    m_post = e0.measure_ternary(vol_post, spacing)
                    rows.append(dict(
                        mode=mode, nominal_target_ratio=tr, struct_seed=seed,
                        achieved_p10_ratio=achieved, n_rounds=n_rounds,
                        p_erode=e0.P_ERODE, expand_vox=e0.EXPAND_VOX,
                        damage_seed=dseed,
                        phi_Ni_pre=m_pre["phi_Ni"], phi_Ni_post=m_post["phi_Ni"],
                        P_span_pre=m_pre["P_span"], P_span_post=m_post["P_span"],
                        P_reach_pre=m_pre["P_reach"], P_reach_post=m_post["P_reach"],
                        tpb_pre=m_pre["tpb_density"], tpb_post=m_post["tpb_density"],
                        lcf_pre=m_pre["largest_component_fraction"],
                        lcf_post=m_post["largest_component_fraction"],
                        retained_P_span=(m_post["P_span"] / m_pre["P_span"]
                                         if m_pre["P_span"] > 0 else np.nan),
                        retained_P_reach=(m_post["P_reach"] / m_pre["P_reach"]
                                          if m_pre["P_reach"] > 0 else np.nan),
                        retained_tpb=(m_post["tpb_density"] / m_pre["tpb_density"]
                                     if m_pre["tpb_density"] > 0 else np.nan),
                        **ddiag,
                    ))
            print(f"  {mode:10s} ratio={tr} seed={seed}: done  "
                  f"[{time.time()-t0:.0f}s]")

    df = pd.DataFrame(rows)
    dest = os.path.join(OUT, "e0b_saturation_bridge.csv")
    df.to_csv(dest, index=False)
    print(f"\n[saved] {dest}  total {time.time()-t0:.0f}s")

    print("\n" + "=" * 78)
    print("SATURATION CHECK (bridge intensities)")
    print("=" * 78)
    for nr in N_ROUNDS_BRIDGE:
        sub = df[df.n_rounds == nr]
        print(f"  n_rounds={nr:2d}: retained_P_span mean={sub.retained_P_span.mean():.3f} "
              f"min={sub.retained_P_span.min():.3f} max={sub.retained_P_span.max():.3f}  "
              f"frac~1.0(>0.98)={100*(sub.retained_P_span>0.98).mean():.0f}%  "
              f"frac~0.0(<0.02)={100*(sub.retained_P_span<0.02).mean():.0f}%  "
              f"frac_intermediate={100*((sub.retained_P_span>=0.02)&(sub.retained_P_span<=0.98)).mean():.0f}%")

    agg = df.groupby(["nominal_target_ratio", "n_rounds"]).agg(
        achieved_p10_ratio=("achieved_p10_ratio", "mean"),
        retained_P_span_mean=("retained_P_span", "mean"),
        retained_P_span_sd=("retained_P_span", "std"),
        retained_tpb_mean=("retained_tpb", "mean"),
        n=("retained_P_span", "count"),
    ).reset_index()
    print()
    with pd.option_context("display.width", 200):
        print(agg.to_string(index=False))

    # combine with original pre-registered run for one combined figure
    orig = pd.read_csv(os.path.join(OUT, "e0_vertical_slice.csv"))
    orig["source"] = "pre-registered"
    df["source"] = "bridge (post-hoc, intensity only)"
    combined = pd.concat([orig, df], ignore_index=True)
    combined.to_csv(os.path.join(OUT, "e0_combined_all_intensities.csv"),
                    index=False)

    fig, ax = plt.subplots(figsize=(9, 6))
    colors = {1.0: "C2", 1.5: "C0", 2.0: "C1"}
    for tr in sorted(combined.nominal_target_ratio.unique()):
        sub = combined[combined.nominal_target_ratio == tr]
        g = sub.groupby("n_rounds")["retained_P_span"].agg(["mean", "std"])
        g = g.sort_index()
        ax.errorbar(g.index, g["mean"], yerr=g["std"], fmt="o-",
                   color=colors.get(tr, "gray"), capsize=3,
                   label=f"nominal {tr}x (achieved ~"
                         f"{sub.achieved_p10_ratio.mean():.2f}x)")
    for nr in e0.N_ROUNDS_LEVELS:
        ax.axvline(nr, color="k", ls=":", lw=0.6, alpha=0.5)
    ax.set_xlabel("D4 damage intensity (n_rounds, p_erode=0.35 fixed)")
    ax.set_ylabel("retained P_span (post/pre), mean +/- sd over 5 seeds x 3 damage seeds")
    ax.set_title("E0 + E0b combined: retained P_span vs damage intensity\n"
                 "EXPLORATORY PIPELINE SPIKE -- NOT CONFIRMATORY", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "e0_combined_retained_pspan_vs_intensity.png"),
               dpi=145)
    print(f"\n[saved] {os.path.join(OUT, 'e0_combined_retained_pspan_vs_intensity.png')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
