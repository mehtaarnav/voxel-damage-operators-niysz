"""
PHASE 6 — assemble the comparison table and state the verdict.

Rows: the three anodes.  Columns: every pristine-state predictor (Phase 4),
every degraded-state outcome (Phase 5), and the matching published values.

STATISTICAL POSITION, STATED UP FRONT AND ENFORCED IN CODE
----------------------------------------------------------
n = 3.  With three items there are only 3! = 6 possible orderings, so ANY
predictor reproduces the outcome ordering with probability 1/6 = 17 % by chance
alone.  No correlation coefficient or p-value is computed or reported here,
because neither would be meaningful.  The only claim this design can support is
a DIRECTIONAL one: whether an ordering matches or does not.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.ground_truth import ground_truth_frame  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase6")
os.makedirs(OUT, exist_ok=True)
P4 = os.path.join(HERE, "out", "phase4")
P5 = os.path.join(HERE, "out", "phase5")
P2 = os.path.join(HERE, "out", "phase2")

GRAINS = ["fine", "medium", "coarse"]

# Predictors: (column, higher_is_better_for_retention, human label)
PREDICTORS = [
    ("lambda2_raw_mean",   True,  "algebraic connectivity lambda2 (raw)"),
    ("lambda2_norm_mean",  True,  "algebraic connectivity lambda2 (normalised)"),
    ("mincut_mean",        True,  "min-cut face-to-face"),
    ("g_eff_mean",         True,  "effective conductance"),
    ("neck_p10_nm_mean",   True,  "10th-percentile neck width"),
    ("Ni_phi",             True,  "Ni volume fraction (conventional)"),
    ("particle_d_volwt_nm", False, "mean Ni particle size, explicit watershed (conventional)"),
    ("chamber_equiv_diam_mean_nm_mean", False,
     "mean Ni chamber size, SNOW watershed (conventional)"),
    ("tpb_density_pre",    True,  "pristine TPB density (conventional)"),
]
# NOTE on direction: for the two particle-size predictors `higher_is_better` is
# set to False, i.e. they are ranked SMALLEST-FIRST.  That encodes the
# conventional expectation that a FINER microstructure is the better electrode.
# The results below show this expectation is inverted for retained Ni
# percolation, which is itself part of the finding -- so the ranking is also
# reported with the opposite direction in the output.

OUTCOMES = [
    ("P_span_retained",  "retained Ni spanning fraction"),
    ("P_reach_retained", "retained Ni reachable fraction"),
    ("tpb_retained",     "retained TPB density"),
    ("P_span_post",      "absolute post-redox spanning fraction"),
]


def read_csv(p):
    return pd.read_csv(p) if os.path.exists(p) else None


def rank_order(series, descending=True):
    """Grain names ordered best-first."""
    s = series.dropna()
    return list(s.sort_values(ascending=not descending).index)


def main():
    print("=" * 78)
    print("PHASE 6 — comparison table and verdict")
    print("=" * 78)

    gt = ground_truth_frame().set_index("sample")
    tbl = pd.DataFrame(index=GRAINS)

    # ---------- pristine predictors from Phase 4c ------------------------
    m = read_csv(os.path.join(P4, "phase4c_metrics_per_anode_8.0um.csv"))
    if m is not None:
        m["grain"] = m["sample"].str.replace("_pre", "", regex=False)
        m = m.set_index("grain")
        for c in m.columns:
            if c.endswith("_mean") or c.endswith("_sd"):
                tbl[c] = m[c]

    # ---------- conventional metrics -------------------------------------
    vf = read_csv(os.path.join(P2, "phase2_volume_fractions.csv"))
    if vf is not None:
        pre = vf[vf.state == "pristine"].set_index("grain")
        tbl["Ni_phi"] = pre["Ni_mine"]
        tbl["Ni_phi_published"] = pre["Ni_published"]

    pa = read_csv(os.path.join(P4, "phase4d_particles.csv"))
    if pa is not None:
        sel = pa[pa.min_distance == 4]
        tbl["particle_d_volwt_nm"] = sel.groupby("grain")["d_volweighted_nm"].mean()

    tp = read_csv(os.path.join(P4, "phase4b_tpb_full_stacks.csv"))
    if tp is not None:
        pre = tp[tp.state == "pristine"].set_index("grain")
        post = tp[tp.state == "degraded"].set_index("grain")
        tbl["tpb_density_pre"] = pre["tpb_density_um-2"]
        tbl["tpb_density_post"] = post["tpb_density_um-2"]
        tbl["tpb_retained"] = post["tpb_density_um-2"] / pre["tpb_density_um-2"]
        tbl["tpb_pre_published"] = pre["tpb_published_total"]
        tbl["tpb_retained_published"] = (post["tpb_published_total"]
                                         / pre["tpb_published_total"])

    # ---------- degraded outcomes from Phase 5 ---------------------------
    pc = read_csv(os.path.join(P5, "phase5_percolation.csv"))
    rt = read_csv(os.path.join(P5, "phase5_retention.csv"))
    if pc is not None:
        pre = pc[pc.state == "pristine"].set_index("grain")
        post = pc[pc.state == "degraded"].set_index("grain")
        tbl["P_span_pre"] = pre["P_span"]
        tbl["P_span_post"] = post["P_span"]
        tbl["P_reach_pre"] = pre["P_reach"]
        tbl["P_reach_post"] = post["P_reach"]
        tbl["percolates_post"] = post["percolates_x"]
        tbl["P_published_pre"] = pre["P_published"]
        tbl["P_published_post"] = post["P_published"]
    if rt is not None:
        rt = rt.set_index("grain")
        tbl["P_span_retained"] = rt["P_span_retained"]
        tbl["P_reach_retained"] = rt["P_reach_retained"]
        tbl["P_published_retained"] = rt["P_pub_retained"]

    tbl = tbl.reindex(GRAINS)
    tbl.to_csv(os.path.join(OUT, "phase6_comparison_table.csv"))

    with pd.option_context("display.width", 250, "display.max_columns", 200,
                           "display.float_format", lambda v: f"{v:,.4g}"):
        print("\nCOMPARISON TABLE")
        print(tbl.T.to_string())

    # ---------- ranking analysis -----------------------------------------
    print("\n" + "=" * 78)
    print("RANKING ANALYSIS  (best-preserved first)")
    print("=" * 78)

    outcome_orders = {}
    for col, lbl in OUTCOMES:
        if col in tbl.columns and tbl[col].notna().any():
            outcome_orders[col] = rank_order(tbl[col], descending=True)
            print(f"  OUTCOME  {lbl:42s}: {outcome_orders[col]}")
    if "P_published_retained" in tbl.columns:
        pub_order = rank_order(tbl["P_published_retained"], descending=True)
        print(f"  OUTCOME  {'published P retained':42s}: {pub_order}")
    if "tpb_retained_published" in tbl.columns:
        print(f"  OUTCOME  {'published TPB_total retained':42s}: "
              f"{rank_order(tbl['tpb_retained_published'], descending=True)}")

    print()
    rows = []
    for col, hib, lbl in PREDICTORS:
        if col not in tbl.columns or tbl[col].isna().all():
            continue
        order = rank_order(tbl[col], descending=hib)
        rev = list(reversed(order))
        rec = {"predictor": lbl, "column": col, "order": " > ".join(order),
               "order_reversed": " > ".join(rev)}
        for ocol, olbl in OUTCOMES:
            if ocol in outcome_orders:
                rec[f"matches::{ocol}"] = (order == outcome_orders[ocol])
                rec[f"matchesREV::{ocol}"] = (rev == outcome_orders[ocol])
        rows.append(rec)
        print(f"  PREDICTOR {lbl:52s}: {order}")
    rank_df = pd.DataFrame(rows)
    rank_df.to_csv(os.path.join(OUT, "phase6_rankings.csv"), index=False)

    print("\n" + "-" * 78)
    print("MATCH MATRIX (does the predictor's ordering equal the outcome's?)")
    print("-" * 78)
    mcols = [c for c in rank_df.columns if c.startswith("matches::")]
    if mcols:
        disp = rank_df[["predictor"] + mcols].copy()
        disp.columns = ["predictor"] + [c.split("::")[1] for c in mcols]
        print(disp.to_string(index=False))

    rcols = [c for c in rank_df.columns if c.startswith("matchesREV::")]
    if rcols:
        print("\n  ...and with the predictor's direction REVERSED "
              "(a True here means the predictor gets the ordering exactly "
              "BACKWARDS in its conventional direction):")
        disp = rank_df[["predictor"] + rcols].copy()
        disp.columns = ["predictor"] + [c.split("::")[1] for c in rcols]
        print(disp.to_string(index=False))

    print("\n" + "=" * 78)
    print("STATISTICAL CAVEAT (enforced)")
    print("=" * 78)
    print("  n = 3 anodes.  There are 3! = 6 possible orderings, so a predictor")
    print("  matches the outcome ordering with probability 1/6 = 17 % by chance.")
    print("  This is a DIRECTIONAL check only.  No correlation coefficient and")
    print("  no p-value is reported, because neither is meaningful at n = 3.")

    print(f"\n[saved] {os.path.join(OUT, 'phase6_comparison_table.csv')}")
    print(f"[saved] {os.path.join(OUT, 'phase6_rankings.csv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
