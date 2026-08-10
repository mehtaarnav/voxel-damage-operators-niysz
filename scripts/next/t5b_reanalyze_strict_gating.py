"""
Re-analyze the already-collected T5b data (t5b_deviations.csv) against the
STRICT PER-SEED gating criteria frozen in preregistration.md #0c amendment B,
which supersedes T5b's original aggregate-mean check (that check silently
passed target_ratio=2.5 via sign cancellation across seeds -- see
out/next/t5b_coupling_decision_report.md for the diagnosis). No aggregate
mean is used as a pass criterion anywhere in this script; every verdict is
per-seed, then rolled up by counting passes.

Does not re-run the experiment -- t5b_coupling_experiment.py already produced
out/next/t5b_deviations.csv; this only re-applies the corrected criteria.
"""
from __future__ import annotations

import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "out", "next")

BASE_N_NODES = 125


def gate(row, target_min):
    checks = {
        "p10_ge_target": row.achieved_p10_ratio >= target_min,
        "phi_le_5pct": abs(row.phi_dev_pct) <= 5.0,
        "phi_le_2pct_target": abs(row.phi_dev_pct) <= 2.0,
        "cpsd_le_5pct": abs(row.cpsd_dev_pct) <= 5.0,
        "p50_ratio_le_1_15": row.neck_p50_ratio <= 1.15,
        "n_nodes_ge_95pct": row.n_nodes >= 0.95 * BASE_N_NODES,
        "P_span_intact": row.P_span >= 0.90,
    }
    full_pass = (checks["p10_ge_target"] and checks["phi_le_5pct"]
                and checks["cpsd_le_5pct"] and checks["p50_ratio_le_1_15"]
                and checks["n_nodes_ge_95pct"] and checks["P_span_intact"])
    return checks, full_pass


def main():
    d = pd.read_csv(os.path.join(OUT, "t5b_deviations.csv"))
    rows = []
    print("=" * 100)
    print("T5b STRICT PER-SEED GATING (preregistration.md #0c amendment B)")
    print("=" * 100)
    for mode in ("lower_tail", "uniform"):
        for tr in (1.5, 2.0, 2.5):
            sub = d[(d["mode"] == mode) & (d["target_ratio"] == tr)]
            print(f"\n--- {mode} target_ratio={tr} ---")
            n_pass = 0
            for _, r in sub.iterrows():
                checks, full = gate(r, tr)
                n_pass += int(full)
                failed = [k for k, v in checks.items() if not v]
                print(f"  seed={r.seed}  p10={r.achieved_p10_ratio:.2f}  "
                      f"phi_dev={r.phi_dev_pct:+.2f}%  "
                      f"cpsd_dev={r.cpsd_dev_pct:+.2f}%  "
                      f"p50_ratio={r.neck_p50_ratio:.2f}  "
                      f"n_nodes={r.n_nodes}  P_span={r.P_span:.3f}  "
                      f"-> {'PASS' if full else 'fail'}"
                      + (f"  (failed: {', '.join(failed)})" if failed else ""))
                rows.append(dict(mode=mode, target_ratio=tr, seed=r.seed,
                                 achieved_p10_ratio=r.achieved_p10_ratio,
                                 phi_dev_pct=r.phi_dev_pct,
                                 cpsd_dev_pct=r.cpsd_dev_pct,
                                 p50_ratio=r.neck_p50_ratio,
                                 n_nodes=r.n_nodes, P_span=r.P_span,
                                 full_pass=full,
                                 failed_checks=",".join(failed)))
            feasible = n_pass >= 4
            print(f"  {mode} ratio={tr}: {n_pass}/5 seeds pass -> "
                  f"{'FEASIBLE' if feasible else 'not feasible'}")

    out_df = pd.DataFrame(rows)
    dest = os.path.join(OUT, "t5b_strict_gating.csv")
    out_df.to_csv(dest, index=False)
    print(f"\n[saved] {dest}")

    print("\n" + "=" * 100)
    print("SUMMARY (per preregistration.md #0c amendment A: primary envelope is 1.5x-2.0x)")
    print("=" * 100)
    summary = out_df.groupby(["mode", "target_ratio"])["full_pass"].sum().reset_index()
    summary["feasible"] = summary["full_pass"] >= 4
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
