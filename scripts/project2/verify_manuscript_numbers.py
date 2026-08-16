"""Check numbers quoted in the manuscript against the tables they come from.

This exists because the failure this project keeps hitting is not a wrong
calculation, it is a correct calculation whose result was updated in one place
and not another: a figure regenerated while its caption kept the old value, a
claim corrected in the primer but not the paper, an annotation asserting a
result that had been retracted.

Each check pairs a value recomputed from a committed CSV with a literal string
that must appear in manuscript.tex. Run by CI on every push.
"""
import os
import re
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P2 = os.path.join(ROOT, "out", "project2")
TEX = open(os.path.join(ROOT, "out", "writeup", "manuscript.tex"),
           encoding="utf8").read()

# The primer, the README and the HTML build quote the same quantities as the
# paper. They drifted once, because only the manuscript was checked here, so
# they are checked too: a value corrected in one document and not the others is
# exactly the failure this script exists to catch.
DOCS = {}
for name, rel in (("primer", ("out", "writeup", "PRIMER_voxel_operators.md")),
                  ("readme", ("README.md",)),
                  ("html", ("out", "writeup", "primer.html"))):
    try:
        DOCS[name] = open(os.path.join(ROOT, *rel), encoding="utf8").read()
    except FileNotFoundError:
        DOCS[name] = ""

checks = []


def check(name, condition, literal=None):
    """condition: recomputed fact is as expected. literal: must appear in tex."""
    in_tex = True if literal is None else (literal in TEX)
    checks.append((name, bool(condition), in_tex, literal))


# ---- trajectory: the gated quantity falls, TPB rises -----------------------
tr = pd.read_csv(f"{P2}/o7_trajectory.csv")
check("TPB rises to 4.99x", round(tr.tpb_ratio.iloc[-1], 2) == 4.99, "4.99")
check("S_spec falls ~3.6%", abs(tr.s_ratio.iloc[-1] - 0.9638) < 5e-4)

# ---- the area-neutral plateau ---------------------------------------------
hs = pd.read_csv(f"{P2}/o7_da_histogram.csv")
neutral = float(hs.loc[hs.dA == 0, "share"].iloc[0])
check("plateau is ~80%", 0.795 < neutral < 0.81)
st = pd.read_csv(f"{P2}/o7_move_stats.csv").iloc[0]
check("plateau count 207369 of 258870",
      int(st.neutral) == 207369 and int(st.accepted) == 258870,
      r"207\,369")

# ---- counterfactual: the contact is the whole effect ----------------------
cf = pd.read_csv(f"{P2}/o7_counterfactual.csv").set_index("case")
check("contact-adjacent alone gives 5.03",
      round(cf.loc["contact_adjacent_only", "tpb_ratio"], 2) == 5.03, "5.03")
check("away-from-contact gives 1.04",
      round(cf.loc["away_from_contact_only", "tpb_ratio"], 2) == 1.04, "1.04")
check("93.0% of removed voxels are next to YSZ",
      abs(st.removed_near_ysz / st.removed - 0.930) < 1e-3,
      r"\SI{93.0}{\percent}")

# ---- threshold sweep: fine is never first to fail -------------------------
th = pd.read_csv(f"{P2}/o7_threshold_transitions.csv")
first = th.groupby(["threshold", "anode"]).transition.mean().unstack().idxmin(axis=1)
check("fine is first to fail at no threshold", (first != "fine").all())
check("nine thresholds were swept", th.threshold.nunique() == 9)

# ---- batched greedy accepts and raises area -------------------------------
ab = pd.read_csv(f"{P2}/o5v2_area_barrier.csv")
g = ab[ab.operator == "greedy_area"]
check("greedy acceptance ~96%",
      abs(float(g[g.n_rounds == 1].acceptance_rate.iloc[0]) - 0.9672) < 1e-3,
      r"\SI{96}{\percent}")
check("greedy S_spec rises to 0.46365",
      float(g[g.n_rounds == 1].S_spec.iloc[0]) == 0.46365, "0.46365")
check("pristine S_spec 0.45052",
      float(g[g.n_rounds == 0].S_spec.iloc[0]) == 0.45052, "0.45052")

# ---- erosion TPB excursion ------------------------------------------------
rg = pd.read_csv(f"{P2}/c1real_rni_gate.csv")
rg["ratio"] = rg.tpb_um2 / rg.tpb_pristine
peak = rg.loc[rg.groupby("anode").ratio.idxmax()].set_index("anode").ratio
check("erosion peak 14.8 (fine)", round(peak["fine"], 1) == 14.8, "14.8")
check("erosion peak 12.3 (medium)", round(peak["medium"], 1) == 12.3, "12.3")
check("erosion peak 7.4 (coarse)", round(peak["coarse"], 1) == 7.4, "7.4")
n8 = rg[rg.n_rounds == 8].ratio
check("collapse range 4.8-10.2%",
      round(n8.min() * 100, 1) == 4.8 and round(n8.max() * 100, 1) == 10.2,
      r"\SIrange{4.8}{10.2}{\percent}")

# ---- 26-connectivity check ------------------------------------------------
cc = pd.read_csv(f"{P2}/o7_connectivity_check.csv")
gain = cc.groupby("anode").gain.mean()
check("fine gains most from 26-connectivity", gain.idxmax() == "fine")
check("fine gain is +0.0004", round(gain["fine"], 4) == 0.0004, "0.0004")
check("medium and coarse gain nothing",
      round(gain["medium"], 4) == 0 and round(gain["coarse"], 4) == 0)

# ---- superseded values must not survive anywhere ---------------------------
for stale, where in (("5.07", "counterfactual contact ratio"),
                     ("208,180", "plateau count"),
                     (r"208\,180", "plateau count"),
                     ("92.8", "removed-adjacent-to-YSZ share"),
                     ("0.00955", "tie-break contrast"),
                     ("factor of 40", "tie-break factor"),
                     ("22.5346", "fine TPB endpoint")):
    for doc, text in list(DOCS.items()) + [("manuscript", TEX)]:
        checks.append((f"superseded {stale} absent from {doc} ({where})",
                       stale not in text, True, None))

# ---- claims that must NOT reappear ----------------------------------------
for dead in ("accepts no move", "accepts zero", "no move accepted",
             "the remaining candidate"):
    checks.append((f"retracted phrase absent: '{dead}'",
                   dead not in TEX, True, None))

# ---------------------------------------------------------------------------
bad = [c for c in checks if not (c[1] and c[2])]
w = max(len(c[0]) for c in checks) + 2
for name, ok_val, ok_tex, lit in checks:
    status = "ok" if (ok_val and ok_tex) else "FAIL"
    detail = "" if (ok_val and ok_tex) else (
        "  <- value mismatch" if not ok_val else f"  <- {lit!r} not in tex")
    print(f"{status:>4}  {name:<{w}}{detail}")
print(f"\n{len(checks) - len(bad)}/{len(checks)} checks pass")
sys.exit(1 if bad else 0)
