"""Transcribe implementation-validation constants that were reported in
markdown rather than written to a table, so every number quoted in the
manuscript has a tabular source.

Nothing is recomputed here. The percolation-threshold recovery is the
finite-size extrapolation of the sweeps in out/phase0 (phase0_sweep_fine.csv,
phase0_sweep_coarse.csv, phase0_sweep_conn.csv), reported in REPORT.md; the
TPB wrap factor is from the analytic single-triple-line unit test recorded in
out/project2/STEP2_REPORT.md. The `source` column records where each came
from.
"""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "out", "project2", "validation_constants.csv")

ROWS = [
    dict(quantity="sc_site_percolation_threshold",
         measured=0.31218,
         reference=0.3116077,
         deviation=0.0006,
         units="site occupation probability",
         note="6-connectivity, finite-size extrapolated; 18- and "
              "26-connectivity controls within 0.0013",
         source="REPORT.md gate 0; sweeps in out/phase0/"),
    dict(quantity="tpb_periodic_wrap_overcount",
         measured=4.0,
         reference=1.0,
         deviation=3.0,
         units="ratio",
         note="analytic single-triple-line case; wrapping gather over-counted "
              "by exactly 4x before the fix",
         source="out/project2/STEP2_REPORT.md section 1"),
    dict(quantity="volume_fraction_worst_deviation",
         measured=0.00196,
         reference=0.0,
         deviation=0.00196,
         units="relative",
         note="worst of 18 values, 6 stacks x 3 phases, against published "
              "volume fractions",
         source="out/phase2/phase2_volume_fractions.csv column worst_rel_diff"),
]

with open(OUT, "w", newline="", encoding="utf8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ROWS[0]))
    w.writeheader()
    w.writerows(ROWS)
print("wrote", OUT)
