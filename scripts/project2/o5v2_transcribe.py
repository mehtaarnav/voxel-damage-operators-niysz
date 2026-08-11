"""Transcribe the O5v2 gate results from the two committed run reports into a
CSV, so that every number appearing in a manuscript figure has a tabular source.

The O5v2 runs wrote their results to markdown rather than to CSV. Nothing here
is recomputed: the values are copied verbatim from

    out/project2/O5v2_RESULTS/o5v2_report.md          (Option A, curvature-ranked)
    out/project2/O5v2_RESULTS/o5v2_optionB_report.md  (Option B, greedy dA<=0)

and the `source` column records which. Pristine S_spec = 0.45052 and pristine
neck volume = 63 voxels are the n = 0 rows, identical in both reports.
"""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "out", "project2", "o5v2_area_barrier.csv")

A = "o5v2_report.md"
B = "o5v2_optionB_report.md"

ROWS = [
    # operator, stencil, n, vol_error, S_spec, neck_voxels, accept_rate, source
    ("curvature_ranked", "6",  0, 0.0, 0.45052, 63, "", A),
    ("curvature_ranked", "6",  1, 0.0, 0.45195, 63, "", A),
    ("curvature_ranked", "6",  3, 0.0, 0.45219, 63, "", A),
    ("curvature_ranked", "6",  5, 0.0, 0.45362, 63, "", A),
    ("curvature_ranked", "26", 0, 0.0, 0.45052, 63, "", A),
    ("curvature_ranked", "26", 1, 0.0, 0.45195, 57, "", A),
    ("curvature_ranked", "26", 3, 0.0, 0.44431, 15, "", A),
    ("curvature_ranked", "26", 5, 0.0, 0.44216,  0, "", A),
    ("greedy_area",      "26", 0, 0.0, 0.45052, 63, 0.0, B),
    ("greedy_area",      "26", 1, 0.0, 0.45052, 63, 0.0, B),
    ("greedy_area",      "26", 3, 0.0, 0.45052, 63, 0.0, B),
    ("greedy_area",      "26", 5, 0.0, 0.45052, 63, 0.0, B),
    ("greedy_area",      "26", 8, 0.0, 0.45052, 63, 0.0, B),
]

with open(OUT, "w", newline="", encoding="utf8") as fh:
    w = csv.writer(fh)
    w.writerow(["operator", "stencil", "n_rounds", "volume_error",
                "S_spec", "neck_voxels", "acceptance_rate", "source_report"])
    w.writerows(ROWS)
print("wrote", OUT)
