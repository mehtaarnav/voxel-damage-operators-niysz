"""
PHASE 4e — is lambda2_raw measuring connectivity, or just network size?

The report claims that lambda2_raw tracks inverse node count.  That claim should
be measured, not asserted.  We fit lambda2 ~ C * N^(-a) across all 21 ROIs (all
three anodes pooled) and report the exponent and the scatter of lambda2 * N^a.

If a ~ 1 and lambda2*N is roughly constant across anodes, then lambda2_raw
carries essentially no information beyond the number of watershed chambers, and
its apparent ranking of the anodes is a restatement of their coarseness.
"""

from __future__ import annotations

import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase4")

df = pd.read_csv(os.path.join(OUT, "phase4c_metrics_per_roi_8.0um.csv"))
df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["lambda2_raw", "n_nodes"])
df = df[df.lambda2_raw > 0]

x = np.log(df["n_nodes"].to_numpy(dtype=float))
y = np.log(df["lambda2_raw"].to_numpy(dtype=float))
a, b = np.polyfit(x, y, 1)          # log l2 = a*log N + b
resid = y - (a * x + b)
r2 = 1.0 - resid.var() / y.var()

print("=" * 74)
print("PHASE 4e — lambda2 vs network size, all 21 ROIs pooled")
print("=" * 74)
print(f"  fit: lambda2 = exp({b:.3f}) * N^({a:.3f})")
print(f"       i.e. exponent a = {a:.3f}   (a = -1 means pure 1/N)")
print(f"       R^2 (in log space) = {r2:.3f}")
print(f"       prefactor C = {np.exp(b):.1f}")

print("\n  Per-anode means of lambda2 * N (constant => lambda2 is just 1/N):")
df["l2_times_N"] = df["lambda2_raw"] * df["n_nodes"]
for s, d in df.groupby("sample"):
    print(f"    {s:12s} n_roi={len(d):2d}  lambda2={d.lambda2_raw.mean():8.3f}  "
          f"N={d.n_nodes.mean():6.1f}  lambda2*N = {d.l2_times_N.mean():8.1f} "
          f"+/- {d.l2_times_N.std(ddof=1):7.1f}")

print("\n  Spread of the per-anode lambda2 means      : "
      f"{df.groupby('sample').lambda2_raw.mean().max() / df.groupby('sample').lambda2_raw.mean().min():.2f}x")
print("  Spread of the per-anode lambda2*N means    : "
      f"{df.groupby('sample').l2_times_N.mean().max() / df.groupby('sample').l2_times_N.mean().min():.2f}x")
print("\n  If the second number is much smaller than the first, multiplying by N")
print("  has removed most of the between-anode difference, i.e. lambda2_raw was")
print("  largely reporting network size rather than connectivity.")

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
for s, d in df.groupby("sample"):
    axes[0].scatter(d.n_nodes, d.lambda2_raw, s=40, label=s)
    axes[1].scatter(d.n_nodes, d.l2_times_N, s=40, label=s)
nn = np.linspace(df.n_nodes.min(), df.n_nodes.max(), 200)
axes[0].plot(nn, np.exp(b) * nn ** a, "k--", lw=1.3,
             label=f"fit  $\\lambda_2 \\propto N^{{{a:.2f}}}$")
axes[0].set_xscale("log"); axes[0].set_yscale("log")
axes[0].set_xlabel("network nodes N"); axes[0].set_ylabel(r"$\lambda_2$ raw")
axes[0].set_title(r"$\lambda_2$ falls with network size")
axes[1].set_xscale("log")
axes[1].set_xlabel("network nodes N"); axes[1].set_ylabel(r"$\lambda_2 \times N$")
axes[1].set_title(r"$\lambda_2 \times N$ — near-constant across anodes")
for ax in axes:
    ax.grid(alpha=0.25); ax.legend(fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "phase4e_lambda2_scaling.png"), dpi=145)
print(f"\n[saved] {os.path.join(OUT, 'phase4e_lambda2_scaling.png')}")
