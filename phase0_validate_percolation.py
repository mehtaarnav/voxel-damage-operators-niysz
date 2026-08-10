"""
PHASE 0 — Validate the percolation-detection code on synthetic random media.

Purpose: before any real tomography is loaded, prove that the connected-component
+ face-spanning check in cmlib.percolation reproduces the known site-percolation
threshold of the simple cubic lattice, p_c = 0.3116077(2)
[Xu, Wang, Lv, Deng, Phys. Rev. E 89, 012120 (2014)].

GATE: detected threshold must be within ~0.02-0.03 of 0.3116.

CONVENTIONS USED HERE (see cmlib/percolation.py docstring for the full rationale)
  * adjacency        : 6-connectivity (face-sharing).  This is REQUIRED for the
                       0.3116 reference to apply; 18- and 26-connectivity have
                       thresholds 0.1372 and 0.09755 and are checked separately
                       as an independent confirmation that the structuring
                       elements are what we think they are.
  * spanning         : one component touching both index-0 and index-(n-1)
                       slices along the tested axis; free boundaries.
  * threshold estimate: the p at which the spanning probability crosses 0.5.
                       At finite L this is NOT p_c; it is displaced by
                       p50(L) - p_c ~ a * L^(-1/nu) with nu = 0.8774 (3D
                       percolation correlation-length exponent).  We therefore
                       ALSO run a finite-size scaling extrapolation to L->inf
                       and report that as the primary number.

Outputs (all written to out/phase0/):
  phase0_sweep_coarse.csv      spanning probability, coarse p grid
  phase0_sweep_fine.csv        spanning probability, fine p grid, 3 box sizes
  phase0_sweep_conn.csv        18- and 26-connectivity confirmation sweeps
  phase0_spanning_probability.png
  phase0_finite_size_scaling.png
  phase0_connectivity_check.png
  phase0_example_slices.png    visual sanity check of the synthetic media
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib.percolation import (  # noqa: E402
    SC_SITE_THRESHOLDS,
    label_phase,
    random_site_medium,
    spanning_labels,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "phase0")
os.makedirs(OUT, exist_ok=True)

P_C_REF = SC_SITE_THRESHOLDS[6]          # 0.3116077
NU_3D = 0.8774                            # 3D percolation correlation-length exponent
                                          # (Wang et al., Phys. Rev. E 87, 052107 (2013))


def spanning_all_axes(mask, connectivity=6):
    """Label once, then test face-to-face spanning along each of the 3 axes."""
    labels, n = label_phase(mask, connectivity)
    if n == 0:
        return (False, False, False)
    return tuple(spanning_labels(labels, ax).size > 0 for ax in (0, 1, 2))


def sweep(p_values, L, n_real, connectivity, seed):
    """Spanning probability vs p.  Returns a tidy DataFrame."""
    rng = np.random.default_rng(seed)
    rows = []
    for p in p_values:
        hits = np.zeros(3, dtype=int)
        for _ in range(n_real):
            m = random_site_medium((L, L, L), float(p), rng)
            hits += np.array(spanning_all_axes(m, connectivity), dtype=int)
        rows.append(
            dict(
                p=float(p),
                L=L,
                connectivity=connectivity,
                n_real=n_real,
                P_span_z=hits[0] / n_real,
                P_span_y=hits[1] / n_real,
                P_span_x=hits[2] / n_real,
                P_span_any_axis_mean=hits.sum() / (3 * n_real),
            )
        )
    return pd.DataFrame(rows)


def _erf_cdf(p, p50, w):
    return 0.5 * (1.0 + erf((p - p50) / (np.sqrt(2.0) * w)))


def estimate_p50(df, col="P_span_any_axis_mean"):
    """Estimate the 50% crossing two ways: linear interpolation and an erf fit."""
    d = df.sort_values("p")
    p = d["p"].to_numpy()
    y = d[col].to_numpy()

    # linear interpolation of the first upward 0.5 crossing
    p50_interp = np.nan
    for i in range(len(p) - 1):
        if y[i] <= 0.5 <= y[i + 1] and y[i + 1] > y[i]:
            p50_interp = p[i] + (0.5 - y[i]) * (p[i + 1] - p[i]) / (y[i + 1] - y[i])
            break

    # erf fit (models the finite-size rounded step)
    p50_fit, w_fit = np.nan, np.nan
    try:
        guess = p50_interp if np.isfinite(p50_interp) else float(np.median(p))
        popt, _ = curve_fit(_erf_cdf, p, y, p0=[guess, 0.01], maxfev=20000)
        p50_fit, w_fit = float(popt[0]), abs(float(popt[1]))
    except Exception as e:  # pragma: no cover
        print(f"    [erf fit failed: {e}]")
    return p50_interp, p50_fit, w_fit


def main():
    t0 = time.time()
    print("=" * 78)
    print("PHASE 0 — percolation-detection validation on synthetic random media")
    print("=" * 78)
    print(f"Reference simple-cubic site threshold (6-connectivity): p_c = {P_C_REF:.6f}")
    print("Conventions: 6-connectivity (face-sharing); face-to-face spanning;")
    print("             free boundaries; independent Bernoulli sites.\n")

    # ---------------------------------------------------------------- coarse
    print("[1/4] Coarse sweep  p = 0.15 .. 0.45 step 0.01,  L = 100,  6 realisations")
    p_coarse = np.round(np.arange(0.15, 0.4501, 0.01), 4)
    df_coarse = sweep(p_coarse, L=100, n_real=6, connectivity=6, seed=20260810)
    df_coarse.to_csv(os.path.join(OUT, "phase0_sweep_coarse.csv"), index=False)
    for _, r in df_coarse.iterrows():
        if 0.26 <= r.p <= 0.37:
            print(f"      p={r.p:.3f}  P_span(z,y,x)="
                  f"({r.P_span_z:.2f},{r.P_span_y:.2f},{r.P_span_x:.2f})")
    print(f"      -> {os.path.join('out','phase0','phase0_sweep_coarse.csv')}")

    # ------------------------------------------------------------------ fine
    print("\n[2/4] Fine sweep near the transition, three box sizes "
          "(finite-size scaling)")
    p_fine = np.round(np.arange(0.290, 0.3351, 0.0025), 5)
    fine_frames = []
    p50_by_L = {}
    for L, n_real in ((50, 40), (100, 24), (160, 12)):
        df = sweep(p_fine, L=L, n_real=n_real, connectivity=6, seed=770000 + L)
        fine_frames.append(df)
        pi, pf, wf = estimate_p50(df)
        p50_by_L[L] = (pi, pf, wf, n_real)
        print(f"      L={L:3d}  n_real={n_real:2d}   p50(interp)={pi:.4f}   "
              f"p50(erf fit)={pf:.4f}   width={wf:.4f}")
    df_fine = pd.concat(fine_frames, ignore_index=True)
    df_fine.to_csv(os.path.join(OUT, "phase0_sweep_fine.csv"), index=False)

    # finite-size extrapolation: p50(L) = p_c_inf + a * L^(-1/nu)
    Ls = np.array(sorted(p50_by_L), dtype=float)
    p50s = np.array([p50_by_L[int(L)][1] for L in Ls])   # use the erf-fit values
    x = Ls ** (-1.0 / NU_3D)
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, p50s, rcond=None)
    p_c_extrap, slope = float(coef[0]), float(coef[1])
    print(f"\n      Finite-size extrapolation  p50(L) = p_c + a*L^(-1/nu), nu={NU_3D}")
    print(f"        fitted a       = {slope:+.4f}")
    print(f"        p_c (L -> inf) = {p_c_extrap:.5f}   "
          f"(reference {P_C_REF:.5f}, difference {p_c_extrap - P_C_REF:+.5f})")

    # -------------------------------------------------- connectivity control
    print("\n[3/4] Independent control: 18- and 26-connectivity thresholds")
    conn_frames = []
    conn_results = {}
    for conn, grid in (
        (18, np.round(np.arange(0.120, 0.1601, 0.0025), 5)),
        (26, np.round(np.arange(0.085, 0.1151, 0.0025), 5)),
    ):
        df = sweep(grid, L=100, n_real=16, connectivity=conn, seed=990000 + conn)
        conn_frames.append(df)
        pi, pf, wf = estimate_p50(df)
        ref = SC_SITE_THRESHOLDS[conn]
        conn_results[conn] = (pf, ref)
        print(f"      {conn:2d}-connectivity: p50(erf fit)={pf:.4f}  "
              f"reference={ref:.4f}  difference={pf - ref:+.4f}")
    pd.concat(conn_frames, ignore_index=True).to_csv(
        os.path.join(OUT, "phase0_sweep_conn.csv"), index=False
    )

    # ----------------------------------------------------------------- plots
    print("\n[4/4] Writing figures")

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.plot(df_coarse.p, df_coarse.P_span_z, "o-", ms=4, lw=1.2,
            label="z axis")
    ax.plot(df_coarse.p, df_coarse.P_span_y, "s--", ms=3, lw=0.9, alpha=0.7,
            label="y axis")
    ax.plot(df_coarse.p, df_coarse.P_span_x, "^--", ms=3, lw=0.9, alpha=0.7,
            label="x axis")
    ax.axvline(P_C_REF, color="crimson", ls=":", lw=1.8,
               label=f"known $p_c$ = {P_C_REF:.4f}")
    ax.axhline(0.5, color="0.6", lw=0.7)
    ax.set_xlabel("site occupation probability $p$")
    ax.set_ylabel("spanning probability")
    ax.set_title("Phase 0: spanning probability vs $p$\n"
                 "simple cubic, 6-connectivity, $L$=100, 6 realisations")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase0_spanning_probability.png"), dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    for L in sorted(p50_by_L):
        d = df_fine[df_fine.L == L].sort_values("p")
        axes[0].plot(d.p, d.P_span_any_axis_mean, "o-", ms=3.5, lw=1.1,
                     label=f"L = {L}")
        pf = p50_by_L[L][1]
        pp = np.linspace(d.p.min(), d.p.max(), 300)
        axes[0].plot(pp, _erf_cdf(pp, pf, p50_by_L[L][2]), "-", lw=0.8, alpha=0.5,
                     color=axes[0].lines[-1].get_color())
    axes[0].axvline(P_C_REF, color="crimson", ls=":", lw=1.8,
                    label=f"$p_c$ = {P_C_REF:.4f}")
    axes[0].axhline(0.5, color="0.6", lw=0.7)
    axes[0].set_xlabel("$p$")
    axes[0].set_ylabel("spanning probability (mean over 3 axes)")
    axes[0].set_title("Fine sweep near the transition")
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].grid(alpha=0.25)

    xs = np.linspace(0, max(x) * 1.08, 100)
    axes[1].plot(x, p50s, "o", ms=7, color="C0", label="$p_{50}(L)$, erf fit")
    axes[1].plot(xs, p_c_extrap + slope * xs, "-", lw=1.2, color="C0",
                 label="linear fit")
    axes[1].axhline(P_C_REF, color="crimson", ls=":", lw=1.8,
                    label=f"known $p_c$ = {P_C_REF:.4f}")
    axes[1].plot(0, p_c_extrap, "*", ms=15, color="k",
                 label=f"extrapolated = {p_c_extrap:.4f}")
    for L, xv, yv in zip(Ls, x, p50s):
        axes[1].annotate(f"L={int(L)}", (xv, yv), textcoords="offset points",
                         xytext=(6, -10), fontsize=8)
    axes[1].set_xlabel(r"$L^{-1/\nu}$   ($\nu$ = %.4f)" % NU_3D)
    axes[1].set_ylabel("$p_{50}$")
    axes[1].set_title("Finite-size scaling extrapolation")
    axes[1].legend(frameon=False, fontsize=8, loc="upper left")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase0_finite_size_scaling.png"), dpi=150)
    plt.close(fig)

    df_conn = pd.concat(conn_frames, ignore_index=True)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for conn, colour in ((18, "C2"), (26, "C4")):
        d = df_conn[df_conn.connectivity == conn].sort_values("p")
        ax.plot(d.p, d.P_span_any_axis_mean, "o-", ms=4, lw=1.2, color=colour,
                label=f"{conn}-connectivity")
        ax.axvline(SC_SITE_THRESHOLDS[conn], color=colour, ls=":", lw=1.5)
    ax.axhline(0.5, color="0.6", lw=0.7)
    ax.set_xlabel("$p$")
    ax.set_ylabel("spanning probability")
    ax.set_title("Phase 0 control: 18- and 26-connectivity thresholds\n"
                 "(dotted lines = literature values 0.1372 and 0.09755)")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase0_connectivity_check.png"), dpi=150)
    plt.close(fig)

    # visual sanity check of the synthetic media themselves
    rng = np.random.default_rng(4242)
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6))
    for a, p in zip(axes, (0.20, 0.29, 0.32, 0.40)):
        m = random_site_medium((100, 100, 100), p, rng)
        lab, _ = label_phase(m, 6)
        span = spanning_labels(lab, 0)
        img = np.zeros(lab.shape[1:], dtype=float)
        sl = lab[:, :, 50]
        img = np.where(np.isin(sl, span), 2.0, np.where(sl > 0, 1.0, 0.0))
        a.imshow(img, cmap="viridis", interpolation="nearest", vmin=0, vmax=2)
        a.set_title(f"p = {p:.2f}\nspanning={span.size > 0}", fontsize=9)
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle("Phase 0: mid-plane slices — dark = empty, mid = non-spanning "
                 "cluster, bright = z-spanning cluster", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase0_example_slices.png"), dpi=150)
    plt.close(fig)

    # ------------------------------------------------------------------ gate
    print("\n" + "=" * 78)
    print("PHASE 0 GATE")
    print("=" * 78)
    dev_L100 = abs(p50_by_L[100][1] - P_C_REF)
    dev_extrap = abs(p_c_extrap - P_C_REF)
    print(f"  p50 at L=100 (raw, no extrapolation) : {p50_by_L[100][1]:.5f}  "
          f"|deviation| = {dev_L100:.5f}")
    print(f"  p_c extrapolated to L -> inf         : {p_c_extrap:.5f}  "
          f"|deviation| = {dev_extrap:.5f}")
    for conn in (18, 26):
        pf, ref = conn_results[conn]
        print(f"  {conn}-connectivity control              : {pf:.5f}  "
              f"|deviation| = {abs(pf - ref):.5f}  (ref {ref})")
    tol = 0.02
    ok = dev_L100 < tol and dev_extrap < tol
    print(f"\n  Gate tolerance: 0.02 (user-specified 0.02-0.03)")
    print(f"  RESULT: {'PASS' if ok else 'FAIL'}")
    print(f"\n  elapsed {time.time() - t0:.1f} s")
    print(f"  figures + csv in: {OUT}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
