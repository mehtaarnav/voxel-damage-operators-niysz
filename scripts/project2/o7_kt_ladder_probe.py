"""Does finite T* produce Rayleigh neck thinning, or just surface roughening?

This adjudicates whether the kT ladder (Option 1) is worth running at all.

Geometry: the O5v2 test case -- two Ni spheres joined by a straight cylinder,
in pore. Two-phase limit (no YSZ), so dE = 2*(nN(a) - nN(b)) exactly, which is
the identity verified in o7_derivation_check.py.

The instability the operator is supposed to model is Rayleigh-type neck
break-up: the NECK thins while the spheres stay compact. The failure mode we
are worried about is that finite temperature instead produces uniform SURFACE
ROUGHENING -- the neck is destroyed by generic thermal erosion, not by the
collective long-wavelength instability.

Discriminator:
  neck_min   -- minimum cross-sectional Ni area along the neck axis
  S_total    -- total Ni surface area (roughening inflates this)
  If the neck thins while S_total stays flat  -> real instability. Ladder works.
  If S_total inflates and the neck follows    -> roughening. Ladder is worthless.
"""
import numpy as np
import scipy.ndimage as ndi

ST7 = ndi.generate_binary_structure(3, 1)     # centre INCLUDED (sum 7)
STRUCT6 = ST7.copy()
STRUCT6[1, 1, 1] = False                      # centre EXCLUDED (sum 6)
# NOTE: cmlib/damage2.py:31 uses the centre-INCLUDED form for nb, which puts a
# +1 offset on Ni sites only. The dE identity requires the centre-excluded
# count. See o7_struct6_offset_check.py.


def make_dumbbell(n=64, R=13.0, r_neck=4.0, sep=26.0):
    """Two spheres joined by a cylinder along x, centred in an n^3 box."""
    z, y, x = np.ogrid[:n, :n, :n]
    c = n / 2.0
    s1 = (z - c) ** 2 + (y - c) ** 2 + (x - (c - sep / 2)) ** 2 <= R ** 2
    s2 = (z - c) ** 2 + (y - c) ** 2 + (x - (c + sep / 2)) ** 2 <= R ** 2
    cyl = (((z - c) ** 2 + (y - c) ** 2) <= r_neck ** 2) & \
          (np.abs(x - c) <= sep / 2)
    return (s1 | s2 | cyl)


def surface_area_vox(ni):
    """Exposed Ni faces = count of Ni-nonNi 6-neighbour bonds."""
    nb = ndi.convolve(ni.astype(np.int16), STRUCT6.astype(np.int16),
                      mode="constant", cval=0)
    return int((6 - nb)[ni].sum())


def neck_profile(ni):
    """Cross-sectional Ni area per x-slice; the neck is the interior minimum."""
    areas = ni.sum(axis=(0, 1))
    n = len(areas)
    mid = slice(n // 2 - 6, n // 2 + 7)
    return int(areas[mid].min()), areas


def run_kmc(ni0, t_star, n_sweeps, seed, moves_per_sweep_frac=0.05):
    """Metropolis single-voxel volume-conserving swap at dimensionless T*.

    dE = 2*(nN(a) - nN(b)) in units of J_NP. Interior sites only -- the
    identity requires exactly 6 neighbours (o7_derivation_check.py).
    """
    ni = ni0.copy()
    rng = np.random.default_rng(seed)
    n, proposed, accepted = ni.shape[0], 0, 0

    interior = np.zeros_like(ni)
    interior[1:-1, 1:-1, 1:-1] = True

    for _ in range(n_sweeps):
        nb = ndi.convolve(ni.astype(np.int16), STRUCT6.astype(np.int16),
                          mode="constant", cval=0)
        surf = ni & ~ndi.binary_erosion(ni, structure=ST7) & interior
        front = (~ni) & ndi.binary_dilation(ni, structure=ST7) & interior
        si = np.flatnonzero(surf.ravel())
        fi = np.flatnonzero(front.ravel())
        if si.size == 0 or fi.size == 0:
            break
        k = max(1, int(moves_per_sweep_frac * si.size))
        a_pick = si[rng.integers(0, si.size, size=k)]
        b_pick = fi[rng.integers(0, fi.size, size=k)]
        nbf = nb.ravel()
        flat = ni.ravel()
        shp = ni.shape
        for ia, ib in zip(a_pick, b_pick):
            ia, ib = int(ia), int(ib)
            if not flat[ia] or flat[ib]:
                continue                      # stale after earlier accepts
            za, ya, xa = np.unravel_index(ia, shp)
            zb, yb, xb = np.unravel_index(ib, shp)
            if abs(za - zb) + abs(ya - yb) + abs(xa - xb) <= 1:
                continue                      # adjacency breaks the algebra
            dE = 2.0 * (int(nbf[ia]) - int(nbf[ib]))
            proposed += 1
            if dE <= 0:
                acc = True
            elif t_star <= 0:
                acc = False
            else:
                acc = rng.random() < np.exp(-dE / t_star)
            if acc:
                flat[ia] = False
                flat[ib] = True
                accepted += 1
    return ni, proposed, accepted


def main():
    ni0 = make_dumbbell()
    v0 = int(ni0.sum())
    s0 = surface_area_vox(ni0)
    neck0, _ = neck_profile(ni0)
    print(f"pristine: V={v0}  S={s0}  neck_min={neck0}\n")

    print(f"{'T*':>8} {'acc_rate':>10} {'dV':>6} {'S/S0':>8} "
          f"{'neck/neck0':>11} {'verdict':>22}")
    print("-" * 70)
    for t_star in [0.0, 0.01, 0.05, 0.10, 0.20, 0.50, 1.0, 2.0]:
        ni, prop, acc = run_kmc(ni0, t_star, n_sweeps=40, seed=7)
        v1 = int(ni.sum())
        s1 = surface_area_vox(ni)
        neck1, _ = neck_profile(ni)
        rate = acc / max(prop, 1)
        s_ratio = s1 / s0
        n_ratio = neck1 / max(neck0, 1)
        # Roughening inflates S. A real instability thins the neck at ~flat S.
        if rate < 1e-6:
            verdict = "FROZEN (no moves)"
        elif s_ratio > 1.10 and n_ratio < 0.9:
            verdict = "ROUGHENING"
        elif s_ratio <= 1.10 and n_ratio < 0.9:
            verdict = "NECK THINNING"
        elif s_ratio > 1.10:
            verdict = "ROUGHENING (neck held)"
        else:
            verdict = "little change"
        print(f"{t_star:8.2f} {rate:10.4f} {v1-v0:6d} {s_ratio:8.3f} "
              f"{n_ratio:11.3f} {verdict:>22}")

    print("\nInterpretation:")
    print("  A usable ladder band requires NECK THINNING at flat surface area.")
    print("  If every moving rung is ROUGHENING, finite T* does not reach the")
    print("  Rayleigh instability and the ladder cannot test the mechanism.")


if __name__ == "__main__":
    main()
