"""Numerical verification of the three-phase Potts dE identity (O7).

Claim:  dE = J_NP * [ 2*(nN(a) - nN(b)) + (1 + cos_t)*(nY(a) - nY(b)) ]

for a volume-conserving swap: Ni at a -> pore, pore at b -> Ni, with a and b
non-adjacent. Verified by brute-force bond enumeration on random 3-phase
lattices. Also checks the reduction to the two-phase O5v2 identity and the
theta = 0 / 90 / 180 limits.
"""
import numpy as np

NI, YSZ, PORE = 0, 1, 2
OFFS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]


def total_energy(lat, J):
    """Brute-force: sum J over all 6-neighbour bonds, each counted once."""
    E = 0.0
    nz, ny, nx = lat.shape
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                for dz, dy, dx in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
                    z2, y2, x2 = z + dz, y + dy, x + dx
                    if z2 < nz and y2 < ny and x2 < nx:
                        E += J[lat[z, y, x], lat[z2, y2, x2]]
    return E


def neigh_counts(lat, site):
    """(nN, nY) among the 6 neighbours of site; out-of-bounds ignored."""
    z, y, x = site
    nz, ny, nx = lat.shape
    nN = nY = 0
    for dz, dy, dx in OFFS:
        z2, y2, x2 = z + dz, y + dy, x + dx
        if 0 <= z2 < nz and 0 <= y2 < ny and 0 <= x2 < nx:
            if lat[z2, y2, x2] == NI:
                nN += 1
            elif lat[z2, y2, x2] == YSZ:
                nY += 1
    return nN, nY


def make_J(J_NP, cos_t):
    """J_YP - J_NY = J_NP*cos_t. One gauge choice: J_NY = J_NP, so
    J_YP = J_NP*(1 + cos_t). Same-phase bonds zero."""
    J_NY = J_NP
    J_YP = J_NP * (1.0 + cos_t)
    J = np.zeros((3, 3))
    J[NI, PORE] = J[PORE, NI] = J_NP
    J[NI, YSZ] = J[YSZ, NI] = J_NY
    J[YSZ, PORE] = J[PORE, YSZ] = J_YP
    return J


def predicted(lat, a, b, J_NP, cos_t):
    nNa, nYa = neigh_counts(lat, a)
    nNb, nYb = neigh_counts(lat, b)
    return J_NP * (2.0 * (nNa - nNb) + (1.0 + cos_t) * (nYa - nYb))


def manhattan(a, b):
    return sum(abs(int(p) - int(q)) for p, q in zip(a, b))


def interior(sites, shape):
    """Sites with all 6 neighbours in bounds. The identity substitutes
    nP = 6 - nN - nY, which requires exactly 6 neighbours; boundary sites
    have fewer and are outside the derivation's scope."""
    nz, ny, nx = shape
    keep = ((sites[:, 0] > 0) & (sites[:, 0] < nz - 1) &
            (sites[:, 1] > 0) & (sites[:, 1] < ny - 1) &
            (sites[:, 2] > 0) & (sites[:, 2] < nx - 1))
    return sites[keep]


def run(n_trials=4000, seed=20260812):
    rng = np.random.default_rng(seed)
    J_NP = 1.0
    worst = {}
    for cos_t, name in [(1.0, "theta=0 (wetting)"),
                        (0.0, "theta=90"),
                        (-np.cos(np.deg2rad(50)), "theta=130 (Bjornsson)"),
                        (-1.0, "theta=180 (dewetting)")]:
        J = make_J(J_NP, cos_t)
        max_err = 0.0
        n_ok = 0
        for _ in range(n_trials):
            lat = rng.integers(0, 3, size=(5, 5, 5)).astype(np.int8)
            ni_sites = interior(np.argwhere(lat == NI), lat.shape)
            pore_sites = interior(np.argwhere(lat == PORE), lat.shape)
            if len(ni_sites) == 0 or len(pore_sites) == 0:
                continue
            a = tuple(ni_sites[rng.integers(len(ni_sites))])
            b = tuple(pore_sites[rng.integers(len(pore_sites))])
            if manhattan(a, b) <= 1:
                continue          # non-adjacency condition, as in O5v2
            E0 = total_energy(lat, J)
            lat2 = lat.copy()
            lat2[a] = PORE
            lat2[b] = NI
            E1 = total_energy(lat2, J)
            actual = E1 - E0
            pred = predicted(lat, a, b, J_NP, cos_t)
            max_err = max(max_err, abs(actual - pred))
            n_ok += 1
        worst[name] = (max_err, n_ok)
        print(f"{name:28s}  n={n_ok:5d}  max|actual-pred| = {max_err:.3e}")

    # Reduction to the two-phase O5v2 identity: no YSZ present at all.
    print("\nReduction check (YSZ absent) vs O5v2  dA = 2*(nb(a)-nb(b)):")
    J = make_J(1.0, -np.cos(np.deg2rad(50)))
    max_err = 0.0
    for _ in range(2000):
        lat = rng.integers(0, 2, size=(5, 5, 5)).astype(np.int8)
        lat[lat == 1] = PORE      # only Ni and pore
        ni_sites = interior(np.argwhere(lat == NI), lat.shape)
        pore_sites = interior(np.argwhere(lat == PORE), lat.shape)
        if len(ni_sites) == 0 or len(pore_sites) == 0:
            continue
        a = tuple(ni_sites[rng.integers(len(ni_sites))])
        b = tuple(pore_sites[rng.integers(len(pore_sites))])
        if manhattan(a, b) <= 1:
            continue
        nNa, _ = neigh_counts(lat, a)
        nNb, _ = neigh_counts(lat, b)
        o5v2 = 2.0 * (nNa - nNb)
        E0 = total_energy(lat, J)
        lat2 = lat.copy()
        lat2[a] = PORE
        lat2[b] = NI
        actual = total_energy(lat2, J) - E0
        max_err = max(max_err, abs(actual - o5v2))
    print(f"  max|actual - 2*(nb(a)-nb(b))| = {max_err:.3e}")

    # theta=180: YSZ must become indistinguishable from pore for Ni.
    print("\ntheta=180 equivalence (YSZ should act exactly like pore):")
    J180 = make_J(1.0, -1.0)
    diffs = 0
    for _ in range(500):
        lat = rng.integers(0, 3, size=(4, 4, 4)).astype(np.int8)
        latP = lat.copy()
        latP[latP == YSZ] = PORE
        ni_sites = interior(np.argwhere(lat == NI), lat.shape)
        pore_sites = interior(np.argwhere(lat == PORE), lat.shape)
        if len(ni_sites) == 0 or len(pore_sites) == 0:
            continue
        a = tuple(ni_sites[rng.integers(len(ni_sites))])
        b = tuple(pore_sites[rng.integers(len(pore_sites))])
        if manhattan(a, b) <= 1:
            continue
        d1 = predicted(lat, a, b, 1.0, -1.0)
        nNa, _ = neigh_counts(latP, a)
        nNb, _ = neigh_counts(latP, b)
        d2 = 2.0 * (nNa - nNb)
        if abs(d1 - d2) > 1e-9:
            diffs += 1
    print(f"  mismatches: {diffs} (expect 0)")

    # Physical scale: is kT_physical/J_NP large enough to cross a barrier?
    print("\nPhysical temperature scale:")
    gamma_Ni = 2.0                      # J/m^2
    for vox_nm in (17.9, 20.0, 29.14):
        A = (vox_nm * 1e-9) ** 2
        J_phys = gamma_Ni * A
        kT = 1.380649e-23 * (950 + 273.15)
        cos_t = -np.cos(np.deg2rad(50))
        smallest_barrier = (1.0 + cos_t)          # in units of J_NP
        ratio = kT / J_phys
        p = np.exp(-smallest_barrier / ratio) if ratio > 0 else 0.0
        print(f"  voxel {vox_nm:5.2f} nm: J_NP = {J_phys:.3e} J, "
              f"kT/J_NP = {ratio:.3e}, P(smallest barrier) = {p:.3e}")


if __name__ == "__main__":
    run()
