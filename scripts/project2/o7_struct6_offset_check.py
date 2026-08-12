"""Does the centre-included STRUCT6 explain O5v2 Option B's zero acceptance?

cmlib/damage2.py:31 sets STRUCT6 = generate_binary_structure(3,1), which has
sum 7 -- it INCLUDES the centre voxel. Convolving the Ni mask with it gives

    nb(site) = (1 if site is Ni else 0) + (# Ni 6-neighbours)

so for a Ni source a and a pore sink b:  nb(a) = nN(a) + 1,  nb(b) = nN(b).

The exact identity (verified in o7_derivation_check.py) is
    dA = 2*(nN(a) - nN(b)),  so dA <= 0  <=>  nN(a) <= nN(b).
The code tests nb(a) <= nb(b), i.e. nN(a) + 1 <= nN(b), i.e. nN(a) < nN(b).

That is strictly stronger: it rejects every AREA-NEUTRAL (dA = 0) move.

Question: was O5v2's "acceptance exactly 0.000" a property of the geometry, or
an artifact of the offset? Compare the two acceptance predicates on the same
structure.
"""
import numpy as np
import scipy.ndimage as ndi

ST7 = ndi.generate_binary_structure(3, 1)          # centre INCLUDED, sum 7
ST6 = ST7.copy()
ST6[1, 1, 1] = False                                # centre EXCLUDED, sum 6


def make_dumbbell(n=64, R=13.0, r_neck=4.0, sep=26.0):
    z, y, x = np.ogrid[:n, :n, :n]
    c = n / 2.0
    s1 = (z - c) ** 2 + (y - c) ** 2 + (x - (c - sep / 2)) ** 2 <= R ** 2
    s2 = (z - c) ** 2 + (y - c) ** 2 + (x - (c + sep / 2)) ** 2 <= R ** 2
    cyl = (((z - c) ** 2 + (y - c) ** 2) <= r_neck ** 2) & \
          (np.abs(x - c) <= sep / 2)
    return s1 | s2 | cyl


def counts(ni, stencil):
    return ndi.convolve(ni.astype(np.int16), stencil.astype(np.int16),
                        mode="constant", cval=0)


def report(ni, label):
    interior = np.zeros_like(ni)
    interior[1:-1, 1:-1, 1:-1] = True
    surf = ni & ~ndi.binary_erosion(ni, structure=ST7) & interior
    front = (~ni) & ndi.binary_dilation(ni, structure=ST7) & interior

    nb7 = counts(ni, ST7).ravel()      # as the committed code computes it
    nb6 = counts(ni, ST6).ravel()      # true Ni-neighbour count
    si = np.flatnonzero(surf.ravel())
    fi = np.flatnonzero(front.ravel())

    nN_a, nN_b = nb6[si], nb6[fi]
    code_a, code_b = nb7[si], nb7[fi]

    print(f"\n--- {label} ---")
    print(f"surface Ni sites: {si.size},  pore-front sites: {fi.size}")
    print(f"true nN(a): min={nN_a.min()} max={nN_a.max()}")
    print(f"true nN(b): min={nN_b.min()} max={nN_b.max()}")

    # Committed predicate: nb7(a) <= nb7(b)
    code_possible = code_a.min() <= code_b.max()
    # Correct predicate: nN(a) <= nN(b)
    true_possible = nN_a.min() <= nN_b.max()
    # Strictly area-lowering only
    strict_possible = nN_a.min() < nN_b.max()

    print(f"committed predicate nb7(a) <= nb7(b) satisfiable: {code_possible}"
          f"   (min {code_a.min()} vs max {code_b.max()})")
    print(f"correct  predicate  nN(a) <= nN(b) satisfiable: {true_possible}"
          f"   (min {nN_a.min()} vs max {nN_b.max()})")
    print(f"strictly area-LOWERING nN(a) < nN(b) satisfiable: {strict_possible}")

    # How many area-neutral pairs exist that the committed code discards?
    neutral = 0
    for v in range(0, 7):
        na = int((nN_a == v).sum())
        nbv = int((nN_b == v).sum())
        if na and nbv:
            neutral += min(na, nbv)
    print(f"area-NEUTRAL (dA = 0) pair capacity discarded by the offset: "
          f"{neutral}")


def main():
    ni = make_dumbbell()
    print("O5v2 Option B geometry: two spheres joined by a straight cylinder.")
    print("O5v2 recorded acceptance = 0.000 and attributed it to geometry")
    print("('nba_min > nbb_max'). Testing that attribution.")
    report(ni, "dumbbell (O5v2 test geometry)")

    rng = np.random.default_rng(3)
    noisy = ni.copy()
    flip = rng.random(ni.shape) < 0.02
    noisy = noisy ^ (flip & ndi.binary_dilation(ni, ST7))
    report(noisy, "same geometry, 2% surface noise (more realistic)")


if __name__ == "__main__":
    main()
