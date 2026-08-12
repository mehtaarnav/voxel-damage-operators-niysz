"""Sequential greedy vs batched greedy on the RECOVERED O5v2 structure.

The dA = 2*(nN(a) - nN(b)) identity is exact for ONE isolated move against a
current neighbour field. apply_o5v2b applies ~59 moves per round, ranked once
against a STALE field, so the per-move guarantees do not compose -- which is
why the corrected batched operator accepts ~96% of moves and still raises
S_spec.

The sequential operator -- one move, recompute the field, repeat -- is the
version the identity actually describes. It has never been run. If its
acceptance is zero, the recorded closure means what it claimed. If it accepts
and reduces area, the closure is wrong on its own terms.

Structure recovered from the original inline heredoc (not previously saved):
    60^3, two R=10 spheres at z=15 and z=45, 3x3 neck, YSZ slab at x>=55.
Verified by S_spec = 0.45052.
"""
import numpy as np
import scipy.ndimage as ndi

ST7 = ndi.generate_binary_structure(3, 1)       # centre included (committed)
ST6 = ST7.copy()
ST6[1, 1, 1] = False                             # centre excluded (frozen spec)


def original_structure():
    z, y, x = np.ogrid[:60, :60, :60]
    big = ((((z - 15) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100) |
           (((z - 45) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100))
    big = np.array(big)
    big[15:46, 29:32, 29:32] = True
    ysz = np.zeros_like(big)
    ysz[:, :, 55:] = True
    return big, ysz


def nbcount(ni, stencil=ST6):
    return ndi.convolve(ni.astype(np.int16), stencil.astype(np.int16),
                        mode="constant", cval=0)


def s_spec(ni):
    nb = nbcount(ni, ST6)
    return float((6 - nb)[ni].sum()) / max(int(ni.sum()), 1)


def neck_metrics(ni):
    """Free span between the spheres is z = 25..35 for the 3x3 column."""
    col = ni[25:36, 29:32, 29:32]
    per_slice = col.sum(axis=(1, 2))
    return int(col.sum()), int(per_slice.min())


def sequential_greedy(ni0, ysz, max_moves, tol_plateau=200):
    """One move at a time, neighbour field recomputed after every accept.

    Greedy: remove at the Ni surface site with fewest Ni neighbours, add at the
    pore-front site with most. Accept iff dA <= 0, i.e. nN(a) <= nN(b).
    """
    ni = ni0.copy()
    interior = np.zeros_like(ni)
    interior[1:-1, 1:-1, 1:-1] = True
    proposed = accepted = 0
    plateau = 0
    for _ in range(max_moves):
        nb = nbcount(ni, ST6).ravel()
        surf = ni & ~ndi.binary_erosion(ni, structure=ST7) & interior
        front = (~ni) & (~ysz) & ndi.binary_dilation(ni, structure=ST7) & interior
        si = np.flatnonzero(surf.ravel())
        fi = np.flatnonzero(front.ravel())
        if si.size == 0 or fi.size == 0:
            break
        a_order = si[np.argsort(nb[si], kind="stable")]
        b_order = fi[np.argsort(-nb[fi], kind="stable")]
        shp = ni.shape
        did = False
        # first admissible non-adjacent pair
        for ia in a_order[:50]:
            for ib in b_order[:50]:
                ia, ib = int(ia), int(ib)
                za, ya, xa = np.unravel_index(ia, shp)
                zb, yb, xb = np.unravel_index(ib, shp)
                if abs(za - zb) + abs(ya - yb) + abs(xa - xb) <= 1:
                    continue
                dA = 2 * (int(nb[ia]) - int(nb[ib]))
                proposed += 1
                if dA <= 0:
                    ni.ravel()[ia] = False
                    ni.ravel()[ib] = True
                    accepted += 1
                    plateau = plateau + 1 if dA == 0 else 0
                    did = True
                break
            if did:
                break
        if not did or plateau > tol_plateau:
            break
    return ni, proposed, accepted, plateau


def main():
    ni0, ysz = original_structure()
    s0 = s_spec(ni0)
    nvol0, nmin0 = neck_metrics(ni0)
    front0 = ((~ni0) & (~ysz) &
              ndi.binary_dilation(ni0, structure=ST7)).sum()
    print("RECOVERED O5v2 STRUCTURE")
    print(f"  S_spec   = {s0:.5f}   (original reported 0.45052)")
    print(f"  match    = {abs(s0 - 0.45052) < 5e-5}")
    print(f"  neck: volume={nvol0}  min cross-section={nmin0}")
    print(f"  pore-front sites = {int(front0)}   Ni voxels = {int(ni0.sum())}")

    # Confirm the stencil comparison the user reported.
    interior = np.zeros_like(ni0)
    interior[1:-1, 1:-1, 1:-1] = True
    surf = ni0 & ~ndi.binary_erosion(ni0, structure=ST7) & interior
    front = (~ni0) & (~ysz) & ndi.binary_dilation(ni0, structure=ST7) & interior
    for name, st in (("committed (7-element)", ST7), ("true 6-neighbour", ST6)):
        nb = nbcount(ni0, st).ravel()
        nba_min = int(nb[np.flatnonzero(surf.ravel())].min())
        nbb_max = int(nb[np.flatnonzero(front.ravel())].max())
        print(f"  {name:24s} nba_min={nba_min} nbb_max={nbb_max} "
              f"reject_at_t0={nba_min > nbb_max}")

    print("\nSEQUENTIAL GREEDY (field recomputed after every accepted move)")
    print(f"{'budget':>8} {'proposed':>9} {'accepted':>9} {'rate':>7} "
          f"{'S_spec':>9} {'dS':>9} {'neckvol':>8} {'neckmin':>8} {'dV':>5}")
    print("-" * 82)
    for budget in (61, 305, 1220):
        ni, prop, acc, plat = sequential_greedy(ni0, ysz, budget)
        s1 = s_spec(ni)
        nv, nm = neck_metrics(ni)
        print(f"{budget:>8} {prop:>9} {acc:>9} {acc/max(prop,1):>7.3f} "
              f"{s1:>9.5f} {s1-s0:>+9.5f} {nv:>8} {nm:>8} "
              f"{int(ni.sum())-int(ni0.sum()):>5}")

    print("\nGate (ii) of PREREG_O5V2_OPTIONB.md: S_spec(1) < S_spec(0), STRICT.")
    print("dS < 0 => sequential greedy PASSES where batched failed.")
    print("dS >= 0 => the closure holds for the sequential operator too.")


if __name__ == "__main__":
    main()
