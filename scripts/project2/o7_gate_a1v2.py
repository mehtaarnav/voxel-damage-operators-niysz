"""Gate A1v2, all five conditions, sequential greedy operator.

PREREG_O5V2_OPTIONB.md, Gate A1v2:
  (i)   |dPhi_Ni| <= 0.005
  (ii)  S_spec(1) < S_spec(0), strict and monotonic
  (iii) TPB(n) <= TPB(0)
  (iv)  R_Ni non-increasing
  (v)   YSZ untouched
  n = 1, 3, 5; seeds 300/301/302.

Run on the RECOVERED original structure (S_spec = 0.45052), i.e. like-for-like
with the run that closed the route.

NOTE on seeds: the committed apply_o5v2b constructs
`rng = np.random.default_rng(seed)` and then never uses it -- the operator is
deterministic despite taking a seed. The sequential operator here uses the seed
for tie-breaking among equal-nb candidates, which is the natural stochastic
element and makes the three seeds meaningful.
"""
import numpy as np
import scipy.ndimage as ndi

from cmlib.damage2 import tpb_density_um2

ST7 = ndi.generate_binary_structure(3, 1)
ST6 = ST7.copy()
ST6[1, 1, 1] = False


def original_structure():
    z, y, x = np.ogrid[:60, :60, :60]
    big = ((((z - 15) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100) |
           (((z - 45) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100))
    big = np.array(big)
    big[15:46, 29:32, 29:32] = True
    ysz = np.zeros_like(big)
    ysz[:, :, 55:] = True
    return big, ysz


def nbcount(ni):
    return ndi.convolve(ni.astype(np.int16), ST6.astype(np.int16),
                        mode="constant", cval=0)


def s_spec(ni):
    return float((6 - nbcount(ni))[ni].sum()) / max(int(ni.sum()), 1)


def phi_ni(ni):
    return float(ni.sum()) / ni.size


def r_ni_facespan(ni, pristine_total, axis=0):
    """Spanning-cluster Ni voxels / PRISTINE Ni voxels (PREREG_RNI_METRIC.md)."""
    lab, n = ndi.label(ni, structure=ST7)
    if n == 0:
        return 0.0
    lo = set(np.unique(lab.take(0, axis=axis))) - {0}
    hi = set(np.unique(lab.take(-1, axis=axis))) - {0}
    span = lo & hi
    if not span:
        return 0.0
    return float(np.isin(lab, list(span)).sum()) / max(pristine_total, 1)


def largest_component_frac(ni, pristine_total):
    lab, n = ndi.label(ni, structure=ST7)
    if n == 0:
        return 0.0
    sizes = np.bincount(lab.ravel())[1:]
    return float(sizes.max()) / max(pristine_total, 1)


def sequential_greedy(ni0, ysz, max_moves, seed, tol_plateau=200):
    ni = ni0.copy()
    rng = np.random.default_rng(seed)
    interior = np.zeros_like(ni)
    interior[1:-1, 1:-1, 1:-1] = True
    proposed = accepted = plateau = 0
    for _ in range(max_moves):
        nb = nbcount(ni).ravel()
        surf = ni & ~ndi.binary_erosion(ni, structure=ST7) & interior
        front = (~ni) & (~ysz) & ndi.binary_dilation(ni, structure=ST7) & interior
        si = np.flatnonzero(surf.ravel())
        fi = np.flatnonzero(front.ravel())
        if si.size == 0 or fi.size == 0:
            break
        # seed-dependent tie-breaking among equal nb
        ja = rng.random(si.size) * 1e-6
        jb = rng.random(fi.size) * 1e-6
        a_order = si[np.argsort(nb[si] + ja, kind="stable")]
        b_order = fi[np.argsort(-(nb[fi] + jb), kind="stable")]
        shp = ni.shape
        did = False
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
    return ni, proposed, accepted


def main():
    ni0, ysz0 = original_structure()
    tot0 = int(ni0.sum())
    s0 = s_spec(ni0)
    phi0 = phi_ni(ni0)
    tpb0 = tpb_density_um2(ni0, ysz0, voxel_nm=20.0)
    rni0 = r_ni_facespan(ni0, tot0)
    lcf0 = largest_component_frac(ni0, tot0)

    print("PRISTINE")
    print(f"  S_spec={s0:.5f}  Phi_Ni={phi0:.6f}  TPB={tpb0:.6f} um^-2")
    print(f"  R_Ni(face-span,z)={rni0:.4f}   largest-component frac={lcf0:.4f}")
    print(f"  Ni touches z=0 plane: {bool(ni0[0].any())}, "
          f"z=-1 plane: {bool(ni0[-1].any())}")
    ni_ysz_contact = int((ni0 & ndi.binary_dilation(ysz0, ST7)).sum())
    print(f"  Ni voxels adjacent to YSZ: {ni_ysz_contact}")

    # moves-per-round calibrated to the committed operator: k = 0.03 * n_surf0
    n_surf0 = int((ni0 & ~ndi.binary_erosion(ni0, structure=ST7)).sum())
    k_round = int(round(0.03 * n_surf0))
    print(f"\n  surface voxels={n_surf0}  moves/round k={k_round}\n")

    hdr = (f"{'seed':>5} {'n':>2} {'prop':>6} {'acc':>6} {'rate':>6} "
           f"{'S_spec':>9} {'dPhi':>9} {'TPB':>8} {'R_Ni':>7} {'LCF':>7} "
           f"{'YSZ ok':>7}")
    print(hdr)
    print("-" * len(hdr))

    results = {}
    for seed in (300, 301, 302):
        for n in (1, 3, 5):
            ni, prop, acc = sequential_greedy(ni0, ysz0, k_round * n, seed)
            s1 = s_spec(ni)
            dphi = abs(phi_ni(ni) - phi0)
            tpb = tpb_density_um2(ni, ysz0, voxel_nm=20.0)
            rni = r_ni_facespan(ni, tot0)
            lcf = largest_component_frac(ni, tot0)
            ysz_ok = not bool((ni & ysz0).any())
            results[(seed, n)] = (s1, dphi, tpb, rni, lcf, ysz_ok)
            print(f"{seed:>5} {n:>2} {prop:>6} {acc:>6} {acc/max(prop,1):>6.3f} "
                  f"{s1:>9.5f} {dphi:>9.2e} {tpb:>8.4f} {rni:>7.4f} "
                  f"{lcf:>7.4f} {str(ysz_ok):>7}")

    print("\nGATE A1v2 VERDICT")
    ok_i = all(r[1] <= 0.005 for r in results.values())
    print(f"  (i)   |dPhi_Ni| <= 0.005                     : "
          f"{'PASS' if ok_i else 'FAIL'}")
    strict = all(results[(s, 1)][0] < s0 for s in (300, 301, 302))
    mono = all(results[(s, 1)][0] >= results[(s, 3)][0] >= results[(s, 5)][0]
               for s in (300, 301, 302))
    print(f"  (ii)  S_spec(1) < S_spec(0), strict+monotonic: "
          f"{'PASS' if strict and mono else 'FAIL'}"
          f"   (strict={strict}, monotonic={mono})")
    ok_iii = all(r[2] <= tpb0 + 1e-12 for r in results.values())
    print(f"  (iii) TPB(n) <= TPB(0)                       : "
          f"{'PASS' if ok_iii else 'FAIL'}"
          f"   [TPB(0)={tpb0:.6f}]")
    ok_iv = all(results[(s, 1)][3] <= rni0 + 1e-12 and
                results[(s, 3)][3] <= results[(s, 1)][3] + 1e-12 and
                results[(s, 5)][3] <= results[(s, 3)][3] + 1e-12
                for s in (300, 301, 302))
    print(f"  (iv)  R_Ni non-increasing                    : "
          f"{'PASS' if ok_iv else 'FAIL'}   [R_Ni(0)={rni0:.4f}]")
    ok_v = all(r[5] for r in results.values())
    print(f"  (v)   YSZ untouched                          : "
          f"{'PASS' if ok_v else 'FAIL'}")

    print("\nDEGENERACY AUDIT — is each condition actually informative here?")
    print(f"  (iii) TPB(0) = {tpb0:.6f}. Ni-YSZ adjacent voxels = "
          f"{ni_ysz_contact}.")
    print(f"  (iv)  R_Ni(0) = {rni0:.4f}. Ni touches neither z face "
          f"=> face-spanning is identically 0.")


if __name__ == "__main__":
    main()
