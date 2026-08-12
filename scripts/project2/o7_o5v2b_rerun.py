"""Re-run O5v2 Option B with the corrected stencil.

O5v2 Option B reported acceptance rate EXACTLY 0.000 at every intensity and
concluded that the agglomeration route is CLOSED, because "spheres joined by a
straight cylinder are already at a local minimum with respect to single-voxel
swaps".

PREREG_O5V2_OPTIONB.md freezes the algorithm as:
    dA = 2*(nb(a) - nb(b)),  nb = "the count of Ni 6-neighbours",
    accept iff dA <= 0  <=>  nb(a) <= nb(b).

cmlib/damage2.py:31 sets STRUCT6 = generate_binary_structure(3,1), which has
sum 7 -- the CENTRE IS INCLUDED. Convolving the Ni mask with it yields

    nb(a) = nN(a) + 1   for a Ni site a      (centre contributes 1)
    nb(b) = nN(b)       for a pore site b    (centre contributes 0)

so the implemented predicate is nN(a) + 1 <= nN(b), i.e. nN(a) < nN(b) --
STRICTLY STRONGER than the frozen spec's nN(a) <= nN(b). Every area-NEUTRAL
(dA = 0) move is silently rejected.

This script runs the committed operator and a corrected copy (identical except
for the stencil) on the same structures and compares acceptance.

The original O5v2 driver is not in the repo -- only a transcription of its
results (scripts/project2/o5v2_transcribe.py). Exact re-execution is therefore
impossible; instead both operators are run on a FAMILY of structures, including
the project's own lattice generator, so the comparison does not depend on
matching one geometry.
"""
import numpy as np
import scipy.ndimage as ndi

from cmlib.damage2 import apply_o5v2b, O5V2_P_COARSEN
from cmlib import synth

ST7 = ndi.generate_binary_structure(3, 1)      # centre included, sum 7
ST6 = ST7.copy()
ST6[1, 1, 1] = False                            # centre excluded, sum 6


def apply_o5v2b_fixed(ni_mask, ysz_mask, n_rounds, seed,
                      p_coarsen=O5V2_P_COARSEN):
    """Byte-for-byte apply_o5v2b except nb uses the centre-EXCLUDED stencil,
    which is what PREREG_O5V2_OPTIONB.md specifies."""
    ni = ni_mask.copy()
    rng = np.random.default_rng(seed)
    n_surf0 = int((ni & ~ndi.binary_erosion(ni, structure=ST7)).sum())
    k_move = int(round(p_coarsen * n_surf0))
    proposed = accepted = 0
    for _ in range(max(0, int(n_rounds))):
        if k_move <= 0:
            break
        nb = ndi.convolve(ni.astype(np.int16), ST6.astype(np.int16),
                          mode="constant", cval=0)          # <-- the fix
        surf = ni & ~ndi.binary_erosion(ni, structure=ST7)
        front = (~ni) & (~ysz_mask) & ndi.binary_dilation(ni, structure=ST7)
        si = np.flatnonzero(surf.ravel())
        fi = np.flatnonzero(front.ravel())
        if si.size == 0 or fi.size == 0:
            break
        nba, nbb = nb.ravel()[si], nb.ravel()[fi]
        a_ord = si[np.argsort(nba, kind="stable")]
        b_ord = fi[np.argsort(-nbb, kind="stable")]
        nba_s = np.sort(nba, kind="stable")
        nbb_s = -np.sort(-nbb, kind="stable")
        k = min(k_move, a_ord.size, b_ord.size)
        flat = ni.ravel()
        shp = ni.shape
        moved_this = 0
        for t in range(k):
            proposed += 1
            if nba_s[t] > nbb_s[t]:
                break
            ia, ib = int(a_ord[t]), int(b_ord[t])
            za, ya, xa = np.unravel_index(ia, shp)
            zb, yb, xb = np.unravel_index(ib, shp)
            if abs(int(za) - int(zb)) + abs(int(ya) - int(yb)) + \
               abs(int(xa) - int(xb)) <= 1:
                continue
            flat[ia] = False
            flat[ib] = True
            accepted += 1
            moved_this += 1
        if moved_this == 0:
            break
    n0, n1 = int(ni_mask.sum()), int(ni.sum())
    return ni, dict(voxels_pre=n0, voxels_post=n1,
                    volume_error=abs(n1 - n0) / max(n0, 1),
                    proposed=proposed, accepted=accepted,
                    acceptance_rate=accepted / max(proposed, 1))


def s_spec(ni):
    """Specific surface area: exposed Ni faces per Ni voxel."""
    nb = ndi.convolve(ni.astype(np.int16), ST6.astype(np.int16),
                      mode="constant", cval=0)
    v = int(ni.sum())
    return float((6 - nb)[ni].sum()) / max(v, 1)


def dumbbell(n=64, R=13.0, r_neck=4.0, sep=26.0):
    z, y, x = np.ogrid[:n, :n, :n]
    c = n / 2.0
    s1 = (z - c) ** 2 + (y - c) ** 2 + (x - (c - sep / 2)) ** 2 <= R ** 2
    s2 = (z - c) ** 2 + (y - c) ** 2 + (x - (c + sep / 2)) ** 2 <= R ** 2
    cyl = (((z - c) ** 2 + (y - c) ** 2) <= r_neck ** 2) & \
          (np.abs(x - c) <= sep / 2)
    return s1 | s2 | cyl


def project_lattice(seed=300):
    """The project's own generator -- closest available proxy to the O5v2
    synthetic structure."""
    rng = np.random.default_rng(seed)
    centres, pairs, shape = synth.jittered_lattice_geometry(
        nlat=4, pitch_vox=22, r_vox=8.0, margin=4, jitter_frac=0.12, rng=rng)
    ni = np.zeros(shape, dtype=bool)
    zz, yy, xx = np.indices(shape)
    for (cz, cy, cx) in centres.values():
        ni |= ((zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2) <= 8.0 ** 2
    for a, b in pairs:                      # straight necks, as in O5v2
        z0, y0, x0 = centres[a]
        z1, y1, x1 = centres[b]
        for t in np.linspace(0, 1, 60):
            pz, py, px = z0 + t * (z1 - z0), y0 + t * (y1 - y0), \
                x0 + t * (x1 - x0)
            ni |= ((zz - pz) ** 2 + (yy - py) ** 2 +
                   (xx - px) ** 2) <= 3.0 ** 2
    return ni


def main():
    cases = [
        ("dumbbell R13 neck4", dumbbell()),
        ("dumbbell R13 neck3", dumbbell(r_neck=3.0)),
        ("dumbbell R16 neck5", dumbbell(R=16.0, r_neck=5.0, sep=32.0)),
        ("project lattice s300", project_lattice(300)),
        ("project lattice s301", project_lattice(301)),
        ("project lattice s302", project_lattice(302)),
    ]

    print("O5v2 Option B: committed (centre-included nb) vs corrected "
          "(centre-excluded nb)")
    print("Frozen spec: nb = count of Ni 6-neighbours. "
          "Committed code includes the centre.\n")
    hdr = (f"{'structure':<24} {'n':>2} "
           f"{'acc_committed':>14} {'acc_corrected':>14} "
           f"{'S_spec0':>9} {'S_spec_fix':>11} {'dV_fix':>7}")
    print(hdr)
    print("-" * len(hdr))

    for name, ni in cases:
        ysz = np.zeros_like(ni)
        s0 = s_spec(ni)
        for n_rounds in (1, 3, 5):
            _, ic = apply_o5v2b(ni, ysz, n_rounds, seed=300)
            nf, if_ = apply_o5v2b_fixed(ni, ysz, n_rounds, seed=300)
            print(f"{name:<24} {n_rounds:>2} "
                  f"{ic['acceptance_rate']:>14.4f} "
                  f"{if_['acceptance_rate']:>14.4f} "
                  f"{s0:>9.5f} {s_spec(nf):>11.5f} "
                  f"{if_['voxels_post'] - if_['voxels_pre']:>7d}")
        print()

    print("Gate (ii) of PREREG_O5V2_OPTIONB.md requires S_spec(1) < S_spec(0),")
    print("STRICT. That gate is what closed the agglomeration route.")


if __name__ == "__main__":
    main()
