"""Validate cmlib/seqgreedy.py against ground truth before any ROI run.

Checks, on the recovered O5v2 structure and on random blobs:
  1. the incrementally maintained nN field equals a fresh convolution
  2. volume is exactly conserved
  3. specific surface area is non-increasing (every accepted move has dA <= 0)
  4. YSZ is never overwritten
"""
import numpy as np
import scipy.ndimage as ndi

from cmlib.seqgreedy import SeqGreedy, _nb_field, ST7, ST6


def s_spec(ni):
    nb = ndi.convolve(ni.astype(np.int16), ST6.astype(np.int16),
                      mode="constant", cval=0)
    return float((6 - nb)[ni].sum()) / max(int(ni.sum()), 1)


def original_structure():
    z, y, x = np.ogrid[:60, :60, :60]
    big = ((((z - 15) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100) |
           (((z - 45) ** 2 + (y - 30) ** 2 + (x - 30) ** 2) < 100))
    big = np.array(big)
    big[15:46, 29:32, 29:32] = True
    ysz = np.zeros_like(big)
    ysz[:, :, 55:] = True
    return big, ysz


def random_blob(n=48, seed=0, phi=0.35, ysz_frac=0.25):
    rng = np.random.default_rng(seed)
    f = ndi.gaussian_filter(rng.random((n, n, n)), 3.0)
    ni = f > np.quantile(f, 1 - phi)
    g = ndi.gaussian_filter(rng.random((n, n, n)), 3.0)
    ysz = (g > np.quantile(g, 1 - ysz_frac)) & (~ni)
    return ni, ysz


def check(name, ni0, ysz, n_moves):
    v0 = int(ni0.sum())
    s0 = s_spec(ni0)
    op = SeqGreedy(ni0, ysz, seed=300)
    op.run(n_moves)
    ni = op.ni

    fresh = _nb_field(ni)
    nb_ok = bool(np.array_equal(fresh, op.nb))
    v1 = int(ni.sum())
    s1 = s_spec(ni)
    ysz_ok = not bool((ni & ysz).any())

    print(f"{name:<34} moves={op.accepted:>6}/{op.proposed:<6} "
          f"neutral={op.neutral:>6}  dV={v1-v0:>3}  "
          f"S {s0:.5f}->{s1:.5f} ({s1-s0:+.5f})")
    print(f"{'':34} nb field exact: {nb_ok}   "
          f"volume conserved: {v1 == v0}   "
          f"area non-increasing: {s1 <= s0 + 1e-12}   "
          f"YSZ intact: {ysz_ok}")
    return nb_ok and (v1 == v0) and (s1 <= s0 + 1e-12) and ysz_ok


def main():
    ok = True
    ni0, ysz = original_structure()
    for m in (61, 305, 1220):
        ok &= check(f"recovered O5v2 struct, {m} moves", ni0, ysz, m)
    for seed in (0, 1, 2):
        ni0, ysz = random_blob(seed=seed)
        ok &= check(f"random blob seed {seed}, 2000 moves", ni0, ysz, 2000)
    print("\nALL CHECKS PASS" if ok else "\nVALIDATION FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
