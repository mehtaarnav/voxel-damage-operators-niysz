"""Unit tests for the Step 2 operators. Run before any scientific run."""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy import ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.damage import apply_d4  # noqa: E402
from cmlib.damage2 import (  # noqa: E402
    apply_o1, apply_o2, apply_o3, survival_mask, tpb_density_um2,
)
from cmlib.project2 import rasterize_ysz  # noqa: E402

S6 = ndi.generate_binary_structure(3, 1)
fails = []


def check(name, cond, detail=""):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name} {detail}")
    if not cond:
        fails.append(name)


print("== survival_mask ==")
n = 20000
s0 = survival_mask(n, 0, 0.25, 1)
s1 = survival_mask(n, 1, 0.25, 1)
s5 = survival_mask(n, 5, 0.25, 1)
check("n=0 -> all survive", s0.all())
check("monotone in n", bool((s5 <= s1).all() and (s1 <= s0).all()))
check("marginal ~ (1-p)^n", abs(s5.mean() - 0.75 ** 5) < 0.01,
      f"{s5.mean():.4f} vs {0.75**5:.4f}")

print("== O1 is bit-identical to the frozen apply_d4 ==")
rng = np.random.default_rng(0)
ni = np.zeros((40, 40, 40), bool)
ni[8:32, 8:32, 8:32] = True
ysz = np.zeros_like(ni)
ysz[0:4] = True
a, _ = apply_o1(ni, ysz, 3, 7)
b, _ = apply_d4(ni, ysz, 3, 0.35, 1, 7)
check("O1 == apply_d4", bool(np.array_equal(a, b)))
z, _ = apply_o1(ni, ysz, 0, 7)
check("O1 n=0 is a no-op", bool(np.array_equal(z, ni)))

print("== O1 leaves YSZ untouched / removes only Ni ==")
y2 = ysz.copy()
a, _ = apply_o1(ni, ysz, 2, 3)
check("YSZ array unchanged", bool(np.array_equal(ysz, y2)))
check("no Ni created inside YSZ", bool(not (a & ysz).any()))

print("== O2 candidate selection at the frozen 25th percentile ==")
from cmlib.damage2 import o2_candidates  # noqa: E402
cand_t, thr_t = o2_candidates(None, np.arange(1.0, 101.0), 25)
check("25th pct selects ~25% of throats", int(cand_t.sum()) == 25,
      f"(n={int(cand_t.sum())}, thr={thr_t:.2f})")

print("== O2 severs a narrow throat and disconnects ==")
# three cubes in a line, two throats of different width; the narrow one is
# the only candidate at pct=50, so the sever is deterministic in identity.
m = np.zeros((46, 30, 30), bool)
m[3:11, 10:20, 10:20] = True
m[15:23, 10:20, 10:20] = True
m[27:35, 10:20, 10:20] = True
m[10:16, 14:16, 14:16] = True      # narrow throat 0-1
m[22:28, 12:18, 12:18] = True      # wide   throat 1-2
regions = np.zeros(m.shape, np.int32)
regions[m] = 2
regions[:13][m[:13]] = 1
regions[25:][m[25:]] = 3
conns = np.array([[0, 1], [1, 2]])
diam = np.array([40.0, 120.0])
check("chain starts connected", ndi.label(m, structure=S6)[1] == 1)
out, d = apply_o2(m, regions, conns, diam, n_rounds=40, seed=5, pct=50)
check("O2 severs the narrow throat", d["n_severed"] == 1,
      f"(cand={d['n_candidates']}, thr={d['neck_threshold_nm']:.0f})")
check("O2 keeps one component only", ndi.label(out, structure=S6)[1] == 1)
check("O2 removed volume", out.sum() < m.sum(), f"{m.sum()} -> {out.sum()}")
out0, d0 = apply_o2(m, regions, conns, diam, n_rounds=0, seed=5, pct=50)
check("O2 n=0 is a no-op", bool(np.array_equal(out0, m)))

print("== O2 monotonicity in n_rounds ==")
prev = None
mono = True
for nr in (0, 1, 2, 4, 8, 16):
    o, dd = apply_o2(m, regions, conns, diam, nr, seed=11, pct=50)
    if prev is not None and dd["n_severed"] < prev:
        mono = False
    prev = dd["n_severed"]
check("severed count non-decreasing", mono)

print("== O3 removes neck volume and updates the graph ==")
c = {0: (15., 15., 15.), 1: (15., 15., 35.), 2: (15., 35., 15.)}
pairs = [(0, 1), (0, 2)]
sint = np.array([True, True])
shape = (50, 50, 50)
nimask = np.zeros(shape, bool)
y_full = rasterize_ysz(c, pairs, sint, 6.0, 5.0, shape)
y0, intact0, i0 = apply_o3(c, pairs, sint, 6.0, 5.0, shape, nimask, 0, 3)
check("O3 n=0 preserves all contacts", i0["n_intact_post"] == 2)
check("O3 n=0 mask == pristine", bool(np.array_equal(y0, y_full)))
yb, intactb, ib = apply_o3(c, pairs, sint, 6.0, 5.0, shape, nimask, 60, 3)
check("O3 fractures at high n", ib["n_fractured"] == 2,
      f"(fractured={ib['n_fractured']})")
check("O3 removes neck volume", yb.sum() < y_full.sum(),
      f"{y_full.sum()} -> {yb.sum()}")
check("O3 returns updated graph", bool((~intactb).all()))
check("O3 fully fractured -> 3 grains", ndi.label(yb, structure=S6)[1] == 3)
check("O3 leaves Ni untouched", bool(not (yb & nimask).any()))

print("== O3 monotonicity ==")
prev, mono = None, True
for nr in (0, 1, 3, 6, 12, 30):
    _, _, dd = apply_o3(c, pairs, sint, 6.0, 5.0, shape, nimask, nr, 21)
    if prev is not None and dd["n_fractured"] < prev:
        mono = False
    prev = dd["n_fractured"]
check("fractured count non-decreasing", mono)

print("== TPB units ==")
# analytic: a single straight triple line of length L in a WxWxW box.
# Build three quadrants meeting along the z axis -> one triple line, length W.
W = 40
zz, yy, xx = np.ogrid[:W, :W, :W]
a_ = (yy >= W // 2) & (xx >= W // 2)
b_ = (yy < W // 2) & (xx >= W // 2)
niT = np.broadcast_to(a_, (W, W, W)).copy()
yszT = np.broadcast_to(b_, (W, W, W)).copy()
vox = 20.0
t = tpb_density_um2(niT, yszT, vox)
expect = (W * vox) / ((W * vox) ** 3) * 1e6
check("single triple line matches analytic", abs(t - expect) / expect < 0.05,
      f"{t:.4f} vs {expect:.4f} um^-2")
check("TPB zero when a phase is absent",
      tpb_density_um2(niT, np.zeros_like(yszT), vox) == 0.0)


print("== O5 volume conservation, monotone volume, YSZ untouched ==")
from cmlib.damage2 import apply_o5  # noqa: E402
big = np.zeros((60, 60, 60), bool)
zz2, yy2, xx2 = np.ogrid[:60, :60, :60]
big |= ((zz2 - 15) ** 2 + (yy2 - 30) ** 2 + (xx2 - 30) ** 2) < 100
big |= ((zz2 - 45) ** 2 + (yy2 - 30) ** 2 + (xx2 - 30) ** 2) < 100
big[15:46, 29:32, 29:32] = True   # thin neck; free span is z=25..35
yszB = np.zeros_like(big)
yszB[:, :, 55:] = True
pre_neck = int(big[27:34, 29:32, 29:32].sum())
o5, d5 = apply_o5(big, yszB, 6, 4)
check("O5 conserves volume <=0.5%", d5["volume_error"] <= 0.005,
      f"({d5['voxels_pre']} -> {d5['voxels_post']}, err={d5['volume_error']:.5f})")
check("O5 actually moved voxels", d5["voxels_moved"] > 0,
      f"({d5['voxels_moved']})")
check("O5 leaves YSZ untouched", bool(not (o5 & yszB).any()))
post_neck = int(o5[27:34, 29:32, 29:32].sum())
check("O5 thins the neck", post_neck < pre_neck, f"{pre_neck} -> {post_neck}")
o5a, da = apply_o5(big, yszB, 0, 4)
check("O5 n=0 is a no-op", bool(np.array_equal(o5a, big)))
errs = []
for nr in (1, 3, 6, 10):
    _, dd = apply_o5(big, yszB, nr, 9)
    errs.append(dd["volume_error"])
check("O5 volume error stays <=0.5% at all n", max(errs) <= 0.005,
      f"max={max(errs):.5f}")

print("\n" + ("ALL TESTS PASS" if not fails else f"FAILURES: {fails}"))
sys.exit(1 if fails else 0)
