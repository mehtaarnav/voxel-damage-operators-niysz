"""R_Ni sanity check + O6 gate on the pruning-invariant metric.
Frozen in PREREG_RNI_METRIC.md (bfc5d89). p_erode not adjusted."""
import os, sys, time
import numpy as np, pandas as pd, tifffile
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cmlib.damage2 import apply_o6, tpb_density_um2
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE
from cmlib.io import label_histogram, slice_paths
from cmlib.percolation import percolation_summary, percolating_mask
from cmlib.phases import assign_labels
from cmlib.roi import tile_rois
OUT = os.path.join(ROOT, "out", "project2"); AXIS, CONN = 2, 6

def load(folder, mapping, r):
    ps = slice_paths(folder)
    sh = (r["z1"]-r["z0"], r["y1"]-r["y0"], r["x1"]-r["x0"])
    ni = np.empty(sh, bool); ysz = np.empty(sh, bool)
    for i, z in enumerate(range(r["z0"], r["z1"])):
        a = tifffile.imread(ps[z])[r["y0"]:r["y1"], r["x0"]:r["x1"]]
        ni[i] = a == mapping["Ni"]; ysz[i] = a == mapping["YSZ"]
    return ni, ysz

def r_ni(mask, pristine_n):
    """Spanning-cluster Ni voxels / PRISTINE Ni voxels."""
    return int(percolating_mask(mask, axis=AXIS, connectivity=CONN).sum()) / pristine_n

rows = []
sample = {s[2]: s for s in SAMPLES if s[3] == "pristine"}
for grain in ("fine", "medium", "coarse"):
    k = sample[grain]
    nz, ny, nx_ = k[6], k[5], k[4]; vz, vy, vx = k[9], k[8], k[7]
    mapping = assign_labels(label_histogram(k[1])["counts"], ZENODO_LABEL_NOTE[k[0]])
    side = 12.0 if grain == "coarse" else 8.0
    r = tile_rois(nz, ny, nx_, vz, vy, vx, side, max_rois=1)[0]
    ni, ysz = load(k[1], mapping, r)
    n0 = int(ni.sum()); vox = float((vz*vy*vx)**(1/3))
    p0 = percolation_summary(ni, axis=AXIS, connectivity=CONN, check_other_axes=False)["P_span"]
    r0 = r_ni(ni, n0); t0v = tpb_density_um2(ni, ysz, vox)
    print(f"{grain:7s} pristine P_span={p0:.6f}  R_Ni(0)={r0:.6f}  "
          f"|diff|={abs(r0-p0):.2e}  {'SANITY OK' if abs(r0-p0)<1e-9 else '** SANITY FAIL **'}  "
          f"TPB={t0v:.3f}", flush=True)
    for n in (1, 3, 5, 8):
        t = time.time(); d, _i = apply_o6(ni, ysz, n, 300)
        rows.append(dict(anode=grain, n_rounds=n, pristine_P_span=p0, R_Ni_0=r0,
                         R_Ni=r_ni(d, n0), vol_loss=1-int(d.sum())/n0,
                         tpb_um2=tpb_density_um2(d, ysz, vox), tpb_pristine=t0v,
                         seconds=round(time.time()-t, 1)))
        print(f"   n={n}: R_Ni={rows[-1]['R_Ni']:.4f} volLoss={rows[-1]['vol_loss']:.4f} "
              f"TPB={rows[-1]['tpb_um2']:.2f} [{rows[-1]['seconds']}s]", flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(OUT, "c1real_rni_gate.csv"), index=False)
    del ni, ysz
df = pd.DataFrame(rows); ok = True
print("\nR_Ni GATE")
for a in df.anode.unique():
    s = df[df.anode == a].sort_values("n_rounds")
    mono = bool((s.R_Ni.diff().dropna() <= 1e-12).all())
    below = bool(s.R_Ni.iloc[0] <= s.R_Ni_0.iloc[0] + 1e-12)
    print(f"  {a:7s} R_Ni monotone non-increasing={mono}  R_Ni(1)<=R_Ni(0)={below}  "
          f"R_Ni: {s.R_Ni_0.iloc[0]:.4f} -> {list(s.R_Ni.round(4))}")
    ok = ok and mono and below
print("\n  R_Ni GATE: " + ("PASS" if ok else "** FAIL **"))
