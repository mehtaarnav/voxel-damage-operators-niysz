"""Relook, step 1: per-sub-block measurements on all six stacks.

Questions it answers, with no damage model involved:
  (a) Is the pristine->post change in Ni spanning fraction and TPB larger than
      the spatial scatter between sub-blocks of one stack?
  (b) Does the post-redox structure look like a dilated version of the
      pristine one?  Dilation predicts: solid fractions scale together, and the
      share of Ni/YSZ contact that touches pore goes up.

Per 5 um cube: phase fractions, Ni P_span (6-conn, mean over 3 axes),
TPB density (cmlib.tpb voxel-edge convention), Ni/YSZ face density, and the
fraction of Ni/YSZ faces where at least one of the two voxels has a pore
6-neighbour ("open" contact) versus none ("sealed" contact).
"""
import os, sys, time
import numpy as np, pandas as pd, tifffile
from scipy import ndimage as ndi
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE
from cmlib.io import label_histogram, slice_paths
from cmlib.phases import assign_labels
from cmlib.roi import tile_rois
from cmlib.tpb import tpb_density_volume

OUT = os.path.join(ROOT, "out", "relook"); os.makedirs(OUT, exist_ok=True)
SIDE_UM = 5.0
S6 = ndi.generate_binary_structure(3, 1)


def pspan_mean(m):
    lab, n = ndi.label(m, S6)
    if n == 0 or m.sum() == 0:
        return 0.0, [0.0, 0.0, 0.0]
    sizes = np.bincount(lab.ravel())
    per = []
    for ax in range(3):
        lo = np.unique(np.take(lab, 0, axis=ax)); hi = np.unique(np.take(lab, -1, axis=ax))
        sp = np.intersect1d(lo, hi); sp = sp[sp > 0]
        per.append(sizes[sp].sum() / m.sum())
    return float(np.mean(per)), per


def faces(a, b, ax):
    """Pairs of 6-neighbours along axis ax with a on one side and b on the other."""
    s0 = [slice(None)] * 3; s1 = [slice(None)] * 3
    s0[ax] = slice(0, -1); s1[ax] = slice(1, None)
    return (a[tuple(s0)] & b[tuple(s1)]) | (b[tuple(s0)] & a[tuple(s1)]), tuple(s0), tuple(s1)


def contact_stats(ni, ysz, pore, spacing):
    near_pore = ndi.binary_dilation(pore, S6)
    n_tot = 0; n_open = 0; area_nm2 = 0.0
    vz, vy, vx = spacing
    face_area = {0: vy * vx, 1: vz * vx, 2: vz * vy}
    for ax in range(3):
        f, s0, s1 = faces(ni, ysz, ax)
        op = f & (near_pore[s0] | near_pore[s1])
        n_tot += int(f.sum()); n_open += int(op.sum())
        area_nm2 += f.sum() * face_area[ax]
    vol_um3 = ni.size * vz * vy * vx / 1e9
    return n_tot, n_open, area_nm2 / 1e6 / vol_um3


def imread_retry(path, tries=5):
    for i in range(tries):
        try:
            return tifffile.imread(path)
        except OSError:
            if i == tries - 1:
                raise
            time.sleep(2)


rows = []
done = set()
part = os.path.join(OUT, "subblock_scan_partial.csv")
if os.path.exists(part):
    prev = pd.read_csv(part); rows = prev.to_dict("records"); done = set(prev["sample"])
for key, folder, grain, state, nx, ny, nz, vx, vy, vz in SAMPLES:
    if key in done:
        print(f"{key}: cached", flush=True); continue
    t0 = time.time()
    mp = assign_labels(label_histogram(folder)["counts"], ZENODO_LABEL_NOTE[key])
    rois = tile_rois(nz, ny, nx, vz, vy, vx, SIDE_UM)
    ps = slice_paths(folder)
    bands = sorted({(r["z0"], r["z1"]) for r in rois})
    for z0, z1 in bands:
        band = np.stack([imread_retry(ps[z]) for z in range(z0, z1)])
        for r in [r for r in rois if r["z0"] == z0]:
            v = band[:, r["y0"]:r["y1"], r["x0"]:r["x1"]]
            ni = v == mp["Ni"]; ysz = v == mp["YSZ"]; pore = v == mp["pore"]
            ps_ni, per = pspan_mean(ni)
            ps_ysz, _ = pspan_mean(ysz)
            tpb = tpb_density_volume(v, {"Ni": mp["Ni"], "YSZ": mp["YSZ"], "pore": mp["pore"]},
                                     (vz, vy, vx))["tpb_density_um-2"]
            ntot, nopen, sv = contact_stats(ni, ysz, pore, (vz, vy, vx))
            rows.append(dict(sample=key, grain=grain, state=state, roi=r["roi"],
                             phi_ni=ni.mean(), phi_ysz=ysz.mean(), phi_pore=pore.mean(),
                             pspan_ni=ps_ni, pspan_ysz=ps_ysz, tpb=tpb,
                             sv_niysz=sv, open_frac=nopen / max(ntot, 1)))
        del band
    pd.DataFrame(rows).to_csv(part, index=False)
    print(f"{key}: {len(rois)} blocks [{time.time()-t0:.0f}s]", flush=True)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "subblock_scan.csv"), index=False)
print(df.groupby("sample")[["phi_ni", "phi_ysz", "phi_pore", "pspan_ni", "pspan_ysz",
                           "tpb", "sv_niysz", "open_frac"]].agg(["mean", "std"]).round(3).to_string())
