"""
Ni Platform v3 -- Option B: extract real Ni particle/neck graphs from the
Holzer/Pecho segmentation and run the FROZEN vulnerability audit on them.

Pre-registered in out/project2/PREREG_NI_PLATFORM_V3.md (commit 50b0bde),
committed before this file was written. Audit metrics are those already frozen
in PREREG_NI_VULNERABILITY_AUDIT.md (267efd3) and are not redefined here.

Read-only with respect to the generator: cmlib/damage.py, cmlib/synth.py and
cmlib/project2.py are untouched. No damage operator is implemented.

Purpose: the lattice platform's minimum cut is exactly one full cross-section
with ZERO seed variance. That is a property of the lattice. This measures what a
REAL Ni network's cut structure looks like, which is both the test of that
diagnosis and the target any future disordered generator must hit.
"""
from __future__ import annotations

import os
import sys
import time

import networkx as nx
import numpy as np
import pandas as pd
import tifffile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE  # noqa: E402
from cmlib.io import label_histogram, slice_paths  # noqa: E402
from cmlib.phases import assign_labels  # noqa: E402
from cmlib.pnm import extract_ni_network  # noqa: E402
from cmlib.roi import tile_rois  # noqa: E402

OUT = os.path.join(ROOT, "out", "project2")
AXIS, CONN = 2, 6          # transport axis = x, real-data convention
SIDE_UM = 8.0
N_ROI = 3
INF = 10 ** 9


def load_roi(folder, ni_label, r):
    """`r` is a tile_rois dict with z0/z1/y0/y1/x0/x1."""
    ps = slice_paths(folder)
    out = np.empty((r["z1"] - r["z0"], r["y1"] - r["y0"], r["x1"] - r["x0"]),
                   dtype=bool)
    for i, z in enumerate(range(r["z0"], r["z1"])):
        a = tifffile.imread(ps[z])
        out[i] = a[r["y0"]:r["y1"], r["x0"]:r["x1"]] == ni_label
    return out


def build_st_graph(G, extras):
    """Add virtual S/T on the two faces of the transport axis, unit capacities.

    S/T edges get infinite capacity and idx=-1 so they are never cut and never
    counted -- identical convention to the synthetic audit.
    """
    H = nx.Graph()
    for n in G.nodes:
        H.add_node(n)
    for i, (u, v, d) in enumerate(G.edges(data=True)):
        H.add_edge(u, v, idx=i, capacity=1, neck_nm=d["neck_nm"],
                   area_nm2=d["area_nm2"])
    H.add_node("S")
    H.add_node("T")
    for n in extras["face_lo"]:
        if n in H:
            H.add_edge("S", n, idx=-1, capacity=INF)
    for n in extras["face_hi"]:
        if n in H:
            H.add_edge(n, "T", idx=-1, capacity=INF)
    return H


def min_cut_edges(H):
    val, (side_s, _) = nx.minimum_cut(H, "S", "T", capacity="capacity")
    idx = {d["idx"] for u, v, d in H.edges(data=True)
           if d["idx"] >= 0 and ((u in side_s) != (v in side_s))}
    return int(val), idx


def spans_after(H, removed):
    K = H.copy()
    K.remove_edges_from([(u, v) for u, v, d in H.edges(data=True)
                         if d["idx"] in removed])
    return nx.has_path(K, "S", "T")


def main():
    rows, curves = [], []
    for key, folder, grain, state, nz, ny, nx_, vz, vy, vx in SAMPLES:
        if state != "pristine":
            continue
        counts = label_histogram(folder)["counts"]
        mapping = assign_labels(counts, ZENODO_LABEL_NOTE[key])
        rois = tile_rois(nz, ny, nx_, vz, vy, vx, SIDE_UM, max_rois=N_ROI)
        for ri, sl in enumerate(rois):
            t0 = time.time()
            mask = load_roi(folder, mapping["Ni"], sl)
            G, diag, extras = extract_ni_network(
                mask, spacing_nm=(vz, vy, vx), axis=AXIS, connectivity=CONN)
            if G is None or G.number_of_edges() == 0:
                print(f"  {grain} roi{ri}: NO NETWORK ({diag.get('reason')})",
                      flush=True)
                continue
            H = build_st_graph(G, extras)
            if not (nx.has_path(H, "S", "T")):
                print(f"  {grain} roi{ri}: no S-T path, skipped", flush=True)
                continue

            # --- positional disorder diagnostics ---
            cen = {n: np.array([G.nodes[n]["z_nm"], G.nodes[n]["y_nm"],
                                G.nodes[n]["x_nm"]]) for n in G.nodes}
            dists = np.array([np.linalg.norm(cen[u] - cen[v])
                              for u, v in G.edges])
            degs = np.array([d for _n, d in G.degree()])
            necks = np.array([d["neck_nm"] for _u, _v, d in G.edges(data=True)])
            # lattice-likeness: coefficient of variation of pair distances
            cv_dist = float(dists.std() / dists.mean())

            n_edges = G.number_of_edges()
            cut_val, idx_cut = min_cut_edges(H)
            f_cut = len(idx_cut) / n_edges
            q25 = float(np.percentile(necks, 25))
            lowq = {i for i, nv in enumerate(necks) if nv < q25}
            overlap = len(idx_cut & lowq) / max(len(idx_cut), 1)

            rows.append(dict(
                anode=grain, roi=ri, n_nodes=G.number_of_nodes(),
                n_edges=n_edges, mincut_edges=len(idx_cut),
                frac_mincut=f_cut, mean_degree=float(degs.mean()),
                sd_degree=float(degs.std()), min_degree=int(degs.min()),
                max_degree=int(degs.max()),
                mean_pair_dist_nm=float(dists.mean()),
                sd_pair_dist_nm=float(dists.std()), cv_pair_dist=cv_dist,
                neck_p25_nm=q25, overlap_mincut_lowerquartile=overlap,
                seconds=round(time.time() - t0, 1)))
            print(f"  {grain:7s} roi{ri} nodes={G.number_of_nodes()} "
                  f"edges={n_edges} mincut={len(idx_cut)} ({f_cut:.4f}) "
                  f"deg={degs.mean():.2f}+-{degs.std():.2f} "
                  f"CV(dist)={cv_dist:.3f} overlap={overlap:.3f} "
                  f"[{rows[-1]['seconds']}s]", flush=True)
            pd.DataFrame(rows).to_csv(
                os.path.join(OUT, "v3_real_ni_graph_audit.csv"), index=False)

            # --- fragility curves: random and low-neck-area removal ---
            order_low = list(np.argsort(necks))
            for frac in (0.02, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50):
                k = int(frac * n_edges)
                low_sp = spans_after(H, set(int(i) for i in order_low[:k]))
                rnd = []
                for rep in range(3):
                    rg = np.random.default_rng(7000 + rep)
                    rnd.append(spans_after(H, set(
                        int(i) for i in rg.choice(n_edges, size=k,
                                                  replace=False))))
                curves.append(dict(anode=grain, roi=ri, frac_removed=frac,
                                   lowarea_spans=low_sp,
                                   random_span_rate=float(np.mean(rnd))))
                pd.DataFrame(curves).to_csv(
                    os.path.join(OUT, "v3_real_fragility_curves.csv"),
                    index=False)
            del mask

    df = pd.DataFrame(rows)
    if not len(df):
        print("no ROIs produced a network")
        return 1
    print("\n" + "=" * 74)
    print("REAL Ni GRAPH AUDIT -- anode medians")
    print("=" * 74)
    g = df.groupby("anode")[["n_nodes", "n_edges", "mincut_edges",
                             "frac_mincut", "mean_degree", "sd_degree",
                             "cv_pair_dist", "overlap_mincut_lowerquartile"]]
    print(g.median().to_string())
    print("\nfrac_mincut per-ROI ranges (lattice had EXACTLY zero variance):")
    for a in df.anode.unique():
        v = df[df.anode == a].frac_mincut
        print(f"  {a:7s} {v.min():.4f} - {v.max():.4f}  (n={len(v)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
