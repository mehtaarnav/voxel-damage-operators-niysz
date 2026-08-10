"""
Connectivity-margin metrics on the Ni network.

THE SIZE-DEPENDENCE PROBLEM (stated up front because it decides what is
legitimately comparable between anodes)
---------------------------------------------------------------------------
The three anodes yield networks with very different node counts in the same
physical volume (roughly 490 / 140 / 60 nodes per 8 um cube for fine / medium /
coarse) -- that IS the microstructure, not an artifact.  But the raw algebraic
connectivity lambda_2 of a weighted graph Laplacian depends on both the edge
weights AND the number of nodes, so comparing raw lambda_2 across anodes mixes
a genuine connectivity difference with a discretisation difference.

We therefore report FOUR quantities and say which are safe to compare:

  lambda2_raw      2nd-smallest eigenvalue of the weighted Laplacian.
                   The quantity named in the hypothesis.  SIZE-DEPENDENT --
                   report it, but do not rank on it alone.
  lambda2_norm     2nd-smallest eigenvalue of the NORMALISED Laplacian,
                   which lies in [0, 2] and is far less sensitive to node
                   count.  Safe(r) to compare.
  mincut           max-flow / min-cut between the two opposite faces of the
                   transport axis, with edge capacity = conductance.
                   Because every ROI is the SAME PHYSICAL SIZE, this is an
                   effective bottleneck conductance and IS comparable.
  g_eff            the actual effective conductance of the resistor network
                   between the two faces, from a Laplacian solve with Dirichlet
                   boundary conditions.  This is the physically correct
                   through-transport quantity; min-cut is only an upper bound
                   on it (min-cut ignores series resistance).  Included because
                   it is nearly free and is what "retained percolation" should
                   actually track.

CHEEGER CONTEXT.  For a weighted graph, lambda_2/2 <= h <= sqrt(2*lambda_2*d_max)
where h is the Cheeger constant.  So lambda_2 and a normalised min-cut are two
views of the same bottleneck; reporting both is a consistency check, not
duplication.

UNITS.  Edge conductance is cross-sectional area / length in nm (the intrinsic
conductivity sigma is dropped; it multiplies every edge identically and cancels
from every ranking).  lambda_2 and mincut therefore carry units of nm.
"""

from __future__ import annotations

import networkx as nx
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

WEIGHT = "cond"


def _laplacian(G, nodes, weight=WEIGHT):
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    rows, cols, vals = [], [], []
    deg = np.zeros(n)
    for u, v, d in G.edges(data=True):
        if u not in idx or v not in idx:
            continue
        w = float(d.get(weight, 1.0))
        i, j = idx[u], idx[v]
        rows += [i, j]; cols += [j, i]; vals += [-w, -w]
        deg[i] += w; deg[j] += w
    rows += list(range(n)); cols += list(range(n)); vals += list(deg)
    L = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    return L, deg, idx


def algebraic_connectivity(G, weight=WEIGHT, normalized=False):
    """lambda_2 of the (optionally normalised) weighted Laplacian."""
    nodes = list(G.nodes)
    n = len(nodes)
    if n < 2 or G.number_of_edges() == 0:
        return float("nan")
    if not nx.is_connected(G):
        return 0.0
    try:
        return float(nx.algebraic_connectivity(
            G, weight=weight, normalized=normalized, method="tracemin_pcg"))
    except Exception:
        # dense fallback for small graphs / solver trouble
        L, deg, _ = _laplacian(G, nodes, weight)
        if normalized:
            with np.errstate(divide="ignore"):
                dm = 1.0 / np.sqrt(np.where(deg > 0, deg, 1.0))
            D = sp.diags(dm)
            L = D @ L @ D
        try:
            vals = np.sort(np.linalg.eigvalsh(L.toarray()))
            return float(vals[1])
        except Exception:
            return float("nan")


def mincut_between_faces(G, face_lo, face_hi, weight=WEIGHT):
    """Max-flow (= min-cut capacity) between two face node-sets.

    A super-source and super-sink are attached with effectively infinite
    capacity, so the cut is forced to fall inside the network rather than on
    the terminal links.
    """
    face_lo = [n for n in face_lo if n in G]
    face_hi = [n for n in face_hi if n in G]
    if not face_lo or not face_hi:
        return float("nan"), "empty face set"
    overlap = set(face_lo) & set(face_hi)
    if overlap:
        # a node touching both faces short-circuits the domain
        return float("inf"), f"{len(overlap)} node(s) touch both faces"

    H = nx.Graph()
    for u, v, d in G.edges(data=True):
        H.add_edge(u, v, capacity=float(d.get(weight, 1.0)))
    if H.number_of_edges() == 0:
        return float("nan"), "no edges"
    big = 1e12 * max(1.0, max(d["capacity"] for _, _, d in H.edges(data=True)))
    S, T = "__S__", "__T__"
    for n in face_lo:
        H.add_edge(S, n, capacity=big)
    for n in face_hi:
        H.add_edge(n, T, capacity=big)
    if not nx.has_path(H, S, T):
        return 0.0, "no path between faces"
    val = nx.maximum_flow_value(H, S, T, capacity="capacity")
    return float(val), "ok"


def effective_conductance(G, face_lo, face_hi, weight=WEIGHT):
    """Effective conductance between the two faces (Dirichlet Laplacian solve).

    Sets V=1 on face_lo, V=0 on face_hi, solves for interior potentials, and
    returns the total current, which equals the effective conductance for a
    unit potential difference.
    """
    face_lo = sorted(set(n for n in face_lo if n in G))
    face_hi = sorted(set(n for n in face_hi if n in G))
    if not face_lo or not face_hi or set(face_lo) & set(face_hi):
        return float("nan")
    nodes = list(G.nodes)
    L, deg, idx = _laplacian(G, nodes, weight)
    n = len(nodes)
    fixed = np.zeros(n, dtype=bool)
    V = np.zeros(n)
    for m in face_lo:
        fixed[idx[m]] = True; V[idx[m]] = 1.0
    for m in face_hi:
        fixed[idx[m]] = True; V[idx[m]] = 0.0
    free = ~fixed
    if free.sum() == 0:
        A = 0.0
    else:
        Lff = L[free][:, free]
        Lfx = L[free][:, fixed]
        rhs = -Lfx @ V[fixed]
        try:
            Vf = spl.spsolve(Lff.tocsc(), rhs)
        except Exception:
            Vf, _ = spl.cg(Lff.tocsr(), rhs, rtol=1e-10, maxiter=20000)
        V[free] = Vf
    # current out of the source set
    I = 0.0
    src = set(face_lo)
    for u, v, d in G.edges(data=True):
        w = float(d.get(weight, 1.0))
        if (u in src) ^ (v in src):
            a, b = (u, v) if u in src else (v, u)
            I += w * (V[idx[a]] - V[idx[b]])
    return float(I)


def neck_quantiles(G, qs=(10, 25, 50, 90)):
    necks = np.array([d["neck_nm"] for _, _, d in G.edges(data=True)])
    if necks.size == 0:
        return {f"neck_p{q}_nm": float("nan") for q in qs}
    return {f"neck_p{q}_nm": float(np.percentile(necks, q)) for q in qs}


def summarise_network(G, face_lo=None, face_hi=None):
    """All connectivity-margin metrics for one network."""
    face_lo = face_lo if face_lo is not None else G.graph.get("face_lo", [])
    face_hi = face_hi if face_hi is not None else G.graph.get("face_hi", [])
    out = {
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "mean_degree": (2.0 * G.number_of_edges() / G.number_of_nodes()
                        if G.number_of_nodes() else float("nan")),
        "n_face_lo": len(face_lo),
        "n_face_hi": len(face_hi),
    }
    out["lambda2_raw"] = algebraic_connectivity(G, normalized=False)
    out["lambda2_norm"] = algebraic_connectivity(G, normalized=True)
    mc, note = mincut_between_faces(G, face_lo, face_hi)
    out["mincut"] = mc
    out["mincut_note"] = note
    out["g_eff"] = effective_conductance(G, face_lo, face_hi)
    out.update(neck_quantiles(G))
    vols = np.array([d.get("volume_nm3", np.nan)
                     for _, d in G.nodes(data=True)], dtype=float)
    if np.isfinite(vols).any():
        # equivalent-sphere diameter of the median chamber
        vmed = float(np.nanmedian(vols))
        out["chamber_vol_median_nm3"] = vmed
        out["chamber_equiv_diam_median_nm"] = float(
            (6.0 * vmed / np.pi) ** (1.0 / 3.0))
    dia = np.array([d.get("equiv_diam_nm", np.nan)
                    for _, d in G.nodes(data=True)], dtype=float)
    if np.isfinite(dia).any():
        out["chamber_equiv_diam_mean_nm"] = float(np.nanmean(dia))
    return out
