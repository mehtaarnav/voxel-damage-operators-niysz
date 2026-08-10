"""
Ni-phase skeleton -> weighted graph, using skan.

CONVENTIONS (all deliberate; a wrong silent choice here would invalidate the
whole study, so each is stated and justified)

1. WHICH VOXELS ENTER THE GRAPH.
   The graph is built on the PERCOLATING Ni cluster only -- the connected
   component(s) that span the domain face-to-face along the transport axis.
   Rationale: the hypothesis is about the margin of the *electron-carrying*
   network.  Isolated Ni islands carry no current and would otherwise add
   spurious disconnected graph components that make algebraic connectivity
   identically zero.
   Connectivity for this step is 6-connectivity (face-sharing), as validated in
   Phase 0.

2. TRANSPORT AXIS.
   x (array axis 2).  In these datasets x is the through-thickness direction:
   the papers' Figures 11 and S2 plot profiles against "Distance x-coord" over
   0-20 um, which matches the x extent of the image windows (e.g. 19.43 um for
   the fine sample), and the connectivity check is described as being made
   "with the inlet plane on the left (x-direction)".

3. SKELETONIZATION.
   skimage.morphology.skeletonize on the 3D binary mask (Lee et al. 1994
   3D thinning).  NOTE: skeletonize has no `spacing` argument and therefore
   assumes isotropic voxels.  For the three PRISTINE stacks the anisotropy is
   2.4-3.0 % (e.g. 19.53 x 19.53 x 20.00 nm), which is negligible.  For the
   post-redox medium/coarse stacks it is 40 % (17.9 x 17.9 x 25 nm) and this is
   reported rather than silently ignored; those stacks are not skeletonized in
   the pristine-state analysis.

4. LOCAL RADIUS.
   scipy.ndimage.distance_transform_edt on the SAME mask that is skeletonized,
   with `sampling` set to the true physical voxel size in nm.  The value at a
   voxel is the Euclidean distance to the nearest voxel outside the mask, i.e.
   the radius of the largest inscribed sphere centred there.
   BIAS NOTE: because the nearest *background voxel centre* is one voxel beyond
   the last foreground voxel, this over-estimates the true half-width by about
   half a voxel.  Verified analytically in probe_skan2.py: a bar 5 voxels wide
   at 10 nm spacing returns 30.0 nm rather than the geometric 25 nm.  The bias
   is a constant additive ~0.5 voxel and does not affect RANKINGS between
   samples of similar voxel size, but absolute neck widths carry it.

5. NECK WIDTH OF A BRANCH (the edge weight).
   w_e = 2 * min(EDT along the branch path), in nm.  The minimum along the path
   is the branch's narrowest point; doubling converts inscribed radius to a
   width/diameter.  This is the "neck width proxy" of the hypothesis.

6. ATTACHING EDT VALUES TO skan BRANCHES.
   skan 0.13.1's `Skeleton(..., source_image=...)` does NOT populate
   `pixel_values` (verified in probe_skan.py: it returns None).  The working
   route, verified against an analytic shape in probe_skan2.py, is to pass the
   skeleton as a FLOAT image whose nonzero values are the EDT.

7. PARALLEL EDGES.
   skan returns a MultiGraph (two junctions can be joined by several distinct
   branches).  Collapsing to a simple graph, parallel branches are combined by
   SUMMING their conductances, which is the correct series/parallel rule for
   resistors in parallel.
"""

from __future__ import annotations

import numpy as np
import networkx as nx
from scipy import ndimage as ndi
from skan import csr
from skimage.morphology import skeletonize

from .percolation import percolating_mask, percolation_report


def build_ni_graph(ni_mask: np.ndarray,
                   spacing_nm: tuple[float, float, float],
                   axis: int = 2,
                   connectivity: int = 6,
                   restrict_to_percolating: bool = True):
    """Ni mask -> (networkx.Graph, diagnostics dict, arrays dict).

    Parameters
    ----------
    ni_mask   : 3D bool, True where Ni
    spacing_nm: physical voxel size (dz, dy, dx) in nm, matching array axes
    axis      : transport axis for the spanning test (2 = x = through-thickness)
    """
    diag = {}
    rep = percolation_report(ni_mask, axis=axis, connectivity=connectivity)
    diag.update({f"perc_{k}": v for k, v in rep.items()})

    if restrict_to_percolating:
        mask = percolating_mask(ni_mask, axis=axis, connectivity=connectivity)
        if not mask.any():
            diag["graph_built"] = False
            diag["reason"] = "no spanning Ni cluster"
            return None, diag, {}
    else:
        mask = ni_mask

    diag["mask_voxels"] = int(mask.sum())
    diag["mask_frac_of_domain"] = float(mask.mean())

    edt = ndi.distance_transform_edt(mask, sampling=spacing_nm)
    skel = skeletonize(mask)
    diag["skeleton_voxels"] = int(skel.sum())
    if skel.sum() < 10:
        diag["graph_built"] = False
        diag["reason"] = "skeleton too small"
        return None, diag, {"edt": edt, "skel": skel, "mask": mask}

    skel_f = skel.astype(np.float64) * edt          # convention 6
    S = csr.Skeleton(skel_f, spacing=spacing_nm)
    summ = csr.summarize(S, separator="-")
    diag["n_paths"] = int(S.n_paths)

    Gm = csr.skeleton_to_nx(S, summ)                # MultiGraph
    diag["multigraph_nodes"] = Gm.number_of_nodes()
    diag["multigraph_edges"] = Gm.number_of_edges()

    # ---- collapse MultiGraph -> weighted simple Graph -------------------
    # edge attributes:
    #   neck_nm    : 2 * min(EDT) along the branch   (neck WIDTH)
    #   length_nm  : branch arc length
    #   cond       : conductance, see weight_edges()
    G = nx.Graph()
    for u, v, k, d in Gm.edges(keys=True, data=True):
        if u == v:
            continue                                 # drop self-loops
        vals = d.get("values")
        if vals is None or len(vals) == 0:
            continue
        neck = 2.0 * float(np.min(vals))
        path = d.get("path")
        if path is not None and len(path) > 1:
            steps = np.diff(np.asarray(path, dtype=float), axis=0)
            steps *= np.asarray(spacing_nm, dtype=float)[None, :]
            length = float(np.sqrt((steps ** 2).sum(axis=1)).sum())
        else:
            length = float(np.mean(spacing_nm))
        length = max(length, float(np.min(spacing_nm)))

        if G.has_edge(u, v):
            e = G[u][v]
            # parallel branches: keep the widest neck, sum conductance later
            e["neck_nm"] = max(e["neck_nm"], neck)
            e["length_nm"] = min(e["length_nm"], length)
            e["multiplicity"] += 1
        else:
            G.add_edge(u, v, neck_nm=neck, length_nm=length, multiplicity=1)

    # node coordinates (physical, nm) for face detection later
    coords = S.coordinates
    for n in G.nodes:
        if n < len(coords):
            c = coords[n]
            G.nodes[n]["z_nm"] = float(c[0])
            G.nodes[n]["y_nm"] = float(c[1])
            G.nodes[n]["x_nm"] = float(c[2])

    diag["graph_nodes"] = G.number_of_nodes()
    diag["graph_edges"] = G.number_of_edges()
    diag["graph_built"] = G.number_of_edges() > 0
    return G, diag, {"edt": edt, "skel": skel, "mask": mask}


def weight_edges(G: nx.Graph) -> nx.Graph:
    """Attach the conductance used as the graph weight.

    CONVENTION (documented, no standard library exists for this):
        A cylindrical conductor of cross-section A and length L has conductance
        G = sigma * A / L.  Taking the branch's narrowest point as the limiting
        cross-section, A = pi * (w/2)^2 with w the neck WIDTH, so

            cond  =  (w/2)^2 / L        [nm, sigma and pi dropped]

        The constant factors sigma and pi are dropped because every metric we
        report is either a ranking or is normalised; they would multiply every
        edge identically.  Units are nm.
    This is the "neck-conductance weighting" the study calls for; it is a
    modelling choice, not a measurement, and the alternative of weighting by
    neck width alone (cond = w) is also computed as `cond_width` so the
    sensitivity of the conclusion to this choice can be checked.
    """
    for _, _, d in G.edges(data=True):
        w = max(d["neck_nm"], 1e-9)
        L = max(d["length_nm"], 1e-9)
        m = d.get("multiplicity", 1)
        d["cond"] = m * (w * 0.5) ** 2 / L
        d["cond_width"] = m * w
    return G


def largest_component(G: nx.Graph) -> nx.Graph:
    if G.number_of_nodes() == 0:
        return G
    comps = list(nx.connected_components(G))
    biggest = max(comps, key=len)
    return G.subgraph(biggest).copy()
