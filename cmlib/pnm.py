"""
Ni-phase network extraction by watershed partitioning (SNOW), replacing the
skeleton route.

WHY THIS REPLACED THE SKELETON ROUTE
------------------------------------
skimage.morphology.skeletonize returned a proper CURVE skeleton for the fine
anode (83 % of skeleton voxels have degree 2) but a MEDIAL SHEET for the medium
and coarse anodes (79 % and 94 % of voxels have degree >= 4, with a spike at
degree exactly 8 -- a planar signature).  Median branch length was 6.7 voxels
for fine but 2.0 for medium and coarse.  Thinning must preserve topology, and
the thick multiply-connected Ni domains of the coarser anodes have a genuinely
2D medial axis.  Because the skeleton's dimensionality varied monotonically
with grain coarseness, any cross-anode ranking derived from it would be an
artifact confounded with exactly the variable the hypothesis is tested against.
See out/phase3/phase3_GATE_FAILURE_skeleton_dimensionality.png.

THE REPLACEMENT
---------------
SNOW (Sub-Network of an Over-segmented Watershed), Gostick, Phys. Rev. E 96,
023307 (2017), as implemented in porespy.networks.snow2.  Applied to the SOLID
(Ni) phase rather than the pore phase -- the algorithm is phase-agnostic.

    nodes  = Ni "chambers", i.e. watershed regions of the Euclidean distance
             transform (one region per distance-transform peak)
    edges  = shared interfaces between adjacent chambers
    weight = measured directly on the shared interface

This is robust to blobby vs strut-like morphology because it never has to
reduce a solid to a curve.

WATERSHED MARKER PARAMETERS (stated explicitly -- these are exactly the
"marker parameters" that must not be chosen silently)
-----------------------------------------------------------------------
    sigma  = 0.4   Gaussian blur applied to the distance transform before peak
                   detection.  Suppresses single-voxel noise peaks that would
                   over-segment.  porespy/SNOW default.
    r_max  = 4     radius (voxels) of the maximum filter used to merge nearby
                   peaks into a single marker.  This is the dominant
                   over/under-segmentation control.  porespy/SNOW default.
    accuracy = 'standard'
    boundary_width = 0   no artificial boundary pores are added; face contact is
                   determined afterwards from the region image, so that the
                   min-cut source/sink sets are ours and not porespy's.
Both sigma and r_max are in VOXELS, so their physical meaning differs between
samples of different voxel size; a sensitivity sweep over r_max is provided in
phase3_snow_sensitivity.py rather than assuming the default is harmless.

EDGE WEIGHT CONVENTION (no standard exists; stated and justified)
-----------------------------------------------------------------
A throat of cross-sectional area A and length L has conductance sigma*A/L.  We
use
        cond = throat.cross_sectional_area / throat.total_length      [m]
dropping the intrinsic conductivity sigma, which multiplies every edge
identically and cancels from every ranking.  porespy measures A directly on the
shared interface, verified exact on an analytic test shape (probe_porespy.py:
a 6-voxel bar at 10 nm returns A = 3600 nm^2 and inscribed diameter 60.0 nm).

The NECK WIDTH used for the "lower quantile of neck widths" metric is
        throat.inscribed_diameter
which is the diameter of the largest disc fitting in the throat cross-section.

ANISOTROPY CAVEAT
-----------------
snow2 takes a SCALAR voxel_size.  The three pristine stacks are 2.4-3.0 %
anisotropic (e.g. 19.53 x 19.53 x 20.00 nm), so we pass the geometric mean
(vx*vy*vz)^(1/3), which preserves volumes exactly and mis-states lengths by at
most ~1.5 %.  The post-redox medium and coarse stacks are 40 % anisotropic
(17.9 x 17.9 x 25 nm) and are NOT processed with this scalar assumption.
"""

from __future__ import annotations

import networkx as nx
import numpy as np
import porespy as ps

from .percolation import percolating_mask, percolation_report

SNOW_SIGMA = 0.4
SNOW_RMAX = 4
SNOW_ACCURACY = "standard"


def geometric_voxel_size_nm(spacing_nm) -> float:
    return float(np.prod(np.asarray(spacing_nm, dtype=float)) ** (1.0 / 3.0))


def extract_ni_network(ni_mask: np.ndarray,
                       spacing_nm,
                       axis: int = 2,
                       connectivity: int = 6,
                       r_max: int = SNOW_RMAX,
                       sigma: float = SNOW_SIGMA,
                       parallel_kw=None):
    """Ni mask -> (networkx.Graph, diagnostics, extras).

    Graph node attributes : volume_nm3, inscribed_diam_nm, coords (nm)
    Graph edge attributes : area_nm2, length_nm, neck_nm, cond
    Diagnostics include the percolation report of the input mask.

    PARALLEL_KW / A PORESPY BUG FOUND DURING THE SYNTHETIC-STUDY WORK
    -------------------------------------------------------------------
    `ps.networks.snow2`'s own default is `parallel_kw={}` (an empty dict, NOT
    None), which is truthy-for-"use parallel" in its own dispatch and silently
    routes every call through `snow_partitioning_parallel` (chunked,
    divs=[2,2,2] by default) UNLESS the caller passes `parallel_kw=None`
    explicitly. This function previously did not pass `parallel_kw` at all, so
    every SNOW extraction in the original real-data study (Phase 3/4, all 21
    pristine ROIs, the r_max sensitivity sweep, and the "Image was cropped to
    ..." warnings visible in that study's logs) silently ran in CHUNKED mode.

    Found via a Phase-0 unit test for the synthetic decoupling study
    (scripts/next/phase0_validate_synthetic_pipeline.py T4): an idealized
    two-cube dumbbell with a small connecting bar, run through the DEFAULT
    (chunked) path, crashes inside porespy's regions_to_network with
    `IndexError: too many indices for array: array is 1-dimensional, but 2
    were indexed` -- even though the two watershed regions are genuinely
    6-connectivity face-adjacent (confirmed directly). Forcing
    `parallel_kw=None` (serial) fixes the crash and recovers the exact
    analytic answer. See probe_snow2_parallel_bug.py.

    QUANTIFIED IMPACT ON THE ALREADY-REPORTED REAL-DATA RESULTS: re-running
    extraction on the coarse_pre z0y0x0 ROI (the exact one in
    out/phase3/phase3_snow_8.0um_rmax4.csv) with parallel_kw=None gives
    65 pores / 153 throats / neck_p10=187.9 nm / neck_p50=529.6 nm, versus the
    REPORTED 62 pores / 145 throats / neck_p10=199 nm / neck_p50=516 nm from
    chunked mode -- a ~5-6% shift. This is SMALLER than both the already-
    reported between-ROI spread for coarse_pre (neck_p10 205+/-45 nm) and the
    already-reported r_max sensitivity swing (177-267 nm), so it does not
    change any conclusion in REPORT.md, but it is a real, previously-
    undisclosed source of variability in that result, disclosed here rather
    than silently absorbed. REPORT.md and the out/phase3, out/phase4 CSVs from
    the original study are NOT retroactively altered.

    Going forward (this function's new default, `parallel_kw=None`): serial
    extraction is used for all new (synthetic-study) work, both because it is
    the one verified exact on analytic cases and because the small, spatially
    localized narrow-neck features the decoupling study specifically designs
    (Family B) are exactly the kind of feature a chunk-boundary seam could
    corrupt or clip -- silently, without necessarily crashing.
    """
    diag = {}
    rep = percolation_report(ni_mask, axis=axis, connectivity=connectivity)
    diag.update({f"perc_{k}": v for k, v in rep.items()})

    mask = percolating_mask(ni_mask, axis=axis, connectivity=connectivity)
    if not mask.any():
        diag["network_built"] = False
        diag["reason"] = "no spanning Ni cluster"
        return None, diag, {}
    diag["mask_voxels"] = int(mask.sum())

    vox_nm = geometric_voxel_size_nm(spacing_nm)
    vox_m = vox_nm * 1e-9

    snow = ps.networks.snow2(mask, voxel_size=vox_m, boundary_width=0,
                             sigma=sigma, r_max=r_max, accuracy=SNOW_ACCURACY,
                             parallel_kw=parallel_kw)
    net = snow.network
    regions = snow.regions

    n_pores = int(net["pore.all"].size)
    n_throats = int(net["throat.all"].size)
    diag["n_pores"] = n_pores
    diag["n_throats"] = n_throats
    diag["snow_sigma"] = sigma
    diag["snow_r_max"] = r_max
    diag["snow_parallel_kw"] = "serial" if parallel_kw is None else "chunked"
    diag["voxel_size_nm_used"] = vox_nm

    if n_throats == 0:
        diag["network_built"] = False
        diag["reason"] = "no throats"
        return None, diag, {"regions": regions, "mask": mask}

    G = nx.Graph()
    coords = net["pore.coords"] * 1e9                     # m -> nm
    pvol = net["pore.region_volume"] * 1e27               # m^3 -> nm^3
    pdia = net["pore.inscribed_diameter"] * 1e9           # m -> nm
    pedia = net["pore.equivalent_diameter"] * 1e9
    for i in range(n_pores):
        G.add_node(int(i),
                   z_nm=float(coords[i, 0]), y_nm=float(coords[i, 1]),
                   x_nm=float(coords[i, 2]),
                   volume_nm3=float(pvol[i]),
                   inscribed_diam_nm=float(pdia[i]),
                   equiv_diam_nm=float(pedia[i]))

    conns = net["throat.conns"]
    area = net["throat.cross_sectional_area"] * 1e18      # m^2 -> nm^2
    length = net["throat.total_length"] * 1e9             # m -> nm
    neck = net["throat.inscribed_diameter"].astype(float) * 1e9

    for t in range(n_throats):
        u, v = int(conns[t, 0]), int(conns[t, 1])
        if u == v:
            continue
        A = float(area[t])
        L = max(float(length[t]), vox_nm)
        w = float(neck[t])
        cond = A / L
        if G.has_edge(u, v):
            e = G[u][v]
            e["area_nm2"] += A
            e["cond"] += cond
            e["neck_nm"] = max(e["neck_nm"], w)
        else:
            G.add_edge(u, v, area_nm2=A, length_nm=L, neck_nm=w, cond=cond)

    G.remove_nodes_from([n for n in list(G.nodes) if G.degree(n) == 0])
    diag["graph_nodes"] = G.number_of_nodes()
    diag["graph_edges"] = G.number_of_edges()
    diag["network_built"] = G.number_of_edges() > 0

    # region labels touching the two faces of `axis`, for min-cut terminals
    reg_lo = np.unique(np.take(regions, 0, axis=axis))
    reg_hi = np.unique(np.take(regions, regions.shape[axis] - 1, axis=axis))
    reg_lo = set(int(r) for r in reg_lo if r > 0)
    reg_hi = set(int(r) for r in reg_hi if r > 0)
    # porespy region labels are pore_index + 1
    lo_nodes = sorted(r - 1 for r in reg_lo if (r - 1) in G)
    hi_nodes = sorted(r - 1 for r in reg_hi if (r - 1) in G)
    diag["n_face_lo"] = len(lo_nodes)
    diag["n_face_hi"] = len(hi_nodes)

    extras = {"regions": regions, "mask": mask,
              "face_lo": lo_nodes, "face_hi": hi_nodes}
    return G, diag, extras


def largest_component(G: nx.Graph) -> nx.Graph:
    if G.number_of_nodes() == 0:
        return G
    return G.subgraph(max(nx.connected_components(G), key=len)).copy()
