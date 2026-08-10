"""
Five reusable wrapper functions for the synthetic decoupling study.

Every function is a THIN wrapper: it does not redefine any scientific
convention, it only adapts the existing, already-validated real-data pipeline
(`cmlib.percolation`, `cmlib.pnm`, `cmlib.tpb`, `cmlib.particles`,
`cmlib.metrics`) to the ternary `ganLabel = {0:pore, 1:YSZ, 2:Ni}` convention
used by `cmlib.synthvol.TernaryVolume`, so Phase 0-6 code and Phase-1-onward
synthetic code share one call surface and cannot silently diverge in
definition. `labels` dicts passed to `compute_percolation`/`compute_tpb` use
the SAME convention already established in `cmlib.tpb` docstrings/call sites:
`{"Ni": <label>, "YSZ": <label>, "pore": <label>}`.

For the synthetic study, `cmlib.synthvol.LABELS` (`{"pore":0,"YSZ":1,"Ni":2}`)
is the fixed convention and should be passed as `labels` everywhere.
"""

from __future__ import annotations

import networkx as nx
import numpy as np

from .metrics import summarise_network
from .particles import cpsd_r50max, size_stats, watershed_particles
from .pnm import extract_ni_network, largest_component
from .percolation import percolation_summary
from .tpb import tpb_density_volume


def extract_network(binary_ni: np.ndarray, spacing_nm, *, axis: int = 2,
                    r_max: int = 4, sigma: float = 0.4,
                    connectivity: int = 6):
    """Ni mask -> (largest-component Graph or None, diagnostics dict).

    Thin wrapper over cmlib.pnm.extract_ni_network + largest_component. The
    returned graph (if any) has `face_lo`/`face_hi` node-id lists attached as
    graph attributes, restricted to the largest connected component, exactly
    as phase3_extract_network.py does for the real data -- so
    `compute_network_metrics` can be called with no further face bookkeeping.
    """
    G, diag, extras = extract_ni_network(binary_ni, spacing_nm, axis=axis,
                                         connectivity=connectivity,
                                         r_max=r_max, sigma=sigma)
    if G is None or G.number_of_edges() == 0:
        return None, diag
    Gc = largest_component(G)
    Gc.graph["face_lo"] = [n for n in extras.get("face_lo", []) if n in Gc]
    Gc.graph["face_hi"] = [n for n in extras.get("face_hi", []) if n in Gc]
    diag["lcc_nodes"] = Gc.number_of_nodes()
    diag["lcc_edges"] = Gc.number_of_edges()
    return Gc, diag


def compute_percolation(ternary_volume: np.ndarray, labels: dict, *,
                        phase: str = "Ni", axis: int = 2,
                        connectivity: int = 6) -> dict:
    """Percolation summary (P_span, P_reach, P_largest, ...) for one phase.

    `labels` : {"Ni": v, "YSZ": v, "pore": v} label-value mapping.
    `phase`  : which of those three phases to test (default "Ni", the
              electronic-percolation quantity this whole study is about).
    Thin wrapper over cmlib.percolation.percolation_summary.
    """
    mask = ternary_volume == labels[phase]
    return percolation_summary(mask, axis=axis, connectivity=connectivity)


def compute_tpb(ternary_volume: np.ndarray, labels: dict, spacing_nm) -> dict:
    """TPB length density. Thin wrapper over cmlib.tpb.tpb_density_volume.

    `labels` must be the {"Ni","YSZ","pore"} mapping (cmlib.tpb's own
    convention); passed straight through, unchanged.
    """
    return tpb_density_volume(ternary_volume, labels, spacing_nm)


def compute_particle_stats(ni_volume: np.ndarray, spacing_nm, *,
                           min_distance: int = 4, sigma: float = 0.4,
                           voxel_nm: float | None = None) -> dict:
    """Ni particle size by BOTH measures: watershed (primary) and c-PSD
    (neck-insensitive; see cmlib.particles module docstring for why both are
    reported for the synthetic study specifically).

    `voxel_nm`: scalar isotropic voxel size for the c-PSD call. If omitted,
    uses the geometric mean of `spacing_nm` (same convention as
    cmlib.pnm.geometric_voxel_size_nm), but c-PSD is only meaningful for
    (near-)isotropic voxels -- see cmlib/particles.py.
    """
    labels_img, edt, n_peaks = watershed_particles(
        ni_volume, spacing_nm, min_distance=min_distance, sigma_vox=sigma)
    ws_stats, _ = size_stats(labels_img, spacing_nm)

    if voxel_nm is None:
        voxel_nm = float(np.prod(np.asarray(spacing_nm, dtype=float))
                         ** (1.0 / 3.0))
    cpsd_stats = cpsd_r50max(ni_volume, voxel_nm)

    out = {"n_peaks": n_peaks, "voxel_nm_used_for_cpsd": voxel_nm}
    out.update({f"ws_{k}": v for k, v in ws_stats.items()})
    out.update(cpsd_stats)
    return out


def compute_network_metrics(G: nx.Graph | None, face_lo=None,
                            face_hi=None) -> dict:
    """All connectivity-margin metrics for one network.

    Thin wrapper over cmlib.metrics.summarise_network. Returns an all-NaN /
    zero-node dict (not an exception) if G is None, so this can be called
    unconditionally in a batch loop over many structures, some of which may
    have failed to yield a network (e.g. non-percolating damaged structures).
    """
    if G is None:
        return {"n_nodes": 0, "n_edges": 0, "mean_degree": float("nan"),
               "n_face_lo": 0, "n_face_hi": 0,
               "lambda2_raw": float("nan"), "lambda2_norm": float("nan"),
               "mincut": float("nan"), "mincut_note": "no graph",
               "g_eff": float("nan"),
               "neck_p10_nm": float("nan"), "neck_p25_nm": float("nan"),
               "neck_p50_nm": float("nan"), "neck_p90_nm": float("nan")}
    return summarise_network(G, face_lo, face_hi)
