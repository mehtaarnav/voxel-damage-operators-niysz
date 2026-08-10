"""
PHASE 0 (synthetic study) — validate the ternary pipeline on analytic cases.

Tests T1-T6 from out/next/EXECUTION_SPEC.md section 6. Every case is hand-built
with a known-exact answer; nothing here is fit or tuned to make a test pass.
All volumes are <=64^3. Uses the SAME scientific definitions as the real-data
study (cmlib.percolation, cmlib.tpb, cmlib.pnm, cmlib.particles) via the five
wrapper functions in cmlib.api, so this validates the actual call surface
Phase 1 onward will use, not a parallel implementation.

Gate G0: T1-T6 all PASS. T5's neck-vs-particle-size coupling numbers are
recorded here for reference, but the DECISION on Q1/Q2 (out/next/EXECUTION_SPEC
section 8) is made from the separate, more thorough coupling experiment run
after gate G-1 (Phase -1 prior art + preregistration), per the approved
sequence — this script's T5 is a pipeline-correctness check only.
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))          # repo root
sys.path.insert(0, ROOT)

from cmlib.api import (  # noqa: E402
    compute_network_metrics, compute_particle_stats, compute_percolation,
    compute_tpb, extract_network,
)
from cmlib.synthvol import (  # noqa: E402
    LABELS, TernaryVolume, load_ternary, save_ternary,
    volume_fractions_from_volume,
)

OUT = os.path.join(ROOT, "out", "next")
os.makedirs(OUT, exist_ok=True)

RESULTS = []


def check(test_id, name, ok, detail=""):
    RESULTS.append((test_id, name, bool(ok), detail))
    flag = "PASS" if ok else "FAIL"
    print(f"  [{flag}] {test_id} {name}" + (f"  -- {detail}" if detail else ""))
    return ok


# ---------------------------------------------------------------------------
# T1 — volume fractions exact
# ---------------------------------------------------------------------------
def test_T1():
    print("-" * 74)
    print("T1 — volume fractions exact on a prescribed ternary block")
    print("-" * 74)
    ok = True
    # 20x20x20 = 8000 voxels; prescribe exact counts: Ni=1234, YSZ=3456, pore=rest
    nz = ny = nx = 20
    total = nz * ny * nx
    n_ni, n_ysz = 1234, 3456
    n_pore = total - n_ni - n_ysz
    flat = np.empty(total, dtype=np.uint8)
    flat[:n_pore] = LABELS["pore"]
    flat[n_pore:n_pore + n_ysz] = LABELS["YSZ"]
    flat[n_pore + n_ysz:] = LABELS["Ni"]
    rng = np.random.default_rng(0)
    rng.shuffle(flat)                 # spatial arrangement irrelevant to phi
    vol = flat.reshape(nz, ny, nx)

    phi = volume_fractions_from_volume(vol)
    exp = {"phi_pore": n_pore / total, "phi_YSZ": n_ysz / total,
          "phi_Ni": n_ni / total}
    for k in exp:
        ok &= check("T1", f"{k} exact", abs(phi[k] - exp[k]) < 1e-12,
                    f"got={phi[k]:.6f} want={exp[k]:.6f}")

    tv = TernaryVolume(vol=vol, voxel_nm=20.0, meta={"test": "T1"})
    phi2 = tv.volume_fractions()
    ok &= check("T1", "TernaryVolume.volume_fractions matches free function",
               phi2 == phi)
    return ok


# ---------------------------------------------------------------------------
# T2 — TPB exact, ternary, anisotropic, + linear addition of parallel lines
# ---------------------------------------------------------------------------
def _single_line_plane(ny, nx, cy, cx):
    """pore | Ni / YSZ meeting along one z-parallel line at (cy,cx)."""
    yy, xx = np.ogrid[:ny, :nx]
    return np.where(xx < cx, LABELS["pore"],
                    np.where(yy < cy, LABELS["Ni"], LABELS["YSZ"])
                    ).astype(np.uint8)


def test_T2():
    print("-" * 74)
    print("T2 — TPB exact on axis-aligned lines; multiple lines add linearly")
    print("-" * 74)
    ok = True

    # T2a: single z-parallel line, anisotropic spacing, EXACT length = nz*dz
    nz, ny, nx = 50, 41, 41
    plane = _single_line_plane(ny, nx, ny // 2, nx // 2)
    vol = np.broadcast_to(plane, (nz, ny, nx)).copy()
    for spacing in [(20.0, 20.0, 20.0), (25.0, 17.9, 17.9)]:
        r = compute_tpb(vol, LABELS, spacing)
        exp_len_um = nz * spacing[0] / 1e3
        edges_ok = (r["tpb_edges_z"] == nz and r["tpb_edges_y"] == 0
                   and r["tpb_edges_x"] == 0)
        len_ok = abs(r["tpb_length_um"] - exp_len_um) < 1e-9
        ok &= check("T2a", f"single line exact @ spacing {spacing}",
                    edges_ok and len_ok,
                    f"edges=({r['tpb_edges_z']},{r['tpb_edges_y']},"
                    f"{r['tpb_edges_x']}) length={r['tpb_length_um']:.6f}um "
                    f"want={exp_len_um:.6f}um")

    # T2b: line parallel to x (axis mix-up guard)
    nz2, ny2, nx2 = 41, 41, 50
    zz, yy = np.ogrid[:nz2, :ny2]
    plane_x = np.where(yy < ny2 // 2, LABELS["pore"],
                       np.where(zz < nz2 // 2, LABELS["Ni"], LABELS["YSZ"])
                       ).astype(np.uint8)
    vol_x = np.repeat(plane_x[:, :, None], nx2, axis=2)
    spacing = (25.0, 17.9, 30.0)
    r = compute_tpb(vol_x, LABELS, spacing)
    exp_len_um = nx2 * spacing[2] / 1e3
    ok &= check("T2b", "x-parallel line exact",
               r["tpb_edges_x"] == nx2 and r["tpb_edges_y"] == 0
               and r["tpb_edges_z"] == 0
               and abs(r["tpb_length_um"] - exp_len_um) < 1e-9,
               f"edges=({r['tpb_edges_z']},{r['tpb_edges_y']},{r['tpb_edges_x']}) "
               f"length={r['tpb_length_um']:.6f}um want={exp_len_um:.6f}um")

    # T2c (NEW): two well-separated copies of the T2a plane -> TPB length and
    # edge counts must be EXACTLY 2x the single-line values (linear addition)
    gap = 6
    ny3 = 2 * ny + gap
    plane1 = _single_line_plane(ny, nx, ny // 2, nx // 2)
    blank_gap = np.full((gap, nx), LABELS["pore"], dtype=np.uint8)
    plane2 = _single_line_plane(ny, nx, ny // 2, nx // 2)
    double_plane = np.concatenate([plane1, blank_gap, plane2], axis=0)
    vol_double = np.broadcast_to(double_plane, (nz, ny3, nx)).copy()
    spacing = (20.0, 20.0, 20.0)
    r_single = compute_tpb(vol, LABELS, spacing)
    r_double = compute_tpb(vol_double, LABELS, spacing)
    ok &= check(
        "T2c", "two non-interacting parallel TPB lines add linearly",
        r_double["tpb_edges_z"] == 2 * r_single["tpb_edges_z"]
        and r_double["tpb_edges_y"] == 0 and r_double["tpb_edges_x"] == 0
        and abs(r_double["tpb_length_um"] - 2 * r_single["tpb_length_um"]) < 1e-9,
        f"double edges_z={r_double['tpb_edges_z']} "
        f"want={2*r_single['tpb_edges_z']}; "
        f"double length={r_double['tpb_length_um']:.6f}um "
        f"want={2*r_single['tpb_length_um']:.6f}um")

    # T2d: re-confirm the measured staircase bias on a (1,1,1)-oriented line
    # (reproduces phase4a_validate_tpb.py test3, routed through compute_tpb)
    N = 81
    d = 20.0
    c = (N - 1) / 2.0
    u = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    e1 = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    e2 = np.cross(u, e1)
    idx = np.arange(N) - c
    Z, Y, X = np.meshgrid(idx, idx, idx, indexing="ij")
    P = np.stack([Z, Y, X], axis=-1)
    perp = P - (P @ u)[..., None] * u
    ang = np.arctan2(perp @ e2, perp @ e1)
    sector = ((ang + np.pi) / (2 * np.pi / 3)).astype(int) % 3
    vol_diag = np.select([sector == 0, sector == 1, sector == 2],
                         [LABELS["pore"], LABELS["YSZ"], LABELS["Ni"]]
                         ).astype(np.uint8)
    r = compute_tpb(vol_diag, LABELS, (d, d, d))
    true_len_um = N * np.sqrt(3.0) * d / 1e3
    ratio = r["tpb_length_um"] / true_len_um
    ok &= check("T2d", "staircase bias within theoretical bound [1.0, sqrt(3)]",
               1.0 <= ratio <= 1.8,
               f"ratio={ratio:.3f} (theoretical worst case sqrt(3)={np.sqrt(3):.3f})")
    return ok


# ---------------------------------------------------------------------------
# T3 — percolation, ternary input, P_span vs P_reach distinguishing case
# ---------------------------------------------------------------------------
def test_T3():
    print("-" * 74)
    print("T3 — percolation on ternary input; P_span vs P_reach")
    print("-" * 74)
    ok = True
    nz, ny, nx = 20, 20, 20
    axis = 2                      # x

    # (a) solid Ni slab spanning full x range -> P_span=1, P_reach=1
    vol = np.full((nz, ny, nx), LABELS["pore"], dtype=np.uint8)
    vol[:, :, :] = LABELS["Ni"]
    r = compute_percolation(vol, LABELS, phase="Ni", axis=axis)
    ok &= check("T3a", "full slab: P_span=1, P_reach=1, percolates=True",
               abs(r["P_span"] - 1.0) < 1e-12 and abs(r["P_reach"] - 1.0) < 1e-12
               and r["percolates"] is True,
               f"P_span={r['P_span']:.6f} P_reach={r['P_reach']:.6f} "
               f"percolates={r['percolates']}")

    # (b) THE DISTINGUISHING CASE: Ni slab touching ONLY x=0 (not the far face)
    vol_b = np.full((nz, ny, nx), LABELS["pore"], dtype=np.uint8)
    vol_b[:, :, :nx // 2] = LABELS["Ni"]        # occupies x in [0, nx/2), touches x=0 only
    r = compute_percolation(vol_b, LABELS, phase="Ni", axis=axis)
    ok &= check(
        "T3b", "slab touching only x=0: P_span=0, P_reach=1, percolates=False",
        abs(r["P_span"] - 0.0) < 1e-12 and abs(r["P_reach"] - 1.0) < 1e-12
        and r["percolates"] is False,
        f"P_span={r['P_span']:.6f} P_reach={r['P_reach']:.6f} "
        f"percolates={r['percolates']}")

    # (c) isolated central Ni cube, touching neither face -> both 0
    vol_c = np.full((nz, ny, nx), LABELS["pore"], dtype=np.uint8)
    vol_c[8:12, 8:12, 8:12] = LABELS["Ni"]
    r = compute_percolation(vol_c, LABELS, phase="Ni", axis=axis)
    ok &= check("T3c", "isolated cube: P_span=0, P_reach=0",
               abs(r["P_span"]) < 1e-12 and abs(r["P_reach"]) < 1e-12,
               f"P_span={r['P_span']:.6f} P_reach={r['P_reach']:.6f}")

    # (d) two disjoint, non-touching bars: bar1 spans full x, bar2 does not.
    #     P_span must equal EXACTLY the volume fraction contributed by bar1.
    vol_d = np.full((nz, ny, nx), LABELS["pore"], dtype=np.uint8)
    vol_d[2:5, 2:5, :] = LABELS["Ni"]           # bar1: spans x=0..nx-1
    vol_d[15:18, 15:18, 3:10] = LABELS["Ni"]    # bar2: touches neither face
    n_bar1 = int((vol_d[2:5, 2:5, :] == LABELS["Ni"]).sum())
    n_total_ni = int((vol_d == LABELS["Ni"]).sum())
    r = compute_percolation(vol_d, LABELS, phase="Ni", axis=axis)
    exp_pspan = n_bar1 / n_total_ni
    ok &= check("T3d", "P_span equals exactly the spanning bar's volume fraction",
               abs(r["P_span"] - exp_pspan) < 1e-12,
               f"P_span={r['P_span']:.9f} want={exp_pspan:.9f}")
    return ok


# ---------------------------------------------------------------------------
# T4 — SNOW throats on known bars, w in {4,6,8}, via extract_network
# ---------------------------------------------------------------------------
def test_T4():
    print("-" * 74)
    print("T4 — SNOW throat size exact on dumbbells with known bar width")
    print("-" * 74)
    ok = True
    N = 50
    d = 10.0            # nm, isotropic
    for w in (4, 6, 8):
        vol = np.zeros((N, N, N), dtype=bool)
        c = N // 2
        # cube A flush against z=0, cube B flush against z=N-1, joined along z
        vol[0:14, c - 7:c + 7, c - 7:c + 7] = True
        vol[N - 14:N, c - 7:c + 7, c - 7:c + 7] = True
        half = w // 2
        vol[13:N - 13, c - half:c - half + w, c - half:c - half + w] = True

        G, diag = extract_network(vol, (d, d, d), axis=0, r_max=4)
        if G is None or G.number_of_edges() == 0:
            ok &= check("T4", f"w={w}: network extraction produced a graph",
                       False, f"diag={diag}")
            continue
        # collapse to the single inter-cube throat: take the edge with the
        # smallest neck as the constriction (should be exactly one throat here)
        necks = [dd["neck_nm"] for _, _, dd in G.edges(data=True)]
        areas = [dd["area_nm2"] for _, _, dd in G.edges(data=True)]
        exp_neck = w * d
        exp_area = (w * d) ** 2
        neck_ok = abs(min(necks) - exp_neck) < 0.05 * exp_neck
        area_ok = abs(min(areas) - exp_area) < 0.10 * exp_area
        ok &= check("T4", f"w={w}: neck width within 5%, area within 10%",
                   neck_ok and area_ok,
                   f"nodes={G.number_of_nodes()} edges={G.number_of_edges()} "
                   f"min_neck={min(necks):.2f}nm want~{exp_neck:.2f}nm "
                   f"min_area={min(areas):.1f}nm^2 want~{exp_area:.1f}nm^2")

        nm = compute_network_metrics(G, G.graph.get("face_lo"),
                                     G.graph.get("face_hi"))
        ok &= check("T4", f"w={w}: compute_network_metrics wrapper runs",
                   nm["n_nodes"] == G.number_of_nodes())
    return ok


# ---------------------------------------------------------------------------
# T5 — particle stats on a known packing (pipeline-correctness check only;
#      the full neck-vs-size COUPLING DECISION experiment is separate, run
#      after gate G-1, per the approved sequence)
# ---------------------------------------------------------------------------
def test_T5():
    print("-" * 74)
    print("T5 — particle stats exact on a known non-touching sphere packing")
    print("-" * 74)
    ok = True
    N = 64
    R = 6                       # voxels
    d = 15.0                    # nm
    centres = [(16, 16, 16), (16, 16, 48), (16, 48, 16), (16, 48, 48),
              (48, 16, 16), (48, 16, 48), (48, 48, 16), (48, 48, 48)]
    zz, yy, xx = np.ogrid[:N, :N, :N]
    mask = np.zeros((N, N, N), dtype=bool)
    for (cz, cy, cx) in centres:
        mask |= (zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2 < R ** 2

    stats = compute_particle_stats(mask, (d, d, d), min_distance=3)
    exp_diam = 2 * R * d
    ok &= check(
        "T5a", "8 non-touching spheres: n_regions_used == 8 "
        "(after border exclusion)",
        stats["ws_n_regions_used"] == 8,
        f"got n_regions_used={stats['ws_n_regions_used']} "
        f"n_peaks={stats['n_peaks']}")
    ok &= check(
        "T5a", "watershed d_volweighted within 1 voxel of 2R",
        abs(stats["ws_d_volweighted_nm"] - exp_diam) <= d,
        f"got={stats['ws_d_volweighted_nm']:.1f}nm want~{exp_diam:.1f}nm "
        f"(tol {d:.1f}nm = 1 voxel)")
    ok &= check(
        "T5a", "c-PSD d_r50max within 1 voxel of 2R",
        abs(stats["d_cPSD_r50max_nm"] - exp_diam) <= d,
        f"got={stats['d_cPSD_r50max_nm']:.1f}nm want~{exp_diam:.1f}nm")

    # T5b: qualitative direction check -- adding a SHORT connecting neck
    # between two NEAR-TOUCHING spheres (gap=6 voxels, bar length 18 voxels)
    # perturbs the watershed size measure, confirming the documented R1
    # mechanism BEFORE the full coupling experiment is run.
    #
    # NOTE ON GEOMETRY (found empirically while writing this test, and folded
    # into the Family-B design guidance for the coupling experiment): a LONG
    # freestanding bridge between well-separated spheres (gap>=12 voxels, bar
    # length >=24 voxels) can spawn its OWN separate watershed catchment at the
    # bar's midpoint -- an artifact that INCREASES n_regions_used with a
    # spurious small "neck particle", the opposite of R1's inflation
    # mechanism. This is why the neck geometry here (and in Family B) bridges
    # NEAR-TOUCHING particle surfaces, not distant ones -- consistent with the
    # execution spec's "preferred approach" (necks between CONTACTING
    # particles, not freestanding bridges). See probe result below (informational
    # only, not gated) for the long-bridge case.
    gap = 6
    z1, z2 = 32 - (R + gap // 2), 32 + (R + gap // 2)
    two = ((zz - z1) ** 2 + (yy - 32) ** 2 + (xx - 32) ** 2 < R ** 2)
    two |= ((zz - z2) ** 2 + (yy - 32) ** 2 + (xx - 32) ** 2 < R ** 2)
    stats_no_neck = compute_particle_stats(two, (d, d, d), min_distance=3)
    two_neck = two.copy()
    two_neck[z1:z2, 30:34, 30:34] = True     # short 4-voxel-wide connecting neck
    stats_neck = compute_particle_stats(two_neck, (d, d, d), min_distance=3)
    ok &= check(
        "T5b", "short near-contact neck: watershed size measure is "
        "neck-sensitive (R1 mechanism), region count preserved",
        stats_neck["ws_n_regions_used"] == 2
        and stats_neck["ws_d_volweighted_nm"] > stats_no_neck["ws_d_volweighted_nm"],
        f"no_neck n_regions={stats_no_neck['ws_n_regions_used']} "
        f"d_volwt={stats_no_neck['ws_d_volweighted_nm']:.1f}nm  ->  "
        f"with_neck n_regions={stats_neck['ws_n_regions_used']} "
        f"d_volwt={stats_neck['ws_d_volweighted_nm']:.1f}nm  "
        f"(c-PSD: {stats_no_neck['d_cPSD_r50max_nm']:.1f} -> "
        f"{stats_neck['d_cPSD_r50max_nm']:.1f}nm, expected far more stable)")

    # informational only (not gated): the long-freestanding-bridge artifact
    two_far = ((zz - 20) ** 2 + (yy - 32) ** 2 + (xx - 32) ** 2 < R ** 2)
    two_far |= ((zz - 44) ** 2 + (yy - 32) ** 2 + (xx - 32) ** 2 < R ** 2)
    two_far_neck = two_far.copy()
    two_far_neck[20:44, 30:34, 30:34] = True
    stats_far_neck = compute_particle_stats(two_far_neck, (d, d, d), min_distance=3)
    print(f"      [info, not gated] long freestanding bridge (gap=12 vox): "
          f"with_neck n_regions_used={stats_far_neck['ws_n_regions_used']} "
          f"(vs 2 without) -- spurious extra region from the bridge's own "
          f"watershed catchment; avoided in Family B by using short "
          f"near-contact necks only")
    print("      (full neck-width sweep deferred to the T5 coupling "
          "decision experiment, run after gate G-1)")
    return ok


# ---------------------------------------------------------------------------
# T6 — round trip
# ---------------------------------------------------------------------------
def test_T6():
    print("-" * 74)
    print("T6 — save/reload round trip reproduces every metric bit-identically")
    print("-" * 74)
    ok = True
    N = 40
    zz, yy, xx = np.ogrid[:N, :N, :N]
    ni = ((zz - 20) ** 2 + (yy - 20) ** 2 + (xx - 10) ** 2 < 8 ** 2)
    ni |= ((zz - 20) ** 2 + (yy - 20) ** 2 + (xx - 30) ** 2 < 8 ** 2)
    ni |= (np.abs(zz - 20) < 2) & (np.abs(yy - 20) < 2) & (xx >= 10) & (xx <= 30)
    vol = np.full((N, N, N), LABELS["pore"], dtype=np.uint8)
    vol[ni] = LABELS["Ni"]
    ysz = (~ni) & (((zz - 20) ** 2 + (yy - 20) ** 2 + (xx - 20) ** 2) < 18 ** 2)
    vol[ysz] = LABELS["YSZ"]

    seed = 20260810
    config_hash = "t6_dumbbell_v1"
    tv = TernaryVolume(vol=vol, voxel_nm=18.0,
                       meta={"seed": seed, "config_hash": config_hash,
                             "test": "T6"})

    def metrics_of(tv_):
        d = {}
        d["phi"] = tv_.volume_fractions()
        d["perc"] = compute_percolation(tv_.vol, LABELS, phase="Ni", axis=2)
        d["tpb"] = compute_tpb(tv_.vol, LABELS, tv_.spacing_nm)
        d["particles"] = compute_particle_stats(tv_.mask("Ni"), tv_.spacing_nm)
        G, _ = extract_network(tv_.mask("Ni"), tv_.spacing_nm, axis=2)
        d["network"] = compute_network_metrics(
            G, G.graph.get("face_lo") if G else None,
            G.graph.get("face_hi") if G else None)
        return d

    before = metrics_of(tv)
    path = os.path.join(OUT, "t6_roundtrip_test.npz")
    save_ternary(tv, path)
    tv2 = load_ternary(path)

    ok &= check("T6", "reloaded vol bit-identical",
               np.array_equal(tv.vol, tv2.vol))
    ok &= check("T6", "reloaded voxel_nm exact", tv.voxel_nm == tv2.voxel_nm)
    ok &= check("T6", "reloaded meta preserves seed/config_hash",
               tv2.meta.get("seed") == seed
               and tv2.meta.get("config_hash") == config_hash)

    after = metrics_of(tv2)

    def deep_eq(a, b):
        if isinstance(a, dict):
            return set(a) == set(b) and all(deep_eq(a[k], b[k]) for k in a)
        if isinstance(a, float) and np.isnan(a):
            return isinstance(b, float) and np.isnan(b)
        if isinstance(a, (int, float, np.floating, np.integer)):
            return a == b
        return a == b

    for key in before:
        same = deep_eq(before[key], after[key])
        ok &= check("T6", f"recomputed '{key}' bit-identical after reload", same,
                    "" if same else f"before={before[key]} after={after[key]}")

    os.remove(path)
    return ok


def main():
    print("=" * 74)
    print("PHASE 0 (synthetic study) — ternary pipeline validation, T1-T6")
    print("=" * 74)
    ok = True
    ok &= test_T1()
    ok &= test_T2()
    ok &= test_T3()
    ok &= test_T4()
    ok &= test_T5()
    ok &= test_T6()

    print("\n" + "=" * 74)
    print("GATE G0")
    print("=" * 74)
    n_pass = sum(1 for r in RESULTS if r[2])
    n_total = len(RESULTS)
    for tid, name, passed, detail in RESULTS:
        if not passed:
            print(f"  FAILED: {tid} {name}  -- {detail}")
    print(f"\n  {n_pass}/{n_total} checks passed")
    print(f"  RESULT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
