"""
Ternary voxel volume container for the synthetic decoupling study.

This is deliberately minimal: a labelled uint8 array plus a physical voxel
size, and one helper to read off volume fractions. It does NOT duplicate
`cmlib.phases`, which resolves labels from the real dataset's own metadata
(brightness ordering) -- that machinery is specific to the six Zenodo stacks
and is not relevant to a synthetic generator, which chooses its own label
values directly.

LABEL CONVENTION (for every array this module or `cmlib.synth`/`cmlib.damage`
produce): 0 = pore, 1 = YSZ, 2 = Ni. Fixed and explicit, unlike the real
dataset where every stack used a different encoding (see cmlib/phases.py).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

LABEL_PORE = 0
LABEL_YSZ = 1
LABEL_NI = 2
LABELS = {"pore": LABEL_PORE, "YSZ": LABEL_YSZ, "Ni": LABEL_NI}


@dataclass
class TernaryVolume:
    """A labelled ternary voxel volume with a physical, ISOTROPIC voxel size.

    Isotropic by construction (see execution-spec Q7 and the porespy unit trap
    documented in cmlib/particles.py): every synthetic-study tool that needs a
    physical length (SNOW, c-PSD) either requires or silently assumes a scalar
    voxel size, so anisotropic synthetic volumes would reintroduce exactly the
    kind of silent unit bug already found and fixed once.
    """

    vol: np.ndarray                  # uint8, values in {0,1,2}
    voxel_nm: float                  # isotropic voxel size, nm
    meta: dict = field(default_factory=dict)   # generator params, seed, etc.

    def __post_init__(self):
        if self.vol.dtype != np.uint8:
            raise ValueError(f"vol must be uint8, got {self.vol.dtype}")
        if self.vol.ndim != 3:
            raise ValueError(f"vol must be 3D, got shape {self.vol.shape}")
        bad = set(np.unique(self.vol).tolist()) - {LABEL_PORE, LABEL_YSZ, LABEL_NI}
        if bad:
            raise ValueError(f"vol contains labels outside {{0,1,2}}: {bad}")
        if self.voxel_nm <= 0:
            raise ValueError(f"voxel_nm must be > 0, got {self.voxel_nm}")

    @property
    def spacing_nm(self) -> tuple[float, float, float]:
        """(dz, dy, dx), all equal -- kept as a 3-tuple for API compatibility
        with cmlib functions that expect anisotropic `spacing_nm`."""
        return (self.voxel_nm, self.voxel_nm, self.voxel_nm)

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(self.vol.shape)

    def mask(self, phase: str) -> np.ndarray:
        return self.vol == LABELS[phase]

    def volume_fractions(self) -> dict:
        return volume_fractions_from_volume(self.vol)


def save_ternary(tv: "TernaryVolume", path: str) -> None:
    """Save a TernaryVolume losslessly (vol + voxel_nm + meta-as-JSON).

    Used for the Phase-0 round-trip test (T6) and, from Phase 1 onward, to
    persist every generated structure so damaged volumes can be regenerated
    from (config_hash, seed) rather than only stored as summary statistics.
    """
    np.savez_compressed(path, vol=tv.vol, voxel_nm=np.array([tv.voxel_nm]),
                        meta_json=np.array(json.dumps(tv.meta)))


def load_ternary(path: str) -> "TernaryVolume":
    with np.load(path, allow_pickle=False) as z:
        vol = z["vol"].astype(np.uint8)
        voxel_nm = float(z["voxel_nm"][0])
        meta = json.loads(str(z["meta_json"]))
    return TernaryVolume(vol=vol, voxel_nm=voxel_nm, meta=meta)


def volume_fractions_from_volume(vol: np.ndarray) -> dict:
    """Exact voxel-count volume fractions of a {0,1,2}-labelled ternary array."""
    total = vol.size
    counts = np.bincount(vol.ravel(), minlength=3)
    return {
        "phi_pore": float(counts[LABEL_PORE] / total),
        "phi_YSZ": float(counts[LABEL_YSZ] / total),
        "phi_Ni": float(counts[LABEL_NI] / total),
    }
