"""Sequential greedy volume-conserving Ni swap, with incremental updates.

Semantics identical to the naive sequential operator in
`scripts/project2/o7_sequential_greedy.py`:

    dA = 2*(nN(a) - nN(b));  accept iff dA <= 0;  a Ni surface site, b a
    pore-front site (pore = not Ni and not YSZ), a and b non-adjacent;
    exactly one voxel out and one in, so dV = 0.

`nN` is the count of Ni **6-neighbours**, CENTRE EXCLUDED, as frozen in
PREREG_O5V2_OPTIONB.md. (cmlib/damage2.py:31 uses a centre-INCLUDED stencil,
which biases Ni sites by +1; see out/project2/O7_O5V2B_RERUN_REPORT.md.)

The naive version recomputes the whole neighbour field after every accepted
move, which is O(moves x volume) and infeasible on a 67 Mvoxel ROI. Here nN is
maintained incrementally: a swap changes nN only for the 6 neighbours of `a` and
the 6 neighbours of `b`. Candidates are held in 7 buckets (nN = 0..6) with lazy
validation, so selecting the extremal candidate is O(1) amortised.

This module must reproduce the naive operator exactly; see
`scripts/project2/o7_seqgreedy_validate.py`.
"""
from __future__ import annotations

import numpy as np
import scipy.ndimage as ndi

ST7 = ndi.generate_binary_structure(3, 1)
ST6 = ST7.copy()
ST6[1, 1, 1] = False


def _nb_field(ni):
    return ndi.convolve(ni.astype(np.int8), ST6.astype(np.int8),
                        mode="constant", cval=0).astype(np.int8)


class SeqGreedy:
    def __init__(self, ni, ysz, seed=0, strict=False):
        # strict=True accepts only dA < 0, forbidding the area-NEUTRAL moves
        # that dominate under dA <= 0. Used to test whether the neutral
        # plateau is the sole driver of TPB manufacture.
        self.strict = bool(strict)
        self.ni = ni.copy()
        self.ysz = ysz
        self.shape = ni.shape
        self.nz, self.ny, self.nx = ni.shape
        self.nb = _nb_field(self.ni)
        self.rng = np.random.default_rng(seed)
        self.proposed = 0
        self.accepted = 0
        self.neutral = 0
        # strides for flat <-> 3d without unravel_index in the hot loop
        self.sy, self.sz = self.nx, self.nx * self.ny
        self._build_buckets()

    # ---------------------------------------------------------------- buckets
    def _is_surf(self, flat):
        """Ni voxel with at least one non-Ni 6-neighbour, strictly interior."""
        if not self.ni.flat[flat]:
            return False
        return self.nb.flat[flat] < 6 and self._interior(flat)

    def _is_front(self, flat):
        """Pore voxel (not Ni, not YSZ) adjacent to Ni, strictly interior."""
        if self.ni.flat[flat] or self.ysz.flat[flat]:
            return False
        return self.nb.flat[flat] >= 1 and self._interior(flat)

    def _interior(self, flat):
        z, rem = divmod(flat, self.sz)
        y, x = divmod(rem, self.sy)
        return (0 < z < self.nz - 1 and 0 < y < self.ny - 1
                and 0 < x < self.nx - 1)

    def _build_buckets(self):
        ni, ysz, nb = self.ni, self.ysz, self.nb
        interior = np.zeros(self.shape, dtype=bool)
        interior[1:-1, 1:-1, 1:-1] = True
        surf = ni & (nb < 6) & interior
        front = (~ni) & (~ysz) & (nb >= 1) & interior
        self.bs = [None] * 7          # surface Ni, by nN
        self.bf = [None] * 7          # pore front, by nN
        for v in range(7):
            self.bs[v] = list(np.flatnonzero((surf & (nb == v)).ravel()))
            self.bf[v] = list(np.flatnonzero((front & (nb == v)).ravel()))
        self.n_surf0 = int(surf.sum())

    def _push(self, flat):
        v = int(self.nb.flat[flat])
        if self._is_surf(flat):
            self.bs[v].append(flat)
        if self._is_front(flat):
            self.bf[v].append(flat)

    # ------------------------------------------------------ tie-breaking
    # FROZEN POLICY (do not change without re-freezing the pre-registration):
    #   Among the valid candidates in the extremal occupancy bucket -- i.e.
    #   among moves with IDENTICAL dA -- one is chosen UNIFORMLY AT RANDOM from
    #   the seeded generator. Ties are broken by chance, never by insertion
    #   order.
    #
    # This is load-bearing, not cosmetic. Area-neutral (dA = 0) moves dominate
    # on real ROIs, so the achieved area reduction is set by which of the many
    # equal-dA moves is taken. A LIFO policy reached dS = -0.00024 on the
    # recovered O5v2 structure where a randomised one reached -0.02817. The
    # policy is therefore part of the operator specification and is frozen here.
    def _sample(self, bucket, v, validator):
        """Uniform random valid entry of `bucket`, purging stale ones."""
        while bucket:
            i = int(self.rng.integers(len(bucket)))
            f = bucket[i]
            if validator(f) and int(self.nb.flat[f]) == v:
                return f
            bucket[i] = bucket[-1]        # O(1) swap-remove
            bucket.pop()
        return None

    def _pop_min_surf(self):
        """Lowest-nN valid surface site, ties broken uniformly at random."""
        for v in range(7):
            f = self._sample(self.bs[v], v, self._is_surf)
            if f is not None:
                return f, v
        return None, None

    def _pop_max_front(self):
        """Highest-nN valid front site, ties broken uniformly at random."""
        for v in range(6, -1, -1):
            f = self._sample(self.bf[v], v, self._is_front)
            if f is not None:
                return f, v
        return None, None

    def _neighbours(self, flat):
        z, rem = divmod(flat, self.sz)
        y, x = divmod(rem, self.sy)
        out = []
        if z > 0:
            out.append(flat - self.sz)
        if z < self.nz - 1:
            out.append(flat + self.sz)
        if y > 0:
            out.append(flat - self.sy)
        if y < self.ny - 1:
            out.append(flat + self.sy)
        if x > 0:
            out.append(flat - 1)
        if x < self.nx - 1:
            out.append(flat + 1)
        return out

    # ------------------------------------------------------------------- move
    def step(self):
        a, va = self._pop_min_surf()
        if a is None:
            return False
        # non-adjacency: the dA algebra assumes a and b share no bond.
        # Resample rather than removing the candidate -- under the frozen
        # random tie-breaking the sampled entry is not at any fixed position,
        # so it cannot be popped by index. Adjacency collisions are rare
        # because the buckets are large; give up after a bounded number of
        # attempts so the loop cannot spin.
        nbrs = self._neighbours(a)
        b = vb = None
        for _ in range(8):
            cand, vcand = self._pop_max_front()
            if cand is None:
                return False
            if cand not in nbrs:
                b, vb = cand, vcand
                break
        if b is None:
            return False
        self.proposed += 1
        dA = 2 * (va - vb)
        if dA > 0 or (self.strict and dA == 0):
            return False              # greedy: extremal pair fails => none pass
        # apply
        self.ni.flat[a] = False
        self.ni.flat[b] = True
        for nbr in self._neighbours(a):
            self.nb.flat[nbr] -= 1
        for nbr in self._neighbours(b):
            self.nb.flat[nbr] += 1
        self.accepted += 1
        if dA == 0:
            self.neutral += 1
        # refresh membership of everything whose state or nN changed
        touched = {a, b}
        touched.update(self._neighbours(a))
        touched.update(self._neighbours(b))
        for f in touched:
            self._push(f)
        return True

    def run(self, n_moves):
        for _ in range(int(n_moves)):
            if not self.step():
                break
        return self.ni
