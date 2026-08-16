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
    def __init__(self, ni, ysz, seed=0, strict=False, tiebreak="random"):
        # strict=True accepts only dA < 0, forbidding the area-NEUTRAL moves
        # that dominate under dA <= 0. Used to test whether the neutral
        # plateau is the sole driver of TPB manufacture.
        self.strict = bool(strict)
        # tiebreak: "random" is the FROZEN policy. "lifo" is retained only so
        # the reported sensitivity of the result to this choice stays
        # reproducible; it is not a valid setting for a production run.
        assert tiebreak in ("random", "lifo")
        self.tiebreak = tiebreak
        self.ni = ni.copy()
        self.ysz = ysz
        self.shape = ni.shape
        self.nz, self.ny, self.nx = ni.shape
        self.nb = _nb_field(self.ni)
        self.rng = np.random.default_rng(seed)
        self.proposed = 0
        self.accepted = 0
        self.neutral = 0
        self.dA_log = None   # set to [] to record dA of every accepted move
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
        self.ver = np.zeros(ni.size, dtype=np.int32)
        for v in range(7):
            self.bs[v] = [(int(f), 0) for f in
                          np.flatnonzero((surf & (nb == v)).ravel())]
            self.bf[v] = [(int(f), 0) for f in
                          np.flatnonzero((front & (nb == v)).ravel())]
        self.n_surf0 = int(surf.sum())

    def _push(self, flat):
        """Register a voxel in its current bucket, superseding any earlier entry.

        Entries carry a version stamp. Pushing bumps the voxel's version, which
        makes every earlier entry for that voxel stale wherever it sits, so at
        most one live entry per voxel exists at any time.

        Without this, a voxel whose neighbour count drifts away and back gains a
        second entry that is still valid, and _sample -- which draws uniformly
        over ENTRIES -- would then favour it in proportion to how often it had
        been touched. That is not the frozen policy, which is uniform over
        valid candidates.
        """
        self.ver[flat] += 1
        stamp = self.ver[flat]
        v = int(self.nb.flat[flat])
        if self._is_surf(flat):
            self.bs[v].append((flat, stamp))
        if self._is_front(flat):
            self.bf[v].append((flat, stamp))

    # ------------------------------------------------------ tie-breaking
    # FROZEN POLICY (do not change without re-freezing the pre-registration):
    #   Among the valid candidates in the extremal occupancy bucket -- i.e.
    #   among moves with IDENTICAL dA -- one is chosen UNIFORMLY AT RANDOM from
    #   the seeded generator. Ties are broken by chance, never by insertion
    #   order.
    #
    # This is load-bearing, not cosmetic. Area-neutral (dA = 0) moves dominate
    # on real ROIs, so the achieved area reduction is set by which of the many
    # equal-dA moves is taken. Measured on the recovered O5v2 structure at a
    # matched budget of 1220 moves, LIFO and random tie-breaking give
    # materially different area reductions; see
    # scripts/project2/o7_tiebreak_sensitivity.py, which reproduces both. The
    # policy is therefore part of the operator specification and is frozen
    # here as "random".
    def _sample(self, bucket, v, validator):
        """Valid entry of `bucket` with its index, purging stale ones.

        Returns (flat_index, position_in_bucket) so the caller can remove the
        entry in O(1) if it has to reject it for a reason other than dA --
        necessary because under a deterministic policy re-sampling would
        otherwise return the same candidate forever.
        """
        while bucket:
            i = (len(bucket) - 1 if self.tiebreak == "lifo"
                 else int(self.rng.integers(len(bucket))))
            f, stamp = bucket[i]
            if (stamp == self.ver[f] and validator(f)
                    and int(self.nb.flat[f]) == v):
                return f, i
            bucket[i] = bucket[-1]        # O(1) swap-remove
            bucket.pop()
        return None, None

    def _pop_min_surf(self):
        """Lowest-nN valid surface site; ties broken by the frozen policy."""
        for v in range(7):
            f, i = self._sample(self.bs[v], v, self._is_surf)
            if f is not None:
                return f, v, i
        return None, None, None

    def _pop_max_front(self):
        """Highest-nN valid front site; ties broken by the frozen policy."""
        for v in range(6, -1, -1):
            f, i = self._sample(self.bf[v], v, self._is_front)
            if f is not None:
                return f, v, i
        return None, None, None

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
    def step(self, _tries=64):
        """One move. Returns False only when no admissible move remains.

        If the extremal source has no non-adjacent partner, that source is
        retired and the next is tried. Returning False there would end the whole
        run while valid moves elsewhere were still available.
        """
        for _ in range(_tries):
            moved = self._try_one()
            if moved is not None:
                return moved
        return False

    def _retire(self, flat):
        """Drop a voxel from consideration until something touches it again."""
        self.ver[flat] += 1

    def _try_one(self):
        """True/False if a move was decided, None if this source was retired."""
        a, va, _ = self._pop_min_surf()
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
        held = []                      # adjacent candidates, restored below
        while True:
            cand, vcand, idx = self._pop_max_front()
            if cand is None:
                break
            if cand not in nbrs:
                b, vb = cand, vcand
                break
            bucket = self.bf[vcand]    # remove so the next draw differs
            bucket[idx] = bucket[-1]
            bucket.pop()
            held.append((cand, vcand, self.ver[cand]))
            if len(held) > 6:          # a has at most six 6-neighbours
                break
        for f, v, stamp in held:
            self.bf[v].append((f, stamp))
        if b is None:
            # This source cannot be paired. Retire it and let the caller try
            # the next one rather than ending the run.
            self._retire(a)
            return None
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
        if self.dA_log is not None:
            self.dA_log.append(dA)
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
