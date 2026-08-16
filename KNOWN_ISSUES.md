# Known issues

Defects found by review that are recorded rather than fixed, with the reason.
Everything that could change a number in the manuscript has been fixed; what
remains is either in code whose results were withdrawn, or is library hardening
for code that has a single caller. Changing those before submission would add
risk without changing a result.

## Fixed

- `SeqGreedy` bucket entries could duplicate, so tie-breaking was weighted by
  how often a voxel had been touched rather than uniform over valid candidates.
  Entries now carry a version stamp; at most one is live per voxel.
- `SeqGreedy.step` ended the whole run when the extremal source had no
  non-adjacent partner. It now retires that source and tries the next.
- `ni_ysz_interface_area_vox` counted faces across opposite domain boundaries
  via `np.roll`. Now uses free-boundary slicing.
- A comment in `apply_o5v2b` asserted the round could not raise area. It can,
  and does; that is the batching result the manuscript reports.
- `_mean_curvature` was unused and used the centre-included stencil. Removed.

## Recorded, not fixed

- **Masks are assumed boolean.** `~mask` is bitwise on integer arrays, so an
  integer mask would silently corrupt phase logic. Every caller in this
  repository builds masks with `a == label`, which is boolean, so no result is
  affected. A library used elsewhere should cast at its entry points.
- **`apply_o5v2b` ignores its `seed`.** Documented in the function. It is
  retained only to reproduce the published batched figure; production runs use
  `cmlib.seqgreedy`.
- **`apply_o5` and `apply_o5v2` are batched within a round.** Intentional for
  `apply_o5`; for `apply_o5v2` it confounds the curvature-rank demonstration,
  which is now stated in the manuscript's limitations.
- **`o2_candidates` uses a strict `<` at the percentile**, so ties can shrink
  the candidate set below the nominal fraction. `apply_o2` also counts a throat
  as severed when its region slices are missing. Both are in the synthetic
  platform whose results were withdrawn.
- **`apply_o5v2(conn26=True)`** switches the curvature stencil but leaves the
  surface and front definitions at 6-connectivity.
- **`apply_o5` uses `uniform_filter`'s default reflective boundary**, which
  differs from the constant-zero convention used elsewhere.
- **Boundary conventions differ across operators**, as now stated in the
  manuscript.
