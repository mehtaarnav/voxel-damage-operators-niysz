# arXiv submission package

## Files

- `arxiv-submission.tar.gz` — upload this. Contains `manuscript.tex` and
  `figs/` (six PDFs). Verified by a clean-room build from the tarball contents
  alone: 0 errors, 0 undefined references, 14 pages.
- `arxiv_abstract.txt` — the abstract for the web form, plain text, 1905
  characters (arXiv's limit is 1920).
- `arxiv_abstract.py` — regenerates a plain-text abstract from the LaTeX source
  if the abstract changes.

No `.bbl` is needed: the bibliography is an inline `thebibliography`
environment. `elsarticle.cls` is part of TeX Live and is not bundled.

## Metadata to enter

**Title.** Validation criteria that exclude their own mechanism: six failure
modes of voxel damage operators applied to Ni–YSZ tomograms

**Abstract.** Paste `arxiv_abstract.txt`. Note this is a shortened form — the
manuscript's own abstract is 2664 characters and does not fit the field. If you
want them identical, adopt the shortened text in the paper too.

**Primary category.** `cond-mat.mtrl-sci` (Materials Science).
**Cross-list.** `physics.comp-ph` (Computational Physics). Both fit: the
subject is a materials system, the contribution is about simulation method.

**Comments field.** Suggested: "14 pages, 6 figures. All code, pre-registration
documents and result tables at <repository URL>."

**License.** CC BY 4.0 unless you have a reason to choose otherwise. It is the
most permissive standard option and does not restrict later journal
submission. arXiv's default (non-exclusive licence to distribute) is more
restrictive for reuse.

**MSC/ACM.** Leave blank.

## Blockers to clear before submitting

1. **The repository is private.** The data-availability section names
   `voxel-damage-operators-niysz` and lists specific CSV paths. An arXiv posting
   is public immediately; a reader who cannot reach the repository will treat
   the availability claim as unmet. Either make it public before posting, or
   change the statement to availability on request.
2. **No repository URL appears in the manuscript.** Once it is public, add the
   full GitHub URL to the availability section — the bare repository name is not
   resolvable by a stranger.
3. **Journal intent.** If the target is *Journal of Power Sources*, posting a
   preprint is permitted under Elsevier policy, but confirm which version may be
   posted and add the DOI later if it is accepted.

## Known open items, not blockers

- Threshold sensitivity for the erosion bisection is reported at 0.50 and 0.10
  only; a sweep over intermediate thresholds would strengthen the sign-reversal
  claim and has not been run.
- A 26-connectivity sensitivity check was not performed; the limitation is now
  stated in the text.
