# arXiv submission package

## Files

- `arxiv-submission.tar.gz` — upload this. Contains `manuscript.tex` and
  `figs/` (nine PDFs). Verified by a clean-room build from the tarball contents
  alone: 0 errors, 0 undefined references, 16 pages (rebuilt 2026-11-13 after
  commit b3cc705).
- `arxiv_abstract.txt` — the abstract for the web form, plain text, 1907
  characters (arXiv's limit is 1920).
- `arxiv_abstract.py` — regenerates a plain-text abstract from the LaTeX source
  if the abstract changes.

No `.bbl` is needed: the bibliography is an inline `thebibliography`
environment. `elsarticle.cls` is part of TeX Live and is not bundled.

## Metadata to enter

**Title.** Validation criteria that exclude their own mechanism: six failure
modes of voxel damage operators applied to Ni–YSZ tomograms

**Abstract.** Paste `arxiv_abstract.txt` without its last line (the character
count). It is generated from the manuscript's own abstract, so the two match.

**Primary category.** `cond-mat.mtrl-sci` (Materials Science).
**Cross-list.** `physics.comp-ph` (Computational Physics). Both fit: the
subject is a materials system, the contribution is about simulation method.

**Comments field.** Suggested: "16 pages, 9 figures. All code, pre-registration
documents and result tables at
https://github.com/mehtaarnav/voxel-damage-operators-niysz"

**License.** CC BY 4.0 unless you have a reason to choose otherwise. It is the
most permissive standard option and does not restrict later journal
submission. arXiv's default (non-exclusive licence to distribute) is more
restrictive for reuse.

**MSC/ACM.** Leave blank.

## Blockers to clear before submitting

Status as of 2026-11-13. The repository is public and the manuscript gives its
full URL and the Zenodo concept DOI, so the two earlier blockers on
availability are cleared.

1. **Push before posting.** Commits made after the last push (including
   b3cc705, which adds the swelling paragraph and `out/relook/subblock_scan.csv`
   that it cites) are not on GitHub. The availability statement points readers
   there, so push first.
2. **Archive a new Zenodo version.** The concept DOI resolves to the latest
   archived release. Cut a GitHub release after pushing so the archive contains
   the tables the corrected text cites.
3. **Journal intent.** If the target is *Journal of Power Sources*, posting a
   preprint is permitted under Elsevier policy, but confirm which version may be
   posted and add the DOI later if it is accepted.

## Checks done before this package was built

- `scripts/project2/verify_manuscript_numbers.py`: 60/60 checks pass.
- The threshold sweep (nine thresholds, Fig. 6) and the 26-connectivity
  sensitivity check are in the manuscript; both were listed here as open in an
  earlier version of this file.
