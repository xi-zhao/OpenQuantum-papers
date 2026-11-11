# QwenQPU

**QwenQPU: A Language Model with an Internal Quantum Circuit Simulator**

Author: **OpenQuantum**. Reference-checked manuscript dated 11 November 2026; 17 pages, 5 figures, 5 main tables and 2 appendix tables. Prepared for arXiv; not submitted.

## Files

- [Manuscript PDF](paper/main.pdf)
- [LaTeX source](paper/main.tex), with the inline bibliography and [five figure PDFs](paper/figures)
- [arXiv source archive](QwenQPU-arXiv-source.zip), containing only `main.tex` and the five figure PDFs
- [Reference audit](reference-audit.tsv), covering all 30 references and the limits of source access
- [Validation record](validation.json)
- [SHA-256 file manifest](SHA256SUMS)

## Build

Run from `paper/` in an environment with pdfLaTeX and the standard LaTeX packages listed in `main.tex`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The delivered PDF was compiled with TeX Live 2026. All fonts are embedded; all 17 pages were visually inspected. The manuscript retains its stated finite-capacity and finite-precision limits, and does not claim physical quantum speedup. The fixed experiment was not rerun as part of this manuscript transfer.

## Evidence and provenance

The manuscript describes the experiment at [QwenQPU code commit caae45327b46458163d0840a4c1466e5575a7b4f](https://github.com/xi-zhao/QwenQPU/tree/caae45327b46458163d0840a4c1466e5575a7b4f). Model identity, experimental records and reconstruction details are documented in the manuscript and its appendices.

The citation audit checked reference identities and claim-to-source alignment, corrected two body citations and enriched nine bibliography entries. Wilson intervals, exact conditional McNemar probabilities and conservative paired-difference intervals were independently recalculated. Some publisher full texts could not be retrieved; the audit distinguishes metadata, author versions and method-documentation checks from full-text inspection.

This branch stores the manuscript source and supporting records. The `main` branch receives only the reviewed PDF and the repository index, following the repository's existing convention. No new publication license or arXiv submission is created by this transfer.
