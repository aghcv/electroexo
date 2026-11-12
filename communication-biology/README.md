# Communications Biology draft

This directory contains a journal-specific rewrite of the electro-exocytosis manuscript. The source manuscript in `/Users/aghorban/code/electro-exocytosis` was not modified.

## Structure

- `main.tex`: submission manuscript in the journal-preferred order: Introduction, Results, Discussion, Methods.
- `content/`: section-level manuscript source.
- `figures/main/`: five main display items.
- `figures/supplementary/`: layer diagnostics, calibration checks, sensitivity analyses and the extended dosimetry figure.
- `supplementary/supplementary_information.tex`: self-contained Supplementary Information source.
- `supplementary/model_specification_appendix.tex`: the former extended appendix tables and equation scaffold.
- `supplementary/original_sage_draft_preserved.tex`: an unmodified preservation copy of the source draft.
- `references/references.bib`: the source bibliography plus journal-neighbor papers used in the rewrite.
- `source-data/`: location for figure-level machine-readable data or a repository manifest.

## Editorial position

The rewrite presents the work as an implemented model with a preliminary, size-resolved observation fit and explicit validation limits. It avoids describing provisional equations as established biology. The main claim is that separating intracellular release, extracellular kinetics and measurement exposes missing dose dependence and experimental priorities.

## Required author actions before submission

1. Add the corresponding-author email and confirm author names, affiliations and order.
2. Confirm experimental metadata identified in `SUBMISSION_GAPS.md`.
3. Add funding and author-contribution statements.
4. Create a versioned public code release and persistent data archive.
5. Regenerate figures using the specifications in `FIGURE_AUDIT.md`.
6. Compile the article and Supplementary Information PDFs with a LaTeX distribution; none is installed in the present environment.
7. Complete the Nature Portfolio Reporting Summary and editorial-policy checklist if requested by the journal.
