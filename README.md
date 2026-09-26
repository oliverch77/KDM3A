# KDM3A analysis scripts

R code for bulk differential expression, volcano plots, and single-cell signature scoring accompanying the KDM3A manuscript accepted by *Science Advances*.

## Scripts

- `01_bulk_DE_volcano.R`: human MDM, mouse BMM, and RAW264.7; KD/KO versus control only, separately without and with LPS.
- `02_single_cell_scores.R`: two top-100 fold-change signatures, module scores, patient-adjusted mixed models, and plots.

Required libraries are loaded explicitly at the top of each script. Use the original package versions; packages are not installed or updated automatically. Each run records `sessionInfo.txt`.

## Inputs

1. Original tab-delimited **RefSeq transcript counts**, with sample headers and column order preserved. HOMER annotation columns are accepted.
2. A prepared **Seurat RDS** with normalized expression, existing cluster identities, `Patient`, and `percent.mt` metadata. Object preparation is outside this repository.

## Repository files

`CITATION.cff` supplies citation metadata; `LICENSE` contains a proposed MIT code license. `.gitignore` excludes inputs/results, `.gitattributes` standardizes text files, and the GitHub workflow checks R syntax without running analyses.

