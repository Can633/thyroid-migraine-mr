# Thyroid disease and migraine: two-sample Mendelian randomization and colocalization

Analysis code and exported result tables for the manuscript

> "Shared polygenic architecture and limited evidence for a causal effect of thyroid liability on migraine: a two-sample Mendelian randomization and colocalization study"

This repository contains the analysis scripts and the complete exported result tables needed to reproduce the reported analyses. It accompanies the submission to *Headache: The Journal of Head and Face Pain*.

## Zenodo archive

A permanent, citable snapshot of this repository is archived on Zenodo:

- **DOI:** `[STILL MISSING — paste the Zenodo DOI here after publishing, e.g. 10.5281/zenodo.12345678]`

## Repository contents

| Folder | Contents |
|---|---|
| `code/` | R (`.R`), Python (`.py`) and shell (`.sh`) analysis scripts |
| `results/` | Main analysis outputs: MR estimates, heterogeneity, pleiotropy, leave-one-out, radial MR, MRlap, LDSC genetic correlation, colocalization, eQTL annotation, power calculations, MHC sensitivity analyses, and harmonised SNP-level data |
| `results_ukb/` | UK Biobank external-replication outputs (MR and LDSC) |

## Software environment

- **R 4.2.0** — packages `TwoSampleMR` 0.7.13, `coloc` 5.2.3, `MRlap` 0.0.3.3, `ieugwasr` 1.2.0
- **Python 3.13** — `pandas` 2.2, `NumPy` 1.26, `SciPy` 1.11

## Data sources

All input data are publicly available summary statistics:

- Thyroid phenotypes: Global Biobank Meta-analysis Initiative (GWAS Catalog GCST90627740, GCST90627744, GCST90627749, GCST90627754, GCST90627763, GCST90627765)
- Migraine: FinnGen data freeze 13 (endpoints G6_MIGRAINE, G6_MIGRAINE_WITH_AURA, G6_MIGRAINE_NO_AURA; https://r13.finngen.fi)
- Replication: UK Biobank whole-genome-sequencing migraine GWAS (GCST90473326)
- eQTL annotation: GTEx Consortium version 8

## Reproducibility

Reporting follows the STROBE-MR statement; the completed checklist is provided as supplementary material with the manuscript. Scripts are grouped by analysis stage: instrument extraction (`extract_ivs*`, `rev0*`), outcome SNP extraction and harmonisation (`*extract_outcome*`, `extract_small_outcome.py`; harmonised SNP-level files are in `results/`), core MR (`run_mr*`, `make_report.R`), overlap correction (`overlap_correct.py`, `mrlap_*`, `ldsc_overlap.py`), genetic correlation (`ldsc_rg*`), colocalization (`coloc_*`, `run_coloc.R`), eQTL annotation (`eqtl_*`), MHC sensitivity (`mhc0*`), and power analyses (`power*`).

## Citation

If you use these code or results, please cite the associated manuscript and its Zenodo archive (DOI above).
