# Genomic and transcriptomic landscape of cutaneous squamous cell carcinoma

## Overview 
This repository contains the code to reproduce the analyses in the publication: Comprehensive molecular characterization of cutaneous squamous cell carcinoma reveals determinants of metastatic progression.

## Requirements

- R 4.4.0
- R packages are loaded in WES_01_main.R

## Installation

```bash
git clone https://github.com/barentroia/onco_cSCC_genomic_characterization.git
cd onco_cSCC_genomic_characterization
```

## Regenerate publication figures

```r
# Run WES analyses:
source("WES_01_main.R")
# Run RNA-seq and SCCore-GEP analyses:
source("RNASeq_01_main.R")
```
## Applying SCCore-GEP to an external cohort

To run the SCCore-GEP model on the Nassir cohort, or on another external cohort, use the Quarto tutorial:

```bash
quarto render tutorial/Tutorial.qmd


You can also open `tutorial/Tutorial.qmd` in RStudio and click **Render**.
To use your own cohort, edit the input file and cohort name at the top of the notebook.
```

## Structure and scripts description
```
|
└── src                                  <- Folder with code to reproduce the analyses in the manuscript and regenerate figures and tables
    |
    ├── clinical/                        <- Folder with scripts to run multiple imputation in D-SQUAME discovery and validation datasets and calculate inverse sampling probability weights.
    ├── functions                        <- Folder with functions re-used in many scripts 
    ├── integration                      <- Folder with scripts to run analysis on samples with both GEP and WES, for SFigure 11
    ├── sccore_gep                       <- Folder with scripts to run analyses and generate figures/tables related to SCCore-GEP development.
    ├── single_cell_spatial              <- Folder with scripts to run spatial transcriptomics analyses
    ├── rna_seq                          <- Folder with analyses on rna-seq data (clustering and differential gene expression)
    └── validation/                      <- Folder with scripts to run analyses and generate figures/tables for related to SCCore-GEP validation analyses.
    ├── RNASeq_00_config.R               <- Script with configurations and filenames to run all analyses on RNA-seq data
    ├── WES_00_required_libraries.R      <- Script with configurations and filenames to run all analyses on WES data
    ├── RNASeq_01_main.R                 <- Main script to run all analyses on RNA-seq data
    ├── WES_01_main.R                    <- Main script to run all analyses on WES data

```
## Data Availability

The raw sequencing data are available at EGAS50000001852.

The code in this repository reproduces the analyses described in the manuscript.

## Contact

Name: Barbara Pacheco

Email: b.rentroiapacheco.at.erasmusmc.nl

Institution: Erasmus Medical Center, Rotterdam, Netherlands / Shain Lab at UCSF, CA, USA
