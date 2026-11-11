#===============================================================================
# Aim:      Generate Figure 4D: GSEA (MSigDB Hallmarks) of cases versus controls.
# Input:    - DGEA results ~Biopsy_excision+Metastasis (results_dir/dgea)
# Output:   - GSEA bar plot (pdf) and fgsea results (csv) in results_dir/publication
# Requires: code_dir, results_dir (main.R); DGEA_discovery.R run first
# Usage:    Sourced from main.R
#===============================================================================

# 1. Load libraries
#-------------------------------------------------------------------------------
library(tidyverse)
library(biomaRt)
library(ggbeeswarm)
library(fgsea)
library(msigdbr)
library(ggpubr)
library(ggnewscale)
library(conflicted)
library(ggrepel)

# Load custom functions and parameters
source(file.path(code_dir,"functions","plot_DEG_GSEA_fun.R"))

conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(base::as.factor)
#-------------------------------------------------------------------------------

# 2. Load data
#-------------------------------------------------------------------------------
output_dgea_dir <- file.path(results_dir, "dgea")
if (!dir.exists(output_dgea_dir)){dir.create(output_dgea_dir, recursive = T)}
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

# Load DEGs
dgea_met <- read.csv(paste0(output_dgea_dir,"/~Biopsy_excision+Metastasis/DGEA_results_DESeq2.csv"))

#-------------------------------------------------------------------------------

# 3. Generate figures
#-------------------------------------------------------------------------------
set.seed(123)

# Hallmark gene sets (Ensembl IDs)
all_gene_sets = msigdbr(species = "Homo sapiens",category = "H")
H_pathways = split(x = all_gene_sets$ensembl_gene, f = all_gene_sets$gs_name)

# Cases versus controls
runGSEA(dgea_met,
        pathways_to_test = H_pathways,
        design = "~Sample_type+Metastasis",
        figure_name = "figure_4d")

#-------------------------------------------------------------------------------


