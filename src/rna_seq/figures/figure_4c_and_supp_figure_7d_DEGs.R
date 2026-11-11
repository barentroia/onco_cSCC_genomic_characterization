#===============================================================================
# Aim:      Generate Figure 4C and Supplementary Figure 7D: volcano plots of
#           the DEGs of cases versus controls, without and with the DvP
#           (Bailey) genes highlighted.
# Input:    - DGEA results ~Biopsy_excision+Metastasis (results_dir/dgea)
#           - DvP signature (p6_bailey)
# Output:   - Volcano plots (pdf) and plotted data (csv) in results_dir/publication
# Requires: code_dir, results_dir (main.R); DGEA_discovery.R run first
# Usage:    Sourced from main.R
#===============================================================================

# 1. Load libraries
#-------------------------------------------------------------------------------
library(tidyverse)
library(ggbeeswarm)
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

# DvP (Bailey) genes: progenitor and differentiated
Bailey_gene_sig <- read.csv(p6_bailey,sep=";")
Bailey_gene_sig_prog <- unique(Bailey_gene_sig$Progenitor)
Bailey_gene_sig_diff <- unique(Bailey_gene_sig$Differentiated)

#-------------------------------------------------------------------------------

# 3. Generate figures
#-------------------------------------------------------------------------------
set.seed(123)

# Plot DEGs of cases versus controls: Figure 4C (all DEGs) and
# Supplementary Figure 7D (DvP genes highlighted)
plotsSignatures(dgea_met,
                design = "~Sample_type+Metastasis",
                figure_name = "figure_4c")

plotsSignatures(dgea_met,
                signature=list("DvP_differentiated"=Bailey_gene_sig_diff,"DvP_progenitor"=Bailey_gene_sig_prog),
                design="~Sample_type+Metastasis",
                figure_name = "figure_7d")

#-------------------------------------------------------------------------------


