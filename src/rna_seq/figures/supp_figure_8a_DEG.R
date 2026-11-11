#===============================================================================
# Aim:      Generate Supplementary Figure 8A: volcano plot of the DEGs of
#           immunosuppressed versus immunocompetent patients.
# Input:    - DGEA results ~Biopsy_excision+Immunosupp_bin (results_dir/dgea)
# Output:   - Volcano plot (pdf) and plotted data (csv) in results_dir/publication
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
dgea_immuno <- read.csv(paste0(output_dgea_dir,"/~Biopsy_excision+Immunosupp_bin/DGEA_results_DESeq2.csv"))

#-------------------------------------------------------------------------------

# 3. Generate figures
#-------------------------------------------------------------------------------
set.seed(123)

# Plot DEGs: immunosuppressed versus immunocompetent
plotsSignatures(dgea_immuno,
                design = "~Sample_type+Immunosupp_bin",
                figure_name = "figure_8a")


#-------------------------------------------------------------------------------
