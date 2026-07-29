#-------------------------------------------------------------------------------
# Aim: Generate figure 8a
# Input: DEG immunocompromised bin
# Output: Volcano plots
#-------------------------------------------------------------------------------


# Load libraries
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


# Load data
#-------------------------------------------------------------------------------
output_dgea_dir <- file.path(results_dir, "dgea")
if (!dir.exists(output_dgea_dir)){dir.create(output_dgea_dir, recursive = T)}
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

# Load DEGs
dgea_immuno <- read.csv(paste0(output_dgea_dir,"/~Biopsy_excision+Immunosupp_bin/DGEA_results_DESeq2.csv"))

#-------------------------------------------------------------------------------


# Generate figures
#-------------------------------------------------------------------------------
set.seed(123)

# Plot DEGs
# Cases vs controls
plotsSignatures(dgea_immuno,
                design = "~Sample_type+Immunosupp_bin",
                figure_name = "figure_8a")


#-------------------------------------------------------------------------------
