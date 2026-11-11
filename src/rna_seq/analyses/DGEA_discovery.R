#===============================================================================
# Aim:      Differential gene expression analysis (DESeq2) in the D-SQUAME
#           discovery dataset (N = 378): metastasis (case vs control) and
#           immunosuppression (yes vs no), both adjusted for sample type.
# Input:    - Discovery_data_set_N378.rds (output_data_dir)
# Output:   - DESeq2 results (csv, RData) per design in results_dir/dgea/<design>
# Requires: code_dir, results_dir, output_data_dir (main.R)
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
source(file.path(code_dir,"functions","dgea_fun.R"))

conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(base::as.factor)
conflicts_prefer(base::setdiff)

#-------------------------------------------------------------------------------

# 2. General setup
#-------------------------------------------------------------------------------
# output folder
output_dgea_dir <- file.path(results_dir, "dgea")
if (!dir.exists(output_dgea_dir)){dir.create(output_dgea_dir, recursive = T)}
#-------------------------------------------------------------------------------

# 3. DGEA
#-------------------------------------------------------------------------------
discovery_data <- readRDS(paste0(output_data_dir,"/Discovery_data_set_N378.rds"))

p6_design <- "~Biopsy_excision+Metastasis"
dgea_met <- dgea_deseq2(rnaseq_obj = discovery_data, 
                              design = p6_design, 
                              output_folder = paste0(output_dgea_dir,"/",p6_design),
                              contrast_design=c("Metastasis","Case","Control"))


p6_design <- "~Biopsy_excision+Immunosupp_bin"
dgea_met <- dgea_deseq2(rnaseq_obj = discovery_data, 
                              design = p6_design, 
                              output_folder = paste0(output_dgea_dir,"/",p6_design),
                              contrast_design=c("Immunosupp_bin","Yes","No"))

#-------------------------------------------------------------------------------

