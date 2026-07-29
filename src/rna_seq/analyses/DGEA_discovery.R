#-------------------------------------------------------------------------------
# Aim: Differentially expressed genes in discovery
# Input: discovery sequencing data object
# Output: DEG
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
source(file.path(code_dir,"functions","dgea_fun.R"))

conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(base::as.factor)
#-------------------------------------------------------------------------------



# General setup
#-------------------------------------------------------------------------------
# output folder
output_dgea_dir <- file.path(results_dir, "dgea")
if (!dir.exists(output_dgea_dir)){dir.create(output_dgea_dir, recursive = T)}
#-------------------------------------------------------------------------------


# DGEA 
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

