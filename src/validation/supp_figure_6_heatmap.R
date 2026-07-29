#-------------------------------------------------------------------------------
# Aim: Generate Figure S6 heatmaps
# Input: cSCC sequencing data validation + Nassir data set
# Output: heatmaps
#-------------------------------------------------------------------------------


# Load libraries
#-------------------------------------------------------------------------------
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)
library(ggnewscale)
library(conflicted)
library(rtracklayer)

# Load custom functions and parameters
source(file.path("code","functions","dgea_deseq2_badqc_filtering_strict.R"))
source(file.path("code","functions","colors.R"))
source(file.path("code","functions","plotFuncs.R"))

conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
#-------------------------------------------------------------------------------


# General setup
#-------------------------------------------------------------------------------

# output folder
output_dir <- file.path("output", "validation")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

#-------------------------------------------------------------------------------


# Heatmap plot
#-------------------------------------------------------------------------------
# Load validation data
validation_data <- readRDS(paste0(output_data_dir,"/Validation_data_set.rds"))
#validation_data_nassir

#-------------------------------------------------------------------------------


