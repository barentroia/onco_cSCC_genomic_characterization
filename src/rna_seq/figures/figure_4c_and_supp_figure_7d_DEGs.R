#-------------------------------------------------------------------------------
# Aim: Generate figure 4c and supp figure 5d DEG (w/wo DvP genes)
# Input: DEG metastasis
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
dgea_met <- read.csv(paste0(output_dgea_dir,"/~Biopsy_excision+Metastasis/DGEA_results_DESeq2.csv"))

# Bailey genes
Bailey_gene_sig <- read.csv(p6_bailey,sep=";")
Bailey_gene_sig_prog <- unique(Bailey_gene_sig$Progenitor)
Bailey_gene_sig_diff <- unique(Bailey_gene_sig$Differentiated)

#-------------------------------------------------------------------------------


# Generate figures
#-------------------------------------------------------------------------------
set.seed(123)

# Plot DEGs
# Cases vs controls
plotsSignatures(dgea_met,
                design = "~Sample_type+Metastasis",
                figure_name = "figure_4c")

plotsSignatures(dgea_met,
                signature=list("DvP_differentiated"=Bailey_gene_sig_diff,"DvP_progenitor"=Bailey_gene_sig_prog),
                design="~Sample_type+Metastasis",
                figure_name = "figure_7d")

#-------------------------------------------------------------------------------


