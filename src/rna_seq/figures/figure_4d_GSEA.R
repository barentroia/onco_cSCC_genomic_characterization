#-------------------------------------------------------------------------------
# Aim: Generate figure 4d
# Input: DEG metastasis
# Output: GSEA barplot
#-------------------------------------------------------------------------------


# Load libraries
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


# Load data
#-------------------------------------------------------------------------------
output_dgea_dir <- file.path(results_dir, "dgea")
if (!dir.exists(output_dgea_dir)){dir.create(output_dgea_dir, recursive = T)}
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

# Load DEGs
dgea_met <- read.csv(paste0(output_dgea_dir,"/~Biopsy_excision+Metastasis/DGEA_results_DESeq2.csv"))

#-------------------------------------------------------------------------------


# Generate figures
#-------------------------------------------------------------------------------
set.seed(123)

##### Hallmarks GSEA
all_gene_sets = msigdbr(species = "Homo sapiens",category = "H")
H_pathways = split(x = all_gene_sets$ensembl_gene, f = all_gene_sets$gs_name)

# Immunosuppressed versus controls
runGSEA(dgea_met,
        pathways_to_test = H_pathways,
        design = "~Sample_type+Metastasis",
        figure_name = "figure_4d")

#-------------------------------------------------------------------------------


