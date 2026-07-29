#---------------------------------------------------
# Aim: Experiment setup for different experiments
# This script aims to setup all the necessary data for the various modeling experiments.
# Author: B. Rentroia Pacheco
# Input: Clinical dataset
#        TPM gene expression matrix
#        Raw count gene expression matrix
#        Tumor purity data
#        QC data
#        Bad QC sample list
#        Sample weights
#        Gene annotation list with biotype LncRNA or protein coding
# Output: 
#        Processed clinical data
#        Prefiltered and processed TPM expression matrix
#        Prefiltered raw count expression matrix
#        QC exclusion list

#--------------------------------------------------

#0a. Libraries
#---------------------------------------------------
suppressPackageStartupMessages({
library(tidyverse)
library(glmnet)
library(biomaRt)
#library(VennDiagram)
library(gridExtra)
#library(randomForest)
library(ranger)
library(caret)
library(randomForestSRC)
#library(survivalsvm)
library(DESeq2)
library(class)
#library(limma)
#library(MatchIt)
library(parsnip)
library(dials)
library(workflows)
#library(bonsai)
library(tune)
#library(recipes)
library(yardstick)
#library(xgboost)
library(Hmisc)
library(survival)
library(WeightedROC)
library(intsurv)
library(spatstat)
library(readxl)
library(readr)
library(dplyr)
library(patchwork)
library(svglite)
library(parallel)
library(purrr)
library(conflicted)
library(stringr)
library(ggplot2)
library(ggtext)
library(rtracklayer)
})
suppressMessages({
conflicts_prefer(base::as.factor)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("setdiff", "base")
conflicted::conflict_prefer("intersect", "base")
conflicted::conflict_prefer("setequal", "base")
conflicted::conflict_prefer("cbind", "base")
conflicted::conflict_prefer("rbind", "base")
conflicted::conflict_prefer("rename", "dplyr")

})

#---------------------------------------------------

# 0b. Functions
#-------------------------------------------------------------------------------
# save_model_bootstrap: This function save the workflow's results
# @param b_model: workflow's output
# @param model_id: lmodel's name
# @param filedir: directory where to save results
save_model_bootstrap <- function(b_model, model_id, filedir){
  # Save R object
  saveRDS(b_model, file.path(filedir, paste0(model_id, "_robj.rds")))
}
#-------------------------------------------------------------------------------

# 0c. Source functions
#-------------------------------------------------------------------------------
if (source_type == "early"){
  source(file.path(dir_scripts_functions, "auxilliary_modeleval.R"))
}else if (source_type == "late"){
  source(file.path(dir_scripts_functions, "auxilliary_combination.R"))
}else if (source_type == "late_wes"){
  source(file.path(dir_scripts_functions, "auxilliary_combination_integration.R"))
}
source(file.path(dir_scripts_functions,  "build_models.r"))
source(file.path(dir_scripts_functions, "apply_models.r"))

#-------------------------------------------------------------------------------

# 1. Load the data
#-------------------------------------------------------------------------------
#########################################
### TO DO CHANGE INPUT DATA LOCATION
# Gene expression data
gene_exp_data <-  suppressMessages(read.csv(p1_tpms))
# Gene expression counts data
x_data_counts <-  suppressMessages(read.csv(p2_counts))
# Tumor purity data
#tp_data <- read.csv(p11_tp_data)
# QCs data
#qcs_data <- read.csv(p12_qcs_data)
# Samples with bad QCs
samples_bad_qcs <-  suppressMessages(read.csv(p7_badQC))
# Clinical data
clinical_data <-  suppressMessages(read.csv(file.path(results_dir, "intermediate", "clinical","merged_clinical_data_imputed_d_squame_discovery.csv")))
clinical_data <- clinical_data[order(clinical_data$SkylineDx.ID),]
# Weights data
# Use 195 for unmatched samples and 183 for matched samples
if (samples_matched == "unmatched"){
  weights_ncc_fn <- p11_weights_unmatched
}else{
  weights_ncc_fn <- p10_weights
}
weights_ncc <-  suppressMessages(read.csv(weights_ncc_fn))

# Genes ids with biotype "protein_coding" or "lncRNA": Retrieved from biomart on 2024-10-14 15:40:04 CEST
genes_to_keep <- readRDS(p12_genes_to_keep)
# Get DvP Bailey signature
bailey_diff_prog <- suppressMessages(read_csv(p6_bailey))
# WES mutation dataset
mutation_data <-  suppressMessages(read_excel(p13_mutation_data))
# Potential pathways to use in WES data
all_pathways <-  suppressMessages(read_excel(p14_all_pathways))
# 138 WES samples with good quality
input138 <-  suppressMessages(read_csv(p15_input138))
# Load GTF
gtf <- import(p9_gtf)
#-------------------------------------------------------------------------------

# 1. Preprocess DvP signature
#-------------------------------------------------------------------------------
# Create Ensembl gene ID -> external gene name mapping
gene_mappings <- as.data.frame(mcols(gtf)) %>%
  filter(type == "gene") %>%
  select(
    ensembl_gene_id = gene_id,
    external_gene_name = gene_name
  ) %>%
  distinct()

# Split Differentiated and Progenitor genes,
# then convert gene names to Ensembl IDs
bailey_gene_ids <- bailey_diff_prog %>%
  separate(
    `Differentiated;Progenitor`,
    into = c("Differentiated", "Progenitor"),
    sep = ";"
  ) %>%
  pivot_longer(
    cols = c(Differentiated, Progenitor),
    names_to = "cell_type",
    values_to = "external_gene_name"
  ) %>%
  filter(
    !is.na(external_gene_name),
    external_gene_name != ""
  ) %>%
  left_join(
    gene_mappings,
    by = "external_gene_name"
  ) %>%
  select(
    ensembl_gene_id,
    external_gene_name
  ) %>%
  distinct()
bailey_gene_ids$ensembl_gene_id <- strp(bailey_gene_ids$ensembl_gene_id )
#-------------------------------------------------------------------------------

# 2. Preprocess clinical
#-------------------------------------------------------------------------------
# Preprocess clinical data
clinical_data <- clinical_data %>%
  ## Combine tumor purity with clinical data
 # full_join(tp_data %>%
#              dplyr::select(SkylineDx.ID, ESTIMATE) %>%
#              dplyr::rename(ESTIMATE_tumorpurity = ESTIMATE),
#            by = "SkylineDx.ID") %>%
  ## Combine weights with clinical data
  full_join(weights_ncc %>%
              dplyr::select(SkylineDx.ID, Weight_rescaled),
            by = "SkylineDx.ID") %>%
  ## Add QCs info
 # left_join(qcs_data,
#            by = c("SkylineDx.ID" = "Sample")) %>%
  mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0),
         # Truncation at 5 years
         Metastasis_num_trunc5 = ifelse(FU_metastasis_years > 5, 0 , Metastasis_num),
         FU_metastasis_years_trunc5 = ifelse(FU_metastasis_years > 5, 5, FU_metastasis_years),
         # Pairs are factors
         Set_id = as.factor(Set_id),
         # Differentiation as factor and rename groups
         Differentiation = ifelse(Differentiation == "Good/moderate", "Good_moderate",
                                  ifelse(Differentiation == "Poor/undifferentiated", "Poor_undifferentiated", Differentiation)),
         Differentiation = factor(Differentiation, levels = c("Good_moderate", "Poor_undifferentiated")),
         # Tissue involvement as factor and rename groups
         Tissue_involvement = ifelse(Tissue_involvement == "Subcutaneous fat", "Subcutaneous_fat",
                                     ifelse(Tissue_involvement == "Beyond subcutaneous fat", "Beyond_subcutaneous_fat", Tissue_involvement)),
         Tissue_involvement = factor(Tissue_involvement, levels = c("Dermis", "Subcutaneous_fat", "Beyond_subcutaneous_fat")),
         # Tumor location as factor and rename groups
         Tumor_location = factor(ifelse(Tumor_location == "Trunk/Extremities", "Trunk_extremities",
                                        ifelse(Tumor_location == "Scalp/neck", "Scalp_neck", Tumor_location)),
                                levels = c( "Trunk_extremities", "Scalp_neck", "Face")),
         # Sex as factor
         Sex = factor(Sex, levels = c("Female", "Male")),
         # PNI_or_LVI as factor
         PNI_or_LVI = factor(PNI_or_LVI, levels = c("No", "Yes")),
         # Number_of_cSCC_before_culprit as integer
         Number_of_cSCC_before_culprit = as.integer(Number_of_cSCC_before_culprit),
         # AJCC as factor
         AJCC_8 = factor(AJCC_8, levels = c("T1", "T2", "T3", "T4")),
         # BWH as factor
         BWH = factor(BWH, levels = c("T1", "T2a", "T2b", "T3")),
         # Binarize AJCC staging
         AJCC_bin = factor(ifelse(AJCC_8 != "T1", "T2_T3_T4", "T1"), levels = c("T1", "T2_T3_T4")),
         # Binarize BWH staging
         BWH_bin = factor(ifelse(BWH != "T1", "T2a_T2b_T3", "T1"), levels = c("T1", "T2a_T2b_T3")),
         # Binarize AJCC staging (2)
         AJCC_bin2 = factor(ifelse(AJCC_8 %in% c("T1", "T2"), "T1_T2", "T3_T4"), levels = c("T1_T2", "T3_T4")),
         # Binarize BWH staging (2)
         BWH_bin2 = factor(ifelse(BWH %in% c("T1", "T2a"), "T1_T2a", "T2b_T3"), levels = c("T1_T2a", "T2b_T3")),
         # Binarize sample type
         Biopsy_excision_2f = factor(ifelse(Biopsy_excision == "Biopsy", "Biopsy", "Excision"),
                                     levels = c("Biopsy", "Excision")),
         # Binarize number of prior cSCC
         Number_of_cSCC_before_culprit_bin = factor(ifelse(Number_of_cSCC_before_culprit == 0, "0", "0+"),
                                                    levels = c("0","0+")),
         # Combine Number_of_cSCC_before_culprit and stages variables
         AJCC_and_Number_of_prior_cSCC = paste0(AJCC_8, "_", Number_of_cSCC_before_culprit_bin),
         BWH_and_Number_of_prior_cSCC = paste0(BWH, "_", Number_of_cSCC_before_culprit_bin),
         # Create immunosuppression variable
         Immunosuppressed = factor(ifelse(HM_at_cSCC == "Yes" | OTR_at_cSCC == "Yes", "Yes","No"),
                                   levels = c("No","Yes")),
         # Divide % trimmed into quartiles
         #Percent_trimmed_quartiles = ifelse(percent_trimmed < quantile(percent_trimmed, probs = 0.25, type = 1), "1st",
        #                                    ifelse(percent_trimmed < quantile(percent_trimmed, probs = 0.5, type = 1), "2nd",
        #                                           ifelse(percent_trimmed < quantile(percent_trimmed, probs = 0.75, type = 1), "3rd",
        #                                                  "4th"))),
         # Divide sequencing depth into quartiles
         #Seq_depth_quartiles = ifelse(Tot.Seqs < quantile(Tot.Seqs, probs = 0.25, type = 1), "1st",
        #                              ifelse(Tot.Seqs < quantile(Tot.Seqs, probs = 0.5, type = 1), "2nd",
        #                                     ifelse(Tot.Seqs < quantile(Tot.Seqs, probs = 0.75, type = 1), "3rd",
        #                                            "4th"))),
         # Divide ESTIMATE tumor purity into quartiles
         #ESTIMATE_tumorpurity_quartiles = ifelse(ESTIMATE_tumorpurity < quantile(ESTIMATE_tumorpurity, probs = 0.25, type = 1), "1st",
        #                                         ifelse(ESTIMATE_tumorpurity < quantile(ESTIMATE_tumorpurity, probs = 0.5, type = 1), "2nd",
        #                                                ifelse(ESTIMATE_tumorpurity < quantile(ESTIMATE_tumorpurity, probs = 0.75, type = 1), "3rd",
        #                                                       "4th"))),
         # Divide samples into 2 groups, based on the year in which material was obtained
         Year_of_obtained_material_bin = ifelse(as.numeric(Year_of_obtained_material) <= 2009, "<=2009", ">2009"),
         # Binarize based on number of NCCN very high risk factors
         hr_t_diam = ifelse(Tumor_diameter > 4, 1, 0),
         hr_differentiation = ifelse(Differentiation == "Poor_undifferentiated", 1, 0),
         hr_pni = ifelse(PNI_bin == "Yes", 1, 0),
         hr_t_involvement = ifelse(Tissue_involvement == "Beyond_subcutaneous_fat" | Depth_of_Invasion > 6, 1, 0),
         hr_lvi = ifelse(Lymphovascular_invasion_bin == "Yes", 1, 0),
         N_very_high_risk_features = hr_t_diam + hr_differentiation + hr_pni + hr_t_involvement + hr_lvi,
         N_very_high_risk_features_bin = ifelse(N_very_high_risk_features == 0, "0", ">=1")) %>%
  dplyr::select(-hr_t_diam, -hr_differentiation, -hr_pni, -hr_t_involvement, -hr_lvi, N_very_high_risk_features)

# Preprocess gene expression data
feature_names <- gene_exp_data[, 1]
gene_exp_data <- t(gene_exp_data[, -1])
colnames(gene_exp_data) <- feature_names
gene_exp_data <- as.data.frame(gene_exp_data)
## Align gene expression and clinical datasets
gene_exp_data <- gene_exp_data[as.character(clinical_data$SkylineDx.ID),]
stopifnot(all(rownames(gene_exp_data) == clinical_data$SkylineDx.ID))
rownames(clinical_data) <- clinical_data$SkylineDx.ID
## If remove_badQCsamples == T, identify samples with bad QC and remove them

# Remove whole sets with bad samples if matched or only bad samples if unmatched
if (samples_matched == "unmatched"){
  sample_ids_bad_qcs <- clinical_data %>%
    dplyr::filter(SkylineDx.ID %in% samples_bad_qcs$Skyline_ID) %>%
    pull(SkylineDx.ID)
}else{
  sample_ids_bad_qcs <- clinical_data %>%
    dplyr::filter(Set_id %in% samples_bad_qcs$Set_id) %>%
    pull(SkylineDx.ID)
}
#-------------------------------------------------------------------------------

# 3. Preprocess gene expression data
#-------------------------------------------------------------------------------

  gene_exp_data_for_filtering <- gene_exp_data[setdiff(as.character(clinical_data$SkylineDx.ID), sample_ids_bad_qcs),]  
## Pre-filtering
feature_exclusion_df <- data.frame("Gene_name" = colnames(gene_exp_data_for_filtering),
                                   "Zero_variance" = 0,
                                   "Low_variance_0_01" = 0,
                                   "Filter_logtpm_p1" = 0,
                                   "Protein_lnc" = 0,
                                   "Low_mean_logtpm" = 0)
## Features with constant variance
variance_GEfeats <- apply(gene_exp_data_for_filtering, 2, var)
feature_exclusion_df$Zero_variance[which(variance_GEfeats == 0)] <- 1
feature_exclusion_df$Low_variance_0_01[which(variance_GEfeats < 0.01)] <- 1

# Filter genes where average logtpm should be 1 in most of the samples, instead of 0.5 in half of the samples
filter_logtpm <- apply(gene_exp_data_for_filtering, 2, function(x) sum(log2(x + 1) > 1)) > (nrow(gene_exp_data_for_filtering) * 0.5)
feature_exclusion_df$Filter_logtpm_p1[!filter_logtpm] <- 1

## Mean log TPM filter: if average of gene in samples has log tpm of 2.65, it will be removed
low_mean_logtpm <- apply(gene_exp_data_for_filtering, 2, function(x) mean(log2(x + 1)) < 2.65)
feature_exclusion_df$Low_mean_logtpm[low_mean_logtpm] <-1

## Identify protein coding genes
feature_exclusion_df$Protein_lnc[!gsub("\\..*","",colnames(gene_exp_data_for_filtering)) %in% genes_to_keep$ensembl_gene_id] <- 1
## Combine filters
feature_exclusion_df$filter_logp_and_protein_lnc <- ifelse(
  feature_exclusion_df$Filter_logtpm_p1 == 1 |
    feature_exclusion_df$Protein_lnc == 1 |
    feature_exclusion_df$Low_mean_logtpm == 1 ,
  1,
  0
)

# Filter count data from 390 to 366 genes
x_data_countsf <- x_data_counts %>%
  column_to_rownames(var = "ID") %>%
  dplyr::select(all_of(clinical_data$SkylineDx.ID))
ds_x_counts_train <- x_data_countsf[, !(colnames(x_data_countsf) %in% sample_ids_bad_qcs)]

## Apply filtering
x_data <- data.matrix(gene_exp_data)[, intersect(which(feature_exclusion_df[, "filter_logp_and_protein_lnc"] == 0),which(feature_exclusion_df$Zero_variance == 0))]
## log2-transformation
x_data <- log2(x_data + 1)


# Preprocess count data
x_data_counts <- x_data_counts %>%
  as.data.frame() %>%
  dplyr::filter(ID %in% colnames(x_data)) %>%
  column_to_rownames(var = "ID") %>%
  dplyr::select(all_of(clinical_data$SkylineDx.ID))
## Align datasets
stopifnot(all(colnames(x_data_counts) == clinical_data$SkylineDx.ID))
stopifnot(all(colnames(x_data_counts) == rownames(x_data)))

#-------------------------------------------------------------------------------

# 4. Preprocess mutation data
#-------------------------------------------------------------------------------
if (experiment_setup == "integration"){
  # Check which pathways can be used
  all_pathways2 <- all_pathways %>%
    group_by(Pathway_publication) %>%
    mutate(n_pub = n()) %>%
    ungroup() %>%
    mutate(result = if_else(n_pub > 1, Pathway_publication, Gene)) %>%
    pull(result)  %>%
    unique()
  pathway_cols <- intersect(all_pathways2, colnames(mutation_data))
  
  # Process CN (chromosome) columns
  # Find CN columns: those starting with 'chr' and ending with 'p' or 'q'
  cn_cols <- grep("^chr.*[pq]$", colnames(mutation_data), value = TRUE)
  
  for (col in cn_cols) {
    # b. Binary for amplifications or deletions only (1 or -1 -> 1)
    mutation_data[[paste0(col, "_amp_bin")]] <- ifelse(mutation_data[[col]] == 1, 1, 0)
    
    mutation_data[[paste0(col, "_del_bin")]] <- ifelse(mutation_data[[col]] == -1, 1, 0)
  }
  # Remove bad wes samples
  mutation_data <- mutation_data[mutation_data$Sample %in% input138$...1,]
  
  # Overlap with GEP samples
  rownames(x_data_counts) <- strp(rownames(x_data_counts))
  colnames(x_data) <- strp(colnames(x_data))
  ids <- mutation_data$SkylineDx.ID
  keep <- colnames(x_data_counts) %in% ids
  x_data_counts <- x_data_counts[, keep]
  x_data <- x_data[keep, ]
  
  # Now align clinical_data to the same samples and order
  common_samples <- intersect(clinical_data$SkylineDx.ID, colnames(x_data_counts))
  clinical_data <- clinical_data[clinical_data$SkylineDx.ID %in% common_samples, ]
  clinical_data <- clinical_data[match(colnames(x_data_counts), clinical_data$SkylineDx.ID), ]
  stopifnot(all(clinical_data$SkylineDx.ID == colnames(x_data_counts)))
  # Bad samples will be removed later
  cat("Number of samples kept:", nrow(clinical_data), "\n")
  
  # Pathways + chr*_amp_bin and chr*_del_bin
  mutation_data_amp_del <- mutation_data %>%
    select(all_of(pathway_cols),SkylineDx.ID,Tumor.cellularity.avg.pct,
           matches("^chr.*_(amp|del)_bin$"))
  
  # Make sure that mutation data has the same samples and order as clinical data
  mutation_data_amp_del <- clinical_data %>%
    select(SkylineDx.ID) %>%
    left_join(mutation_data_amp_del, by = "SkylineDx.ID")
  rownames(mutation_data_amp_del) <- mutation_data_amp_del$SkylineDx.ID
  mutation_data_amp_del$SkylineDx.ID  <- NULL
  all(clinical_data$SkylineDx.ID == rownames(mutation_data_amp_del))
  
}


