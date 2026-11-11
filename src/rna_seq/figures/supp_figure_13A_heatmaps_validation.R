#===============================================================================
# Aim:      Generate Supplementary Figure 13A (heatmaps): scaled log2(TPM + 1)
#           expression of the 23 SCCore-GEP genes in the Nassir et al. and
#           D-SQUAME validation datasets, samples ordered by SCCore-GEP
#           prediction and genes by coefficient.
# Input:    - Final SCCore-GEP model (coxnet_Late_feature_integration_robj.rds)
#           - SCCore-GEP predictions in both validation datasets
#             (results_dir/intermediate/validation/<dataset>)
#           - D_SQUAME_validation_data_set.rds, Nassir_validation_data_set.rds
# Output:   - Heatmaps (pdf) and source data (xlsx) in results_dir/publication
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
library(rtracklayer)
library(ComplexHeatmap)

# Load custom functions and parameters
source(file.path(code_dir,"functions","colors.R"))

conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
#-------------------------------------------------------------------------------

# 2. General setup
#-------------------------------------------------------------------------------
# output folder
output_data_dir <- file.path(results_dir,  "data_rnaseq")
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# 3. Final SCCore-GEP model (23 genes)
#-------------------------------------------------------------------------------
sccoregep_fn <- file.path(results_dir, "intermediate", "sccore_gep","Discovery_models",
                          "coxnet_Late_feature_integration_robj.rds")
GEP_final_model <- readRDS(sccoregep_fn)
final_model_genes <- names(GEP_final_model$Final_model2$coefficients)
stopifnot(length(final_model_genes)==23)
#-------------------------------------------------------------------------------

# 4. SCCore-GEP predictions in D-SQUAME validation
#-------------------------------------------------------------------------------
p_pred_val_emc <- file.path(results_dir2, "intermediate","validation", "d_squame", "SCCoreGEP_predictions.csv")
predictions_validation_emc <- read.csv(p_pred_val_emc)
p_clinical_imputed_val <- file.path(results_dir, "intermediate", "clinical",
                                    "merged_clinical_data_imputed_d_squame_validation.csv")
#-------------------------------------------------------------------------------

# 5. SCCore-GEP predictions in Nassir et al.
#-------------------------------------------------------------------------------
p_pred_val_nassir <- file.path(results_dir, "intermediate","validation", "nassir", "SCCoreGEP_predictions.csv")
predictions_validation_nassir <- read.csv(p_pred_val_nassir)
#-------------------------------------------------------------------------------

# 6. Heatmap plots
#-------------------------------------------------------------------------------
# Load validation data
validation_data <- readRDS(paste0(output_data_dir,"/D_SQUAME_validation_data_set.rds"))
nassir_validation_data <- readRDS(paste0(output_data_dir,"/Nassir_validation_data_set.rds"))
out <- validation_data@out

# Genes and coefficients from the final model
coeff_final_model <- GEP_final_model$Final_model2$coefficients
final_model_ensg <- out[match(names(coeff_final_model),do.call(rbind,strsplit(out$gene_id,"[.]"))[,1]),]$gene_id
names(coeff_final_model) <- out[match(names(coeff_final_model),do.call(rbind,strsplit(out$gene_id,"[.]"))[,1]),]$gene_name


# 6.1 Heatmap of the Nassir et al. dataset
tpm_data_nassir_model <- nassir_validation_data@sequencing_data$tpms %>%
  filter(rownames(nassir_validation_data@sequencing_data$tpms) %in% final_model_ensg)
stopifnot(nrow(tpm_data_nassir_model) > 0)            # model genes found in the expression data
tpm_data_nassir_model <- log2(tpm_data_nassir_model + 1)
rownames(tpm_data_nassir_model) <- out[match(rownames(tpm_data_nassir_model), out$gene_id),]$gene_name

# Clinical data: outcome as Case/Control
clinical_nassir <- nassir_validation_data@clinical
clinical_nassir$Metastasis <- ifelse(clinical_nassir$MetNoMet == "Met", "Case", "Control")

# Scale expression per gene
all_models_scaled <- t(scale(t(tpm_data_nassir_model)))
stopifnot(ncol(all_models_scaled) > 0)                # samples present

# Order samples by SCCore-GEP prediction (highest risk first)
stopifnot(nrow(predictions_validation_nassir) > 0)    # predictions file read
predictions_nassir_filt <- predictions_validation_nassir[predictions_validation_nassir$Sample_id %in% colnames(all_models_scaled),]
stopifnot(nrow(predictions_nassir_filt) > 0)          # predictions matched to expression samples
ordered_predictions <- predictions_nassir_filt[order(predictions_nassir_filt$sccore_gep, decreasing = T),]
all_models_scaled_ordered <- all_models_scaled[, ordered_predictions$Sample_id]

# Clinical data in the same order as the heatmap columns
clinical_data_nassir_ordered <- clinical_nassir[match(ordered_predictions$Sample_id, clinical_nassir$Run),]
stopifnot(nrow(clinical_data_nassir_ordered) == nrow(ordered_predictions))
stopifnot(identical(clinical_data_nassir_ordered$Run, ordered_predictions$Sample_id))

# BWH as used in the Nassir et al. analysis (SCCoreGEP_BWH_performances.R):
# stage taken from the sample name; samples without a stage are T2b
clinical_data_nassir_ordered$BWH <- str_extract(clinical_data_nassir_ordered$SampleID, "T[0-9]+[a-z]?") %>%
  replace_na("T2b")
stopifnot(!any(is.na(clinical_data_nassir_ordered$BWH)))

clinical_data_nassir_ordered$rank <- c(1:nrow(clinical_data_nassir_ordered))
col_ha = HeatmapAnnotation(Metastasis = clinical_data_nassir_ordered$Metastasis, 
                           BWH = clinical_data_nassir_ordered$BWH,
                           rank = clinical_data_nassir_ordered$rank,
                           col=c(cols_clinical_data,rank = circlize::colorRamp2(c(0,ncol(all_models_scaled_ordered)),colors = c("white","#0e7983"))))

# order rows by coefficient
coeff_final_model_ordered <- coeff_final_model[order(coeff_final_model)]
all_models_scaled_ordered <- all_models_scaled_ordered[match(names(coeff_final_model_ordered),rownames(all_models_scaled_ordered)),]
colnames(all_models_scaled_ordered) <- NULL


ht <-Heatmap(all_models_scaled_ordered, top_annotation = col_ha, name="scaled", 
             cluster_columns = FALSE,cluster_rows = FALSE)
pdf_file <- file.path(paste0(output_dir,"/Supp_figure_13A_heatmap_logTPMs_scaled_Nassir_23genes_model_Romy_ungrouped_prediction_genes_rank.pdf"))
pdf(file = pdf_file, width = 7, height = 5, colormodel = "cmyk")
draw(ht)
dev.off()

pdf_file
# Source data
df_nassir <- as.data.frame(t(all_models_scaled_ordered))
df_nassir$rank <- clinical_data_nassir_ordered$rank
df_nassir$BWH <- clinical_data_nassir_ordered$BWH
df_nassir$Metastasis <- clinical_data_nassir_ordered$Metastasis
df_nassir$Sample_id <- clinical_data_nassir_ordered$SampleID

writexl::write_xlsx(df_nassir,paste0(output_dir,"/Supp_figure_13A_heatmap_nassir_source_data.xlsx"))

# 6.2 Heatmap of the D-SQUAME validation dataset
tpm_data_emc_model <- validation_data@sequencing_data$tpms %>% filter(rownames(validation_data@sequencing_data$tpms) %in% final_model_ensg)
tpm_data_emc_model <- log2(tpm_data_emc_model + 1) #log scale
rownames(tpm_data_emc_model) <- out[match(rownames(tpm_data_emc_model),out$gene_id),]$gene_name

# Scale data per gene
all_models_scaled <- t(scale(t(tpm_data_emc_model)))

# Fix sample names
predictions_validation_emc$tpm_sample_id <- validation_data@clinical[match(predictions_validation_emc$Skyline_ID,validation_data@clinical$Skyline_ID),]$tpm_sample_id

# Order samples by prediction (highest risk first)
#predictions_emc_filt <- predictions_validation_emc[predictions_validation_emc$tpm_sample_id %in% colnames(all_models_scaled),]
predictions_emc_filt <- predictions_validation_emc[predictions_validation_emc$Skyline_ID  %in% colnames(all_models_scaled),]
ordered_predictions <- predictions_emc_filt[order(predictions_emc_filt$sccore_gep, decreasing = T),]
#all_models_scaled_ordered <- all_models_scaled[,ordered_predictions$tpm_sample_id]
all_models_scaled_ordered <- all_models_scaled[,ordered_predictions$Skyline_ID]
clinical_data_emc_ordered <- validation_data@clinical[match(ordered_predictions$Skyline_ID,validation_data@clinical$Skyline_ID),]
# BWH from the multiple imputation (as in the published figure): the raw
# clinical data has missing BWH for patients with missing staging variables
clinical_imputed_val <- read.csv(p_clinical_imputed_val)
clinical_data_emc_ordered$BWH <- as.character(
  clinical_imputed_val$BWH_staging[match(clinical_data_emc_ordered$Patient_ID_SKY,
                                         clinical_imputed_val$Patient_ID_SKY)])
stopifnot(!any(is.na(clinical_data_emc_ordered$BWH)))

clinical_data_emc_ordered$Sample_type <- ifelse(clinical_data_emc_ordered$Sample_type == "Biopsy","Biopsy","Excision")
clinical_data_emc_ordered$rank <- c(1:nrow(clinical_data_emc_ordered))
col_ha = HeatmapAnnotation(Metastasis = clinical_data_emc_ordered$Metastasis, 
                           BWH = clinical_data_emc_ordered$BWH,
                           Biopsy_excision = clinical_data_emc_ordered$Sample_type,
                           rank = clinical_data_emc_ordered$rank,
                           col=c(cols_clinical_data,rank = circlize::colorRamp2(c(0,ncol(all_models_scaled_ordered)),colors = c("white","#0e7983"))))

# Order rows (genes) by coefficient
coeff_final_model_ordered <- coeff_final_model[order(coeff_final_model)]
all_models_scaled_ordered <- all_models_scaled_ordered[match(names(coeff_final_model_ordered),rownames(all_models_scaled_ordered)),]
colnames(all_models_scaled_ordered) <- NULL

pdf(file = paste0(output_dir, "/Supp_figure_13A_heatmap_logTPMs_scaled_EMC_validation_23genes_model_Romy_ungrouped_prediction_genes_rank.pdf"), width = 9, height = 5)
draw(Heatmap(all_models_scaled_ordered, top_annotation = col_ha, name = "scaled",
             cluster_columns = FALSE, cluster_rows = FALSE))
dev.off()

# Source data
df_val <- as.data.frame(t(all_models_scaled_ordered))
df_val$rank <- clinical_data_emc_ordered$rank
df_val$Biopsy_excision <- clinical_data_emc_ordered$Sample_type
df_val$BWH <- clinical_data_emc_ordered$BWH
df_val$Metastasis <- clinical_data_emc_ordered$Metastasis
#df_val$Sample_id <- gsub("[AB]", "", df_val$Sample_id)
df_val$Sample_id <- gsub("[AB]$", "", clinical_data_emc_ordered$Skyline_ID)

writexl::write_xlsx(df_val,paste0(output_dir,"/Supp_figure_13A_heatmap_validationN102_source_data.xlsx"))

#-------------------------------------------------------------------------------

