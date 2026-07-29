#-------------------------------------------------------------------------------
# Aim: Generate supp figure 12a
# Input: validation sequencing data objects
# Output: heatmap
#-------------------------------------------------------------------------------


# Load libraries
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


# General setup
#-------------------------------------------------------------------------------
# output folder
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------


# Final model
# TODO replace with final model in results folder!
#-------------------------------------------------------------------------------
GEP_final_model <- readRDS(m1_final_model)
final_model_genes <- names(GEP_final_model$Final_model2$coefficients)
stopifnot(length(final_model_genes)==23)
#-------------------------------------------------------------------------------


# Predictions validation EMC
# TODO replace with predictions in results folder!
#-------------------------------------------------------------------------------
p_pred_val_emc <- file.path(output_dir,"emc", "table_5_validation_performance_predictions.csv")
predictions_validation_emc <- read.csv(p_pred_val_emc)
#-------------------------------------------------------------------------------


# Predictions validation Nassir et al
# TODO replace with predictions in results folder!
#-------------------------------------------------------------------------------
p_pred_val_nassir <- file.path(output_dir,"nassir", "table_5_validation_performance_predictions.csv")
predictions_validation_nassir <- read.csv(p_pred_val_nassir)
#-------------------------------------------------------------------------------


# Heatmap plot
#-------------------------------------------------------------------------------
# Load validation data
validation_data <- readRDS(paste0(output_data_dir,"/Validation_data_set.rds"))
nassir_validation_data <- readRDS(paste0(output_data_dir,"/Nassir_validation_data_set.rds"))
out <- validation_data@out

# Genes and coefficients from the final model
coeff_final_model <- GEP_final_model$Final_model2$coefficients
final_model_ensg <- out[match(names(coeff_final_model),do.call(rbind,strsplit(out$gene_id,"[.]"))[,1]),]$gene_id
names(coeff_final_model) <- out[match(names(coeff_final_model),do.call(rbind,strsplit(out$gene_id,"[.]"))[,1]),]$gene_name

# Heatmap from the Nassir et al data set
tpm_data_nassir_model <- nassir_validation_data@sequencing_data$tpms %>% filter(rownames(nassir_validation_data@sequencing_data$tpms) %in% final_model_ensg)
tpm_data_nassir_model <- log2(tpm_data_nassir_model + 1) #log scale
rownames(tpm_data_nassir_model) <- out[match(rownames(tpm_data_nassir_model),out$gene_id),]$gene_name
nassir_validation_data@clinical$Metastasis <- ifelse(nassir_validation_data@clinical$MetNoMet == "Met","Case","Control")
nassir_validation_data@clinical$BWH <- ifelse(is.na(nassir_validation_data@clinical$BWH),"T2b",nassir_validation_data@clinical$BWH)

# scale data
all_models_scaled <- t(scale(t(tpm_data_nassir_model)))

# order by prediction rank
predictions_nassir_filt <- predictions_validation_nassir[predictions_validation_nassir$Sample_id %in% colnames(all_models_scaled),]
ordered_predictions <- predictions_nassir_filt[order(predictions_nassir_filt$gep23, decreasing = T),]
all_models_scaled_ordered <- all_models_scaled[,ordered_predictions$Sample_id]
clinical_data_nassir_ordered <- nassir_validation_data@clinical[match(ordered_predictions$Sample_id,nassir_validation_data@clinical$Run),]

clinical_data_nassir_ordered$rank <- c(1:nrow(clinical_data_nassir_ordered))
col_ha = HeatmapAnnotation(Metastasis = clinical_data_nassir_ordered$Metastasis, 
                           BWH = clinical_data_nassir_ordered$BWH,
                           rank = clinical_data_nassir_ordered$rank,
                           col=c(cols_clinical_data,rank = circlize::colorRamp2(c(0,ncol(all_models_scaled_ordered)),colors = c("white","#0e7983"))))

# order rows by coefficient
coeff_final_model_ordered <- coeff_final_model[order(coeff_final_model)]
all_models_scaled_ordered <- all_models_scaled_ordered[match(names(coeff_final_model_ordered),rownames(all_models_scaled_ordered)),]
colnames(all_models_scaled_ordered) <- NULL

pdf(file=paste0(output_dir,"/Supp_figure_12a_heatmap_logTPMs_scaled_Nassir_23genes_model_Romy_ungrouped_prediction_genes_rank.pdf"), width=7, height=5)
Heatmap(all_models_scaled_ordered, top_annotation = col_ha, name="scaled", 
        cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

df_nassir <- as.data.frame(t(all_models_scaled_ordered))
df_nassir$rank <- clinical_data_nassir_ordered$rank
df_nassir$BWH <- clinical_data_nassir_ordered$BWH
df_nassir$Metastasis <- clinical_data_nassir_ordered$Metastasis
df_nassir$Sample_id <- clinical_data_nassir_ordered$SampleID

writexl::write_xlsx(df_nassir,paste0(output_dir,"/Supp_figure_12a_heatmap_nassir_source_data.xlsx"))



# Heatmap from the EMC validation data
tpm_data_emc_model <- validation_data@sequencing_data$tpms %>% filter(rownames(validation_data@sequencing_data$tpms) %in% final_model_ensg)
tpm_data_emc_model <- log2(tpm_data_emc_model + 1) #log scale
rownames(tpm_data_emc_model) <- out[match(rownames(tpm_data_emc_model),out$gene_id),]$gene_name

# scale data
all_models_scaled <- t(scale(t(tpm_data_emc_model)))

# fix sample names
predictions_validation_emc$tpm_sample_id <- validation_data@clinical[match(predictions_validation_emc$Skyline_ID,validation_data@clinical$Skyline_ID),]$tpm_sample_id

# order by prediction rank
predictions_emc_filt <- predictions_validation_emc[predictions_validation_emc$tpm_sample_id %in% colnames(all_models_scaled),]
ordered_predictions <- predictions_emc_filt[order(predictions_emc_filt$gep23, decreasing = T),]
all_models_scaled_ordered <- all_models_scaled[,ordered_predictions$tpm_sample_id]
clinical_data_emc_ordered <- validation_data@clinical[match(ordered_predictions$Skyline_ID,validation_data@clinical$Skyline_ID),]

clinical_data_emc_ordered$Sample_type <- ifelse(clinical_data_emc_ordered$Sample_type == "Biopsy","Biopsy","Excision")
clinical_data_emc_ordered$rank <- c(1:nrow(clinical_data_emc_ordered))
col_ha = HeatmapAnnotation(Metastasis = clinical_data_emc_ordered$Metastasis, 
                           BWH = clinical_data_emc_ordered$BWH,
                           Biopsy_excision = clinical_data_emc_ordered$Sample_type,
                           rank = clinical_data_emc_ordered$rank,
                           col=c(cols_clinical_data,rank = circlize::colorRamp2(c(0,ncol(all_models_scaled_ordered)),colors = c("white","#0e7983"))))

# order rows by coefficient
coeff_final_model_ordered <- coeff_final_model[order(coeff_final_model)]
all_models_scaled_ordered <- all_models_scaled_ordered[match(names(coeff_final_model_ordered),rownames(all_models_scaled_ordered)),]
colnames(all_models_scaled_ordered) <- NULL

pdf(file=paste0(output_dir,"/Supp_figure_12a_heatmap_logTPMs_scaled_EMC_validation_23genes_model_Romy_ungrouped_prediction_genes_rank.pdf"), width=9, height=5)
Heatmap(all_models_scaled_ordered, top_annotation = col_ha, name="scaled", 
        cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

df_val <- as.data.frame(t(all_models_scaled_ordered))
df_val$rank <- clinical_data_emc_ordered$rank
df_val$Biopsy_excision <- clinical_data_emc_ordered$Sample_type
df_val$BWH <- clinical_data_emc_ordered$BWH
df_val$Metastasis <- clinical_data_emc_ordered$Metastasis
df_val$Sample_id <- gsub("[AB]", "", df_val$Sample_id)

writexl::write_xlsx(df_val,paste0(output_dir,"/Supp_figure_12a_heatmap_validationN102_source_data.xlsx"))

#-------------------------------------------------------------------------------
