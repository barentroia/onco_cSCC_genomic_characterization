#-------------------------------------------------------------------------------
# Aim: Generate Figure 5a heatmap discovery
# Input: discovery sequencing data object
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
#-------------------------------------------------------------------------------
sccoregep_fn <- file.path(results_dir, "intermediate", "sccore_gep", "Discovery_models", "coxnet_Late_feature_integration_robj.rds") 
GEP_final_model <- readRDS(sccoregep_fn)
final_model_genes <- names(GEP_final_model$Final_model2$coefficients)
stopifnot(length(final_model_genes)==23)

#-------------------------------------------------------------------------------


# Heatmap plot
#-------------------------------------------------------------------------------
# Load validation data
discovery_data <- readRDS(paste0(output_data_dir,"/Discovery_data_set_N366.rds"))
out <- discovery_data@out

# Genes and coefficients from the final model
coeff_final_model <- GEP_final_model$Final_model2$coefficients
final_model_ensg <- out[match(names(coeff_final_model),do.call(rbind,strsplit(out$gene_id,"[.]"))[,1]),]$gene_id
names(coeff_final_model) <- out[match(names(coeff_final_model),do.call(rbind,strsplit(out$gene_id,"[.]"))[,1]),]$gene_name

# Select genes from final model
tpm_data_discovery_model <- discovery_data@sequencing_data$tpms %>% filter(rownames(discovery_data@sequencing_data$tpms) %in% final_model_ensg)
tpm_data_discovery_model <- log2(tpm_data_discovery_model + 1) #log scale
rownames(tpm_data_discovery_model) <- out[match(rownames(tpm_data_discovery_model),out$gene_id),]$gene_name
all_models_scaled <- t(scale(t(tpm_data_discovery_model)))

# Predictions of the final model
predictions_final_model_disc <- GEP_final_model$Final_model2$predictions[colnames(all_models_scaled)]

# order by prediction score
ordered_predictions <- predictions_final_model_disc[order(predictions_final_model_disc, decreasing = T)]
all_models_scaled_ordered <- all_models_scaled[,names(ordered_predictions)]
clinical_data_disc_ordered <- discovery_data@clinical[match(names(ordered_predictions),discovery_data@clinical$SkylineDx.ID),]

# order rows by coefficient
coeff_final_model_ordered <- coeff_final_model[order(coeff_final_model)]
all_models_scaled_ordered <- all_models_scaled_ordered[match(names(coeff_final_model_ordered),rownames(all_models_scaled_ordered)),]

clinical_data_disc_ordered$Biopsy_excision <- ifelse(clinical_data_disc_ordered$Biopsy_excision == "Biopsy","Biopsy","Excision")
clinical_data_disc_ordered$rank <- c(1:nrow(clinical_data_disc_ordered))
col_ha <- HeatmapAnnotation(Metastasis = clinical_data_disc_ordered$Metastasis, 
                            BWH = clinical_data_disc_ordered$BWH,
                            Biopsy_excision = clinical_data_disc_ordered$Biopsy_excision,
                            rank = clinical_data_disc_ordered$rank,
                            col=c(cols_clinical_data,rank = circlize::colorRamp2(c(0,ncol(all_models_scaled_ordered)),colors = c("white","#0e7983"))))


colnames(all_models_scaled_ordered) <- NULL
pdf(file=paste0(output_dir,"/Figure_5a_heatmap_logTPMs_scaled_Discovery_23genes_model_Romy_ungrouped_prediction_genes_annotation.pdf"), width=10, height=5)
Heatmap(all_models_scaled_ordered, top_annotation = col_ha, name="scaled", 
        cluster_columns = FALSE, cluster_rows = FALSE)
dev.off()

# save_disc <- list(final_model_scaled=all_models_scaled_ordered,
#                   clinical_data_disc=clinical_data_disc_ordered %>% select(SkylineDx.ID, Metastasis, BWH, Biopsy_excision),
#                   coeff_final_model=coeff_final_model,
#                   predictions_final_model=ordered_predictions)
# saveRDS(save_disc,paste0(output_dir,"/Raw_data_final_model_discovery_heatmap.rds"))

df_disc <- as.data.frame(t(all_models_scaled_ordered))
df_disc$rank <- clinical_data_disc_ordered$rank
df_disc$Biopsy_excision <- clinical_data_disc_ordered$Biopsy_excision
df_disc$BWH <- clinical_data_disc_ordered$BWH
df_disc$Metastasis <- clinical_data_disc_ordered$Metastasis
df_disc$Sample_id <- clinical_data_disc_ordered$Sample_id
df_disc$Sample_id <- gsub("[AB]", "", df_disc$Sample_id)
df_disc$Sample_id <- paste0("d",df_disc$Sample_id)

writexl::write_xlsx(df_disc,paste0(output_dir,"/Figure_5a_heatmap_discoveryN366_source_data.xlsx"))


#-------------------------------------------------------------------------------



