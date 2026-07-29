#-------------------------------------------------------------------------------
# Aim: Apply SCCore-GEP in D-SQUAME validation dataset
# Author: L.Pozza, adapted from J. Traets
# Input: TPM & imputed clinical data of D-SQUAME validation dataset + SCCore-GEP model
# Output: SCCore-GEP predictions
#-------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(openxlsx)
library(tidyverse)
library(survival)
library(glmnet)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(pROC::roc)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Clinical data
clinical_fn <- e3_clinical
# Weights
lr_weights_fn <- e7_weights
# Expression data
tpms_fn <- e2_tpms
# SCCore-GEP
sccoregep_fn <- file.path(results_dir, "intermediate", "sccore_gep","Discovery_models",
                          "coxnet_Late_feature_integration_robj.rds")
# Time point at which c-index is computed
tp <- 5
# Output folder
output_dir <- file.path(results_dir, "intermediate", "validation", "d_squame")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
# Clinical data
clinical <- read.csv(clinical_fn)
## Keep only biopsy for patients with duplicated sample
clinical <- clinical %>%
      filter(!(Set_info == "complete_set_double_sampletype" & Type_of_material_bin != "Biopsy"))
# Weights
lr_weights <- read.csv(lr_weights_fn)
# Gene expression data
tpms <- read.table(tpms_fn, sep = "\t", header = T)
# SCCore-GEP model
sccoregep <- readRDS(sccoregep_fn)
#-------------------------------------------------------------------------------

# Preprocess gene expression data
#-------------------------------------------------------------------------------
# Fix sample names to match expression data
add_x_if_number <- function(vec) {ifelse(grepl("^[0-9]", vec), paste0("X", vec), vec)}
clinical$tpm_sample_id <- add_x_if_number(clinical$Skyline_ID)
clinical$tpm_sample_id <- gsub("-",".",clinical$tpm_sample_id)

# Missing 7 samples, those have been repeated so they have a "r" in their name
# "24-0330" "24-0521" "24-0609" "24-0610" "S9624"   "S9625"   "S9623"  
diff_samples <- base::setdiff(clinical$tpm_sample_id, colnames(tpms))

# Generate new column with correct names for the tpm/count data
new_sample_names <- unlist(lapply(diff_samples, function(sample_x) {colnames(tpms)[grepl(sample_x, colnames(tpms))]}))
new_sample_names <- new_sample_names[!grepl("_1$|_2$|_2_val_1$|1_val_1$|_1_val_2$|_2_val_2$",new_sample_names)]
names(new_sample_names) <- diff_samples
clinical[match(names(new_sample_names),clinical$tpm_sample_id),]$tpm_sample_id <- new_sample_names

# Transform TPM into log2TPM
rownames(tpms) <- tpms$gene_id
tpms_filtered <- tpms %>% select(clinical$tpm_sample_id) %>% t() %>% as.data.frame()
tpms_filtered_log <- log2(tpms_filtered + 1)
#-------------------------------------------------------------------------------

# Apply SCCore-GEP model
#-------------------------------------------------------------------------------
# Prepare data
## Removing .# behind gene ids, not in feats_to_keep/x_fitted_model?
newx_data <- tpms_filtered_log
colnames(newx_data) <- sub("\\..*", "", colnames(newx_data))
final_fitted_model_sccoregep <- sccoregep$Final_model2$fitted_model
feats_to_keep_sccoregep <- colnames(sccoregep$Final_model2$x_fitted_model)
newx_data_sccoregep <- as.matrix(newx_data[,feats_to_keep_sccoregep])
# Apply model
sccore_gep_preds <- 1 - as.vector(summary(survfit(final_fitted_model_sccoregep, 
                                     s = "lambda.1se", 
                                     x = sccoregep$Final_model2$x_fitted_model, 
                                     y = sccoregep$Final_model2$y_fitted_model, 
                                     weights = sccoregep$Final_model2$w_fitted_model, 
                                     newx = as.matrix(newx_data_sccoregep)), times = tp)$surv)
names(sccore_gep_preds) <- rownames(newx_data)
# Save predictions in dataframe
model_predictions_df <- as.data.frame(sccore_gep_preds) %>%
        rownames_to_column(var = "tpm_sample_id") %>%
        `colnames<-`(gsub("_preds", "", colnames(.))) %>%
        full_join(clinical %>%
                     select(tpm_sample_id, Skyline_ID, Patient_ID_SKY),
                  by = "tpm_sample_id") %>%
        select(Skyline_ID, Patient_ID_SKY, sccore_gep)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(model_predictions_df, file.path(output_dir, "SCCoreGEP_predictions.csv"), row.names = F)
#-------------------------------------------------------------------------------