#-------------------------------------------------------------------------------
# Aim: Perform recalibration of SCCore-GEP predictions on D-SQUAME validation dataset
# Author: L.Pozza, recalibration functions were coded by B. Rentroia-Pacheco
# Input: Uncalibrated SCCore-GEP predictions
# Output: Cox recalibrated SCCore-GEP predictions
#-------------------------------------------------------------------------------

# Functions
#-------------------------------------------------------------------------------
# recalib_preds function: This function recalibrates the predictions using Cox method
# @param predictions: numeric vector with predictions
# @param y: dataset with outcomes of interest and matching information
# @param outcome_n: character corresponding to the column name of the outcome of interest in the y dataframe
# @param tp: number indicating time point
# @param weights_n: character corresponding to the column name of weights in the y dataframe
recalib_preds <- function(predictions, y, outcome_n, fup_n, tp, weights_n){
  
  # 1. Prepare data
  # Avoid having 0 or 1 predictions
  predictions <- ifelse(1 - predictions < 1e-16, 1 - 1e-16, predictions)
  predictions <- ifelse(predictions < 1e-16, 1e-16, predictions)
  
  # 2. Recalibrate predictions
  loglog_predictions <- log(-log(predictions))
  recalib_model <- coxph(Surv(y[, fup_n], y[, outcome_n])~loglog_predictions, weights = y[, weights_n])
  lp <- predict(recalib_model, type = "lp")
  baseline_surv <- summary(survfit(recalib_model), times = tp)$surv
  recalib_predictions <- 1-baseline_surv^exp(lp)
  recalib_predictions <- as.vector(recalib_predictions)
  
  # 3. Output results
  output <- list(uncalibrated_predictions = predictions,
                 recalibrated_predictions = recalib_predictions,
                 recalibration_cox_lp = lp,
                 recalibration_model = recalib_model)
  return(output)
}
#------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(openxlsx)
library(tidyverse)
library(survival)
library(glmnet)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Clinical data
clinical_fn <- e3_clinical
# Weights
lr_weights_fn <- e7_weights
# SCCore-GEP predictions
sccoregep_preds_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "SCCoreGEP_predictions.csv")
# Time point
tp <- 5
# Output folder
output_dir <- file.path(results_dir, "intermediate", "validation", "d_squame")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
# Clinical data
clinical <- read.csv(clinical_fn)
# Weights
lr_weights <- read.csv(lr_weights_fn)
# SCCore-GEP predictions
sccoregep_preds <- read.csv(sccoregep_preds_fn)
#-------------------------------------------------------------------------------

# Preprocess data
#-------------------------------------------------------------------------------
# Keep only biopsy for patients with double sample type, add column with weights and with predictions
clinical <- clinical %>%
  filter(!(Set_info == "complete_set_double_sampletype" & Type_of_material_bin != "Biopsy")) %>%
  mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0)) %>%
  left_join(lr_weights %>% select(Patient_ID_SKY, Weight_rescaled) %>% filter(!is.na(Weight_rescaled)) %>% unique(), by = "Patient_ID_SKY") %>%
  left_join(sccoregep_preds, by = c("Patient_ID_SKY", "Skyline_ID")) %>%
  select(Skyline_ID, Patient_ID_SKY, Metastasis_num, FU_metastasis_years, Weight_rescaled, Type_of_material_bin, sccore_gep)
#-------------------------------------------------------------------------------

# Perform Cox recalibration
#-------------------------------------------------------------------------------
outcome_var <- "Metastasis_num"
time_var <- "FU_metastasis_years"
weights_var <- "Weight_rescaled"
recalib_cox_res <- recalib_preds(clinical$sccore_gep, clinical, outcome_var, time_var, tp, weights_var)
recalib_preds <- data.frame(Skyline_ID = clinical$Skyline_ID,
                            Patient_ID_SKY = clinical$Patient_ID_SKY,
                            sccore_gep_uncalibrated = clinical$sccore_gep,
                            sccore_gep_cox_calibrated = recalib_cox_res$recalibrated_predictions)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
# Dataframe with recalibrated predictions
write.csv(recalib_preds, file.path(output_dir, "recalibrated_SCCoreGEP_predictions.csv"), row.names = F)
#-----------------------
