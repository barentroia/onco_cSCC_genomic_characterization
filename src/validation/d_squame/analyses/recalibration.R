#===============================================================================
# Aim:      Recalibrate the SCCore-GEP predictions in the D-SQUAME validation
#           dataset with the Cox method (weighted Cox model on
#           log(-log(prediction))), as input for the multivariable analysis.
# Author:   L. Pozza; recalibration function coded by B. Rentroia Pacheco
# Input:    - D-SQUAME validation clinical data (e3_clinical) and NCC weights
#             (e7_weights)
#           - Uncalibrated SCCore-GEP predictions (output of apply_SCCoreGEP.R)
# Output:   - recalibrated_SCCoreGEP_predictions.csv in
#             results_dir/intermediate/validation/d_squame
# Requires: results_dir (main.R) and config.R paths
# Usage:    Sourced from main.R
#===============================================================================

# 1. Functions
#-------------------------------------------------------------------------------
#' Recalibrate predictions with the Cox method
#'
#' Fits a weighted Cox model of the outcome on log(-log(predictions)) and
#' returns the recalibrated event probabilities at tp.
#'
#' @param predictions numeric vector. Predicted event probabilities.
#' @param y data.frame. Outcome, follow-up and weights.
#' @param outcome_n character. Outcome column in y.
#' @param fup_n character. Follow-up time column in y.
#' @param tp numeric. Time point (years).
#' @param weights_n character. Weights column in y.
#'
#' @return list with uncalibrated_predictions, recalibrated_predictions,
#'   recalibration_cox_lp and recalibration_model.
recalib_preds <- function(predictions, y, outcome_n, fup_n, tp, weights_n){
  
  # 1.1 Prepare data
  # Avoid having 0 or 1 predictions
  predictions <- ifelse(1 - predictions < 1e-16, 1 - 1e-16, predictions)
  predictions <- ifelse(predictions < 1e-16, 1e-16, predictions)
  
  # 1.2 Recalibrate predictions
  loglog_predictions <- log(-log(predictions))
  recalib_model <- coxph(Surv(y[, fup_n], y[, outcome_n])~loglog_predictions, weights = y[, weights_n])
  lp <- predict(recalib_model, type = "lp")
  baseline_surv <- summary(survfit(recalib_model), times = tp)$surv
  recalib_predictions <- 1-baseline_surv^exp(lp)
  recalib_predictions <- as.vector(recalib_predictions)
  
  # 1.3 Output results
  output <- list(uncalibrated_predictions = predictions,
                 recalibrated_predictions = recalib_predictions,
                 recalibration_cox_lp = lp,
                 recalibration_model = recalib_model)
  return(output)
}
#-------------------------------------------------------------------------------

# 2. Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(openxlsx)
library(tidyverse)
library(survival)
library(glmnet)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# 3. Filenames and settings
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

# 4. Read data
#-------------------------------------------------------------------------------
# Clinical data
clinical <- read.csv(clinical_fn)
# Weights
lr_weights <- read.csv(lr_weights_fn)
# SCCore-GEP predictions
sccoregep_preds <- read.csv(sccoregep_preds_fn)
#-------------------------------------------------------------------------------

# 5. Preprocess data
#-------------------------------------------------------------------------------
# Keep only biopsy for patients with double sample type, add column with weights and with predictions
clinical <- clinical %>%
  filter(!(Set_info == "complete_set_double_sampletype" & Type_of_material_bin != "Biopsy")) %>%
  mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0)) %>%
  left_join(lr_weights %>% select(Patient_ID_SKY, Weight_rescaled) %>% filter(!is.na(Weight_rescaled)) %>% unique(), by = "Patient_ID_SKY") %>%
  left_join(sccoregep_preds, by = c("Patient_ID_SKY", "Skyline_ID")) %>%
  select(Skyline_ID, Patient_ID_SKY, Metastasis_num, FU_metastasis_years, Weight_rescaled, Type_of_material_bin, sccore_gep)
#-------------------------------------------------------------------------------

# 6. Perform Cox recalibration
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

# 7. Save results
#-------------------------------------------------------------------------------
# Dataframe with recalibrated predictions
write.csv(recalib_preds, file.path(output_dir, "recalibrated_SCCoreGEP_predictions.csv"), row.names = F)
#-------------------------------------------------------------------------------
