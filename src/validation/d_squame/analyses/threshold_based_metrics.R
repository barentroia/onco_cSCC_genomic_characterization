#===============================================================================
# Aim:      Compute threshold-based metrics (weighted SE, SP, PPV, NPV, ...) of
#           the recalibrated SCCore-GEP predictions in the D-SQUAME validation
#           dataset, at every observed risk threshold, overall and per BWH and
#           AJCC8 stage group.
# Author:   L. Pozza
# Input:    - D-SQUAME validation clinical data (e3_clinical), merged imputed
#             clinical data and NCC weights (e7_weights)
#           - Recalibrated SCCore-GEP predictions (output of recalibration.R)
# Output:   - threshold_based_metrics.csv in
#             results_dir/intermediate/validation/d_squame
# Requires: results_dir (main.R) and config.R paths
# Usage:    Sourced from main.R
#===============================================================================

# 1. Functions
#-------------------------------------------------------------------------------
#' (Weighted) threshold-based metrics for a binary classification
#'
#' Samples with prediction > cutoff are classified as high risk.
#'
#' @param cutoff numeric. Threshold on the predictions.
#' @param outcome_num numeric vector. Outcome (1 = event, 0 = no event).
#' @param predictions numeric vector. Predicted risks.
#' @param weights numeric vector or NULL. Sample weights (NULL = unweighted).
#'
#' @return Named numeric vector: weighted FP, FN, TP, TN, SE, SP, PPV, NPV,
#'   likelihood ratios, proportions/numbers at low and high risk, positives
#'   and negatives.
computeDiscMetricsBinary <- function(cutoff, outcome_num, predictions, weights = NULL){
  
  # Obtain labels for the desired cutoff
  labels_cutoff <- ifelse(predictions > cutoff, 1, 0)
  
  # Weights are all equal to 1 if they are not specified
  if (is.null(weights)){
    weights <- rep(1, length(predictions))
  }
  
  # Compute TP, TN, FP and FN 
  correctly_assigned <- which(outcome_num == labels_cutoff)
  wrongly_assigned <-  which(outcome_num != labels_cutoff)
  positives <- which(outcome_num == 1)
  negatives <- which(outcome_num == 0)
  ## TN
  TN <- sum(weights[intersect(correctly_assigned, negatives)])
  ## TP
  TP <- sum(weights[intersect(correctly_assigned, positives)])
  ## FN
  FN <- sum(weights[intersect(wrongly_assigned, positives)])
  ## FP
  FP <- sum(weights[intersect(wrongly_assigned, negatives)])
  ## Positives
  Pos <- sum(weights[positives])
  ## Negatives
  Neg <- sum(weights[negatives])
  
  # Compute Sensitivity
  if ((TP + FN) > 0){
    SE <- TP / (TP + FN)
  } else {
    SE <- NA
  }
  # Compute Specificity
  if ((TN + FP) > 0){
    SP <- TN / (TN + FP)
  } else {
    SP <- NA
  }
  # Compute PPV
  if ((TP + FP) > 0){
    PPV <- TP / (TP + FP)
  } else {
    PPV <- NA
  }
  # Compute NPV
  if ((TN + FN) > 0){
    NPV <- TN / (TN + FN)
  } else {
    NPV <- NA
  }
  # Compute positive and negative likelihood ratio
  if (!is.na(SE) & !is.na(SP)){
    LR_pos <- SE /(1 - SP)
    LR_neg <- (1 - SE) / SP
  }
  # Proportion of HR and LR patients
  low_risk_n <- FN + TN
  high_risk_n <- FP + TP
  low_risk_prop <- low_risk_n / (low_risk_n + high_risk_n)
  high_risk_prop <- 1 - low_risk_prop
  
  # Summarize all metrics
  performance_metrics <-c("FP_w" = FP,
                          "FN_w" = FN,
                          "TP_w" = TP,
                          "TN_w" = TN,
                          "SE_w" = SE,
                          "SP_w" = SP,
                          "PPV_w" = PPV,
                          "NPV_w" = NPV,
                          "LR_pos_w" = LR_pos,
                          "LR_neg_w" = LR_neg,
                          "low_risk_w" = low_risk_prop,
                          "high_risk_w" = high_risk_prop,
                          "positives_w" = Pos,
                          "negatives_w" = Neg,
                          "n_low_risk_w" = low_risk_n,
                          "n_high_risk_w" = high_risk_n)
  
  return(performance_metrics)
}

#' Threshold-based metrics per stage group
#'
#' Computes the metrics of computeDiscMetricsBinary() at each threshold
#' within each stratum of BWH and AJCC8 (and their groupings T1-T2a,
#' T1, T1-T2, ...).
#'
#' @param df data.frame. Merged imputed clinical data with BWH_staging and
#'   AJCC_staging.
#' @param weights_df data.frame. Outcome, weights and predictions per patient.
#' @param ids_var character. Patient ID column used to join the data.
#' @param status_var character. Outcome column.
#' @param time_var character. Follow-up time column.
#' @param weights_var character. Weights column.
#' @param preds_var character. Predictions column.
#' @param tp numeric. Time point (not used in the computation).
#' @param thresholds numeric vector. Thresholds at which metrics are computed.
#'
#' @return Long data.frame with stratification_system, stratum, metric,
#'   value and threshold; only thresholds that occur in the stratum are kept.
compute_thresh_based_metrics_strata <- function(df, weights_df, ids_var, status_var, time_var, weights_var, preds_var, tp, thresholds){
  df_gathered <- df %>%
    full_join(weights_df %>% select(-all_of(time_var)), by = ids_var) %>%
    mutate(BWH_staging_T1T2a_T2bT3 = ifelse(BWH_staging %in% c("T1", "T2a"), "T1-T2a", "T2b-T3"),
           BWH_staging_T1_T2aT2bT3 = ifelse(BWH_staging %in% c("T1"), "T1", "T2a-T2b-T3"),
           AJCC_staging_T1T2_T3T4 = ifelse(AJCC_staging %in% c("T1", "T2"), "T1-T2", "T3-T4"),
           AJCC_staging_T1_T2T3T4 = ifelse(AJCC_staging %in% c("T1"), "T1", "T2-T3-T4")) %>%
    select(all_of(c(ids_var, status_var, time_var, weights_var, preds_var,
                    "BWH_staging", "BWH_staging_T1T2a_T2bT3", "BWH_staging_T1_T2aT2bT3",
                    "AJCC_staging", "AJCC_staging_T1T2_T3T4", "AJCC_staging_T1_T2T3T4"))) %>%
    gather(stratification_system, stratum, -all_of(c(ids_var, status_var, time_var, weights_var, preds_var)))
  
  res <- vector("list", length(thresholds))
  for (i in 1:length(thresholds)){
    threshold <- thresholds[[i]]
    res_i <- df_gathered %>%
      group_by(stratification_system, stratum) %>%
      summarise(metrics = list(computeDiscMetricsBinary(threshold, .data[[status_var]], .data[[preds_var]], .data[[weights_var]])), .groups = "keep") %>%
      unnest_wider(metrics) %>%
      ungroup() %>%
      gather(metric, value, -stratification_system, -stratum)
    res[[i]] <- res_i %>%
      mutate(threshold = threshold)
  }
  res <- res %>% bind_rows() %>%
    inner_join(df_gathered %>%
                 select(all_of(c("stratification_system", "stratum", preds_var))) %>%
                 `colnames<-`(c("stratification_system", "stratum", "threshold")),
               by = c("stratification_system", "stratum", "threshold"))
  res
}
#-------------------------------------------------------------------------------

# 2. Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(openxlsx)
library(tidyverse)
library(survival)
library(mice)
library(boot)
library(spatstat)
library(prodlim)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# 3. Filenames and settings
#-------------------------------------------------------------------------------
# Clinical data
clinical_fn <- e3_clinical
clinical_mimp_merged_fn <- file.path(results_dir, "intermediate", "clinical", "merged_clinical_data_imputed_d_squame_validation.csv")
# Weights
lr_weights_fn <- e7_weights
# Recalibrated SCCore-GEP predictions
sccoregep_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "recalibrated_SCCoreGEP_predictions.csv")
# Time point
tp <- 5
# Output folder
output_dir <- file.path(results_dir, "intermediate", "validation", "d_squame")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# 4. Read data
#-------------------------------------------------------------------------------
# Clinical data
clinical <- read.csv(clinical_fn)
clinical_mimp_merged <- read.csv(clinical_mimp_merged_fn)
# Weights
lr_weights <- read.csv(lr_weights_fn)
# SCCore-GEP recalibrated predictions
sccoregep <- read.csv(sccoregep_fn)
#-------------------------------------------------------------------------------

# 5. Preprocess clinical and predictions data
#-------------------------------------------------------------------------------
# Keep only biopsy for patients with double sample type, add column with weights
df <- clinical %>%
  filter(!(Set_info == "complete_set_double_sampletype" & Type_of_material_bin != "Biopsy")) %>%
  mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0)) %>%
  left_join(lr_weights %>% select(Patient_ID_SKY, Weight_rescaled) %>% unique(), by = "Patient_ID_SKY") %>%
  select(Skyline_ID, Patient_ID_SKY, Metastasis_num, FU_metastasis_years, Weight_rescaled) %>%
  left_join(sccoregep %>% select(Skyline_ID, sccore_gep_cox_calibrated), by = "Skyline_ID")
#-------------------------------------------------------------------------------

# 6. Compute threshold-based metrics at different thresholds
#-------------------------------------------------------------------------------
# Define variables
outcome_var <- "Metastasis_num"
time_var <- "FU_metastasis_years"
weights_var <- "Weight_rescaled"
preds_var <- "sccore_gep_cox_calibrated"
thresholds <- sort(unique(df[[preds_var]]))

# Compute threshold-based metrics in entire dataset
thresh_based_metrics_all <- sapply(thresholds, computeDiscMetricsBinary, df[[outcome_var]], df[[preds_var]], df[[weights_var]], simplify = F) %>%
  bind_rows() %>%
  cbind(data.frame(threshold = thresholds))

# Compute threshold-based metrics in groups
thresh_based_metrics_mimp_merged_res <- compute_thresh_based_metrics_strata(clinical_mimp_merged,
                                                                            df,
                                                                            "Patient_ID_SKY",
                                                                            outcome_var,
                                                                            time_var,
                                                                            weights_var,
                                                                            preds_var,
                                                                            tp,
                                                                            thresholds)
thresh_based_metrics_mimp_merged <- thresh_based_metrics_mimp_merged_res %>%
  spread(metric, value)

# Combine results
thresh_based_metrics_mimp_merged_comb <- thresh_based_metrics_all %>%
  mutate(stratification_system = "All",
         stratum = "All") %>%
  rbind(thresh_based_metrics_mimp_merged %>%
          select(all_of(c("stratification_system", "stratum", colnames(thresh_based_metrics_all))))) %>%
  relocate(stratification_system, stratum, threshold) %>%
  arrange(stratification_system, stratum, threshold)
#-------------------------------------------------------------------------------

# 7. Save results
#-------------------------------------------------------------------------------
# Dataframe with results
write.csv(thresh_based_metrics_mimp_merged_comb, file.path(output_dir, "threshold_based_metrics.csv"), row.names = F)
#-------------------------------------------------------------------------------
